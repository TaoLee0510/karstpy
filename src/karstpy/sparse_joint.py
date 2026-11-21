"""Fit and predict sparse joint karyotype fitness landscapes."""

from dataclasses import dataclass, replace
import math
from time import perf_counter
import warnings

import numpy as np
import pandas as pd

from . import _core
from ._problem import (
    KarstDomainCoverageError, build_problem, design_matrix, domain_and_coverage, major_distance,
)
from ._validate import integer, scalar, validate_beta
from .convergence import optimise_scales, warn_if_unconverged
from .data import KarstData, karst_data_times, prepare_karst_data
from .karyotypes import format_karyotypes, parse_karyotypes
from .options import get_option
from .rng import normal_array
from .tmb import KarstInnerFailure


_DEFAULT_CALIBRATION = {"floor": 0.009, "deficit": 0.6, "deficit_unobserved": 1.1}
_FITNESS_COLUMNS = ["fitness_mean", "fitness_sd", "conf_low", "conf_high",
                    "relative_fitness", "relative_sd", "relative_sd_model", "shell"]


class KarstSamplesLeftOutWarning(UserWarning):
    """One or more input time points are outside the part of the domain they inform."""


@dataclass
class KarstFit:
    """Fitted landscape, diagnostics, and the stored latent posterior mode."""

    data: KarstData
    summary: pd.DataFrame
    dropped: list[str]
    settings: dict
    posterior: dict
    trend: dict
    diagnostics: dict
    input: KarstData | None = None

    def __repr__(self):
        return (
            f"KarstFit({self.diagnostics['n_states']} states, {self.diagnostics['n_major']} major karyotypes, "
            f"{self.data.counts.shape[1]} time points, {self.diagnostics['message']})"
        )

    def predict(self, karyotypes):
        """Predict relative fitness and uncertainty inside or outside the fitted domain."""
        labels = _prediction_labels(karyotypes)
        out = self.summary.set_index("karyotype")[_FITNESS_COLUMNS].reindex(labels)
        inside = out.index.isin(self.summary.karyotype)
        if not inside.all():
            outside = np.flatnonzero(~inside)
            predicted = _predict_outside(self, out.index[outside])
            out.iloc[outside] = predicted[_FITNESS_COLUMNS].to_numpy()
        out["in_domain"] = inside
        return out.rename_axis("karyotype").reset_index()


@dataclass
class Posterior:
    objective: object
    point: np.ndarray
    expected: np.ndarray
    factor: object
    value: float
    restart: bool
    stored: bool

    @property
    def mode(self):
        return self.objective.unpack(self.point)

    def solve(self, rhs):
        return self.factor.solve(rhs).reshape(np.shape(rhs))

    def draw(self, size, stream):
        return self.factor.draw(normal_array(stream, (self.objective.random.size, size)))


def _prediction_labels(karyotypes):
    values = np.asarray(karyotypes)
    if values.ndim == 0 and isinstance(values.item(), str):
        values = values.reshape(1)
    if values.ndim == 1 and all(isinstance(item, str) for item in values):
        return [] if values.size == 0 else format_karyotypes(parse_karyotypes(values.tolist(), unique=False))
    return format_karyotypes(values)


def _calibration(calibration):
    if calibration is None:
        return dict.fromkeys(_DEFAULT_CALIBRATION, 0.0)
    if not isinstance(calibration, dict) or not {"floor", "deficit"} <= calibration.keys() <= _DEFAULT_CALIBRATION.keys():
        raise ValueError("`calibration` needs floor and deficit with optional deficit_unobserved.")
    result = {name: scalar(value, name, 0) for name, value in calibration.items()}
    result.setdefault("deficit_unobserved", result["deficit"])
    return result


def calibrate(sd, relative, shell, birth_rate, calibration, observed):
    deficit = np.where(observed, calibration["deficit"], calibration.get("deficit_unobserved", calibration["deficit"]))
    return np.sqrt(sd**2 + (deficit * relative)**2 + (calibration["floor"] * birth_rate * (1 + shell))**2)


def self_probability(karyotypes, beta):
    stay = {cn: math.log(_core.pij(int(cn), int(cn), beta)) for cn in np.unique(karyotypes)}
    return np.exp([sum(stay[cn] for cn in row) for row in karyotypes])


def _fitness_columns(self_factor, log_birth, sd):
    critical = 1.9599639845400536
    fitness = self_factor * np.exp(log_birth)
    return dict(
        log_birth_mean=log_birth, log_birth_sd=sd, fitness_mean=fitness, fitness_sd=fitness * sd,
        conf_low=self_factor * np.exp(log_birth - critical * sd),
        conf_high=self_factor * np.exp(log_birth + critical * sd),
    )


def _predict_outside(fit, labels):
    trend = fit.trend
    cn = parse_karyotypes(labels, unique=False)
    if cn.shape[1] != trend["centre"].size:
        raise ValueError("Karyotypes have a different chromosome dimension from the fit.")
    design, _, shell = design_matrix(cn, trend["centre"], trend["major"], trend["shell_power"])
    birth = design @ trend["coefficients"]
    sd = np.sqrt(np.maximum(((design @ trend["covariance"]) * design).sum(axis=1), 0) + trend["scales"]["residual"]**2)
    reference = trend["reference"]
    covariance = design @ reference["coef_covariance"]
    observed = np.isin(labels, fit.dropped) & (trend["observed_level"] is not None)
    if observed.any():
        birth[observed] = trend["observed_level"]["log_birth"]
        sd[observed] = trend["observed_level"]["sd"]
        covariance[observed] = 0
    values = _fitness_columns(2 * self_probability(cn, fit.data.beta) - 1, birth, sd)
    fitness = values["fitness_mean"]
    relative = fitness - reference["fitness"]
    model_sd = np.sqrt(np.maximum(fitness**2 * sd**2 + reference["variance"] - 2 * fitness * covariance, 0))
    values.update(relative_fitness=relative, relative_sd_model=model_sd, shell=shell,
                  relative_sd=calibrate(model_sd, relative, shell, trend["birth_rate"], trend["calibration"], observed))
    return pd.DataFrame(values)


def _expected(problem, objective, point, data):
    report = objective.report(point)
    mode = objective.unpack(point)
    expected = dict(data, balance_offset=np.asarray(report["balance_residual"]))
    if problem.model == "sparse_joint":
        expected["y0"], expected["y1"] = np.exp(mode["zeta"]), np.exp(report["log_intensity1"])
        counts = np.column_stack([expected["y0"], expected["y1"]])
    elif problem.model == "sparse_joint_series":
        counts = np.column_stack([np.exp(mode["zeta"]), np.exp(report["log_intensity"])])
        expected["y"] = counts
    elif problem.model == "sparse_joint_pairwise":
        member = problem.domain.member
        start = np.exp(mode["zeta"]) * member
        end = np.exp(report["log_intensity"]) * member
        expected.update(y_start=start, y_end=end)
        counts = np.column_stack([np.where(member[:, 0], start[:, 0], np.nan), np.where(member, end, np.nan)])
    else:
        counts = np.column_stack([np.exp(mode["zeta"]), np.exp(mode["zeta_later"])])
        expected.update(y=counts, transition_offset=report["transition_residual"])
    return expected, counts


def _posterior(problem, objective, settings, fixed, start=None, stored=None, stored_value=np.nan, tmb_data=None):
    data = problem.tmb_data if tmb_data is None else tmb_data
    restart = False
    if stored is None:
        try:
            found = objective.find_mode(fixed, start)
        except KarstInnerFailure:
            found = objective.find_mode(fixed, objective.initial[objective.random])
            restart = True
        point, value = found.point, found.laplace
    else:
        point, value = stored, stored_value
    expected, counts = _expected(problem, objective, point, data)
    information = problem.objective(settings, expected)
    hessian = information.hessian(point)
    if not information.factor.factorize(hessian):
        raise KarstInnerFailure("Expected information is not positive definite.")
    return Posterior(objective, point, counts, information.factor, value, restart, stored is not None)


def rebuild_posterior(fit, problem=None, tmb_data=None):
    """Rebuild at the stored mode, or refit the latent variables for modified data."""
    if problem is None:
        problem = build_problem(fit.data, fit.settings)
    objective = problem.objective(fit.settings, tmb_data)
    stored = objective.pack(fit.posterior)
    fixed = stored[objective.fixed]
    if tmb_data is None:
        return _posterior(problem, objective, fit.settings, fixed, stored=stored,
                          stored_value=fit.diagnostics["posterior_objective"])
    return _posterior(problem, objective, fit.settings, fixed, start=stored[objective.random], tmb_data=tmb_data)


def _latent_labels(objective, domain):
    rows = []
    kinds = {"zeta": "initial abundance", "zeta_later": "later abundance", "g": "birth rate",
             "coef": "trend coefficient", "theta": "depth", "theta1": "depth", "balance_level": "balance level"}
    for name, full in objective.slices.items():
        index = full[np.isin(full, objective.random)]
        for i in index - full[0]:
            state = i % len(domain.labels)
            associated = name in ("zeta", "zeta_later", "g")
            rows.append({
                "variable": kinds[name], "karyotype": domain.labels[state] if associated else None,
                "cells": " -> ".join(map(str, domain.y[state])) if associated else None,
                "major": (state < domain.n_major) if associated else None,
            })
    return pd.DataFrame(rows)


def _fit_single(data, settings, control, domain, coverage):
    started = perf_counter()
    problem = build_problem(data, settings, domain)
    objective = problem.objective(settings)
    labels = _latent_labels(objective, domain)
    fixed, value, convergence, message, check = optimise_scales(objective, problem.bounds, control, labels)
    solver = _posterior(problem, objective, settings, fixed)
    posterior = solver.mode
    posterior = {name: posterior[name] for name in [*problem.random, "log_scale"]}
    g_index, coef_index = objective.latent_indices("g"), objective.latent_indices("coef")
    variance = solver.factor.inverse_diagonal()
    sd = np.sqrt(np.maximum(0, variance[g_index]))
    selector = np.zeros((objective.random.size, len(coef_index)))
    selector[coef_index, np.arange(len(coef_index))] = 1
    coef_covariance = solver.solve(selector)[coef_index]
    scales = dict(zip(problem.scale_names, np.exp(posterior["log_scale"])))
    fitness = _fitness_columns(2 * np.exp(problem.kernel.log_self) - 1, posterior["g"], sd)
    major = np.arange(len(domain.labels)) < domain.n_major
    total = problem.used.sum(axis=1)
    factors = (total * major) / total[major].sum()
    contribution = factors * fitness["fitness_mean"]
    rhs = np.zeros(objective.random.size)
    rhs[g_index] = contribution
    covariance = solver.solve(rhs)
    reference = dict(fitness=contribution.sum(), variance=rhs @ covariance, coef_covariance=covariance[coef_index])
    relative = fitness["fitness_mean"] - reference["fitness"]
    relative_sd = np.sqrt(np.maximum(0, fitness["fitness_sd"]**2 + reference["variance"] -
                                    2 * fitness["fitness_mean"] * covariance[g_index]))
    observed = total > 0
    summary = pd.DataFrame({
        "karyotype": domain.labels,
        **{f"count_{name}": col for name, col in zip(data.counts.columns, domain.y.T)},
        **{f"expected_{name}": col for name, col in zip(data.counts.columns, solver.expected.T)},
        "major": major, "shell": problem.shell, **fitness, "relative_fitness": relative,
        "relative_sd": calibrate(relative_sd, relative, problem.shell, settings["birth_rate"], settings["calibration"], observed),
        "relative_sd_model": relative_sd,
    })
    observed_minor = observed & ~major
    observed_level = None
    if observed_minor.any():
        g = posterior["g"][observed_minor]
        observed_level = dict(log_birth=g.mean(), sd=np.sqrt(np.mean(sd[observed_minor]**2) + (g.var(ddof=1) if g.size > 1 else 0)))
    trend = dict(
        coefficients=posterior["coef"], covariance=coef_covariance, centre=problem.centre,
        major=domain.major, shell_power=settings["shell_power"], scales=scales,
        observed_level=observed_level, reference=reference, calibration=settings["calibration"],
        birth_rate=settings["birth_rate"],
    )
    diagnostics = dict(
        convergence=convergence, message=message, objective=value, posterior_objective=solver.value,
        optimisation=check, finite_uncertainty=np.isfinite(relative_sd).all() and np.isfinite(sd).all(),
        n_states=len(domain.labels), n_major=domain.n_major, n_edges=problem.kernel.source.size,
        dropped_karyotypes=len(domain.dropped), dropped_cells=domain.dropped_cells,
        domain_coverage=coverage, birth_rate=settings["birth_rate"],
        fitted_birth_rate=np.exp(posterior["coef"][0]), initial_state=settings["initial_state"],
        history=settings["history"], min_cn=settings["min_cn"], max_cn=settings["max_cn"],
        seconds=perf_counter() - started,
    )
    if "balance_level" in posterior:
        diagnostics["balance_level"] = float(posterior["balance_level"])
        diagnostics["balance_level_sd"] = float(np.sqrt(variance[objective.latent_indices("balance_level")][0]))
    if solver.restart:
        diagnostics["posterior_restart"] = solver.value
    if data.counts.shape[1] > 2:
        diagnostics.update(series=settings["series"], times=karst_data_times(data))
        columns = ["sample", "time", *coverage.columns] if settings["series"] == "anchored" else ["sample", "time", "timepoint", "cells", "uncovered"]
        diagnostics["dropped_samples"] = pd.DataFrame(columns=columns)
        for name in ("pair_history", "intermediate_weight"):
            if settings.get(name) is not None:
                diagnostics[name] = settings[name]
        if "drift" in scales:
            diagnostics["drift"] = scales["drift"]
    return KarstFit(data, summary, domain.dropped, settings.copy(), posterior, trend, diagnostics)


def _coverage_mask(data, settings):
    counts = data.counts
    n_time = counts.shape[1]
    kept = list(range(n_time))
    dropped = []
    times = karst_data_times(data)
    while True:
        sub = counts.iloc[:, kept]
        sub = sub.loc[sub.any(axis=1)]
        domain, coverage = domain_and_coverage(sub, settings)
        if settings["series"] == "anchored" or len(kept) == 2:
            if (coverage.iloc[[0, -1]].cells_in_domain <= 0).any():
                raise KarstDomainCoverageError(coverage)
            column = "majors" if settings["intermediate"] == "majors" else "cells_in_domain"
            for subindex in np.flatnonzero(coverage[column].to_numpy()[1:-1] == 0) + 1:
                i = kept[subindex]
                dropped.append({"sample": i + 1, "time": times[i], **coverage.iloc[subindex].to_dict()})
            excluded = {row["sample"] - 1 for row in dropped}
            kept = [i for i in kept if i not in excluded]
            break
        empty = coverage.index[coverage.cells_in_domain == 0]
        if empty.empty:
            break
        index = empty[0]
        pair = index // 2
        uncovered = pair + index % 2
        remove = (pair + 1 if uncovered == 0 else pair) if uncovered in (0, len(kept) - 1) else uncovered
        i = kept[remove]
        dropped.append({
            "sample": i + 1, "time": times[i], "timepoint": counts.columns[i], "cells": counts.iloc[:, i].sum(),
            "uncovered": counts.columns[kept[uncovered]],
        })
        kept.pop(remove)
    if not dropped:
        return data, domain, coverage, None
    dropped_frame = pd.DataFrame(dropped).sort_values("sample", ignore_index=True)
    sub = counts.iloc[:, kept]
    sub = sub.loc[sub.any(axis=1)]
    times = times[kept]
    reduced = KarstData(sub, sub.index.tolist(), times[-1], data.beta, data.metadata, times if len(kept) > 2 else None)
    domain, coverage = domain_and_coverage(sub, settings)
    return reduced, domain, coverage, dropped_frame


def _dispersion(fit):
    summary = fit.summary
    y = summary.loc[summary.major, summary.columns.str.startswith("count_")].to_numpy()[:, 1:-1]
    e = summary.loc[summary.major, summary.columns.str.startswith("expected_")].to_numpy()[:, 1:-1]
    e = e * (y.sum(axis=0) / e.sum(axis=0))
    return np.mean((y - e)**2 / np.maximum(e, 1e-8))


def _fit_candidates(data, settings, histories, control, domain, coverage):
    if len(histories) == 1:
        return _fit_single(data, dict(settings, history=histories[0]), control, domain, coverage)
    fits, records = [], []
    for history in histories:
        try:
            fit = _fit_single(data, dict(settings, history=history), control, domain, coverage)
        except KarstInnerFailure:
            fits.append(None)
            records.append(dict(history=history, objective=np.nan, posterior_objective=np.nan,
                                convergence=np.nan, inner_failure=True, posterior_restart=np.nan))
            continue
        fits.append(fit)
        records.append(dict(history=history, inner_failure=False, **{
            name: fit.diagnostics.get(name, np.nan)
            for name in ("objective", "posterior_objective", "convergence", "posterior_restart")
        }))
    profile = pd.DataFrame(records)
    if profile.inner_failure.all():
        raise KarstInnerFailure("No candidate history has a usable posterior mode.")
    usable = np.isfinite(profile.objective) & np.isfinite(profile.posterior_objective)
    if not usable.any():
        raise ValueError("No candidate history has a finite marginal likelihood.")
    accepted = usable & (profile.convergence == 0)
    if accepted.any() and get_option("choose_accepted"):
        usable = accepted
    chosen = profile.loc[usable, "posterior_objective"].idxmin()
    fit = fits[chosen]
    if not profile.inner_failure.any():
        profile = profile.drop(columns="inner_failure")
    if profile.posterior_restart.isna().all():
        profile = profile.drop(columns="posterior_restart")
    if settings.get("balance_level") == "free":
        profile["balance_level"] = [candidate.diagnostics.get("balance_level", np.nan) if candidate else np.nan for candidate in fits]
    fit.diagnostics["history_profile"] = profile
    return fit


def fit_sparse_joint(
    counts, dt=None, beta=None, birth_rate=1, birth_rate_log_sd=0.05, major_frequency=0.01,
    major_count=3, max_states=5000, min_cn=1, max_cn=None, shell_power=0.85, shell_sd=0.5,
    initial_frequency=(math.log(1e-5), 3), initial_state="balance", auto_level=1, history=None,
    balance_sd=1, balance_level="fixed", balance_level_sd=3, calibration=_DEFAULT_CALIBRATION,
    control=None, inner_control=None, times=None, series="anchored", intermediate="majors",
    intermediate_weight=None, *, labels=None,
):
    """Fit one fitness landscape to sparse karyotype counts from two or more time points.

    Supply a pandas DataFrame indexed by karyotype or a matrix with explicit
    labels. Each column is a time point. Use consecutive intervals in `dt` or
    absolute observation times in `times`. The fixed missegregation rate is
    `beta`. The default initial abundance prior uses mutation-selection
    balance. A list of positive values for `history` profiles that prior over
    several histories. The result includes relative fitness, fitted count
    expectations, diagnostics and the inputs needed to rebuild its posterior.
    """
    started = perf_counter()
    if isinstance(counts, KarstData):
        data = counts if beta is None else replace(counts, beta=validate_beta(beta))
    else:
        data = prepare_karst_data(counts, dt, 5e-5 if beta is None else beta, times=times, labels=labels)
    count_matrix = data.counts.to_numpy()
    n_time = count_matrix.shape[1]
    if (count_matrix.sum(axis=0) <= 0).any():
        raise ValueError("Every time point needs at least one cell.")
    cn = parse_karyotypes(data.counts.index)
    min_cn = integer(min_cn, "min_cn", 0)
    max_cn = cn.max() + 1 if max_cn is None else integer(max_cn, "max_cn", min_cn)
    if ((cn < min_cn) | (cn > max_cn)).any():
        raise ValueError("Observed karyotypes must lie within the copy-number bounds.")
    choices = {
        "initial_state": (initial_state, ("balance", "free", "auto")),
        "balance_level": (balance_level, ("fixed", "free")),
        "series": (series, ("anchored", "pairwise", "state_space")),
        "intermediate": (intermediate, ("majors", "all")),
    }
    for name, (value, allowed) in choices.items():
        if value not in allowed:
            raise ValueError(f"`{name}` must be one of {allowed}.")
    template = get_option("sparse_joint_template")
    two_time_template = n_time == 2 and template != "series"
    if initial_state == "auto" and not two_time_template:
        initial_state = "balance"
    if balance_level == "free" and not two_time_template and initial_state != "free":
        raise ValueError("`balance_level = 'free'` is implemented for two time points only.")
    initial_frequency = np.asarray(initial_frequency)
    if initial_frequency.shape != (2,):
        raise ValueError("`initial_frequency` must give a mean and standard deviation.")
    eps = np.finfo(float).eps
    settings = dict(
        birth_rate=scalar(birth_rate, "birth_rate", eps),
        birth_rate_log_sd=scalar(birth_rate_log_sd, "birth_rate_log_sd", eps),
        major_frequency=scalar(major_frequency, "major_frequency", 0, 1),
        major_count=scalar(major_count, "major_count", 2),
        max_states=integer(max_states, "max_states", 1), min_cn=min_cn, max_cn=int(max_cn),
        shell_power=scalar(shell_power, "shell_power", eps), shell_sd=scalar(shell_sd, "shell_sd", eps),
        initial_frequency=(scalar(initial_frequency[0], "initial_frequency mean"),
                           scalar(initial_frequency[1], "initial_frequency sd", eps)),
        initial_state=initial_state, balance_sd=scalar(balance_sd, "balance_sd", eps),
        balance_level=balance_level, balance_level_sd=scalar(balance_level_sd, "balance_level_sd", eps),
        calibration=_calibration(calibration), inner_control=dict(inner_control or {}),
        intermediate=intermediate, intermediate_weight=intermediate_weight,
        series=series if n_time > 2 else "anchored", template=template,
    )
    if intermediate_weight is not None:
        settings["intermediate_weight"] = scalar(intermediate_weight, "intermediate_weight", eps, 1)
    if series == "pairwise" and n_time > 2:
        settings["pair_history"] = get_option("pair_history") or "age"
    if series == "state_space" and n_time > 2:
        drift = get_option("state_space_drift")
        settings["drift"] = None if drift is None else scalar(drift, "drift", 0, lower_open=True)
    auto_level = scalar(auto_level, "auto_level")
    control = dict(control or {})
    original = data
    data, domain, coverage, dropped = _coverage_mask(data, settings)
    if initial_state == "auto":
        try:
            check = _fit_single(data, dict(settings, initial_state="balance", history=data.dt, balance_level="free"),
                                control, domain, coverage)
        except KarstInnerFailure:
            check = None
        diag = check.diagnostics if check else {}
        level = diag.get("balance_level", np.nan)
        settings["initial_state"] = "free" if level > auto_level else "balance"
    if history is None:
        divisions = np.array([30, 100, 300, 1000, 3000]) / birth_rate
        histories = np.r_[data.dt, divisions[divisions > data.dt]] if settings["initial_state"] == "balance" else [data.dt]
    else:
        histories = [scalar(h, "history", eps) for h in np.atleast_1d(history)]
        if not histories:
            raise ValueError("`history` must contain positive times.")
        if settings["initial_state"] == "free":
            histories = histories[:1]
    dispersion = np.nan
    if data.counts.shape[1] > 2 and series == "anchored" and intermediate == "majors" and intermediate_weight is None:
        try:
            initial_fit = _fit_single(data, dict(settings, history=min(histories), intermediate_weight=1),
                                      control, domain, coverage)
            dispersion = _dispersion(initial_fit)
        except KarstInnerFailure:
            pass
        settings["intermediate_weight"] = 1 / max(1, dispersion) if np.isfinite(dispersion) else 1
    fit = _fit_candidates(data, settings, histories, control, domain, coverage)
    if initial_state == "auto":
        fit.diagnostics["initial_state_check"] = pd.DataFrame([dict(
            level=level, level_sd=diag.get("balance_level_sd", np.nan), threshold=auto_level,
            chosen=settings["initial_state"], check_convergence=diag.get("convergence", np.nan),
        )])
    if np.isfinite(dispersion):
        fit.diagnostics["intermediate_dispersion"] = dispersion
    if dropped is not None:
        fit.input = original
        fit.diagnostics["dropped_samples"] = dropped
        known = set(fit.data.labels) | set(fit.summary.karyotype)
        missing = [label for label in original.labels if label not in known]
        fit.dropped.extend(missing)
        fit.diagnostics["dropped_cells"] += int(original.counts.loc[missing].to_numpy().sum())
        fit.diagnostics["dropped_karyotypes"] = len(fit.dropped)
        warnings.warn("Left out time points " + ", ".join(dropped.timepoint), KarstSamplesLeftOutWarning, stacklevel=2)
    fit.diagnostics["seconds"] = perf_counter() - started
    warn_if_unconverged(fit.diagnostics["optimisation"])
    return fit
