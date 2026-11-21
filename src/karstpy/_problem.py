"""Build the domain, priors and inputs for the native likelihoods."""

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .data import karst_data_times
from .forward import DivisionKernel, sparse_division_kernel
from .karyotypes import format_karyotypes, parse_karyotypes
from .tmb import LaplaceObjective


@dataclass
class Domain:
    labels: list[str]
    karyotypes: np.ndarray
    major: np.ndarray
    y: np.ndarray
    dropped: list[str]
    dropped_cells: int
    member: np.ndarray | None = None
    pair_major: np.ndarray | None = None

    @property
    def n_major(self):
        return len(self.major)


class KarstDomainCoverageError(ValueError):
    """A required time point has no cell inside the model domain."""

    def __init__(self, coverage):
        self.coverage = coverage
        reasons = []
        for row in coverage.loc[coverage.cells_in_domain <= 0].itertuples():
            unit = f"Pair {row.pair}, time point" if hasattr(row, "pair") else "Time point"
            other = coverage[coverage.timepoint != row.timepoint]
            if hasattr(row, "pair"):
                other = other[other.pair == row.pair]
            source = ", ".join(other.timepoint[other.majors > 0])
            shared = "the time points share no karyotype" if row.shared == 0 else f"it shares {row.shared} karyotype(s)"
            reasons.append(
                f"{unit} {row.timepoint}: none of its {int(row.cells)} cells lies within one step of a major "
                f"karyotype. The majors come from {source}. {shared.capitalize()} and the median distance "
                f"to the nearest major is {row.median_distance_outside:g}."
            )
        super().__init__("A time point has no cell inside the model domain. " + " ".join(reasons))


def one_step_neighbours(karyotypes, lower, upper):
    karyotypes = np.asarray(karyotypes, dtype=np.int64)
    batches = []
    for chromosome in range(karyotypes.shape[1]):
        for direction in (-1, 1):
            moved = karyotypes.copy()
            moved[:, chromosome] += direction
            keep = (moved[:, chromosome] >= lower) & (moved[:, chromosome] <= upper)
            batches.append(moved[keep])
    return list(dict.fromkeys(format_karyotypes(np.concatenate(batches))))


def major_distance(karyotypes, majors):
    cn = np.asarray(karyotypes, dtype=np.int64)
    nearest = np.full(len(cn), np.inf)
    for major in majors:
        nearest = np.minimum(nearest, abs(cn - major).sum(axis=1))
    return nearest


def domain_from_labels(counts, labels, major_labels, **kwargs):
    cn = parse_karyotypes(labels)
    dropped = counts.index[~counts.index.isin(labels)].tolist()
    return Domain(
        labels, cn, cn[:len(major_labels)], counts.reindex(labels, fill_value=0).to_numpy(),
        dropped, int(counts.loc[dropped].to_numpy().sum()), **kwargs,
    )


def build_domain(counts, settings, major_columns=None):
    source = counts if major_columns is None else counts.iloc[:, major_columns]
    values = source.to_numpy()
    share = (values / values.sum(axis=0)).max(axis=1)
    totals = values.sum(axis=1)
    eligible = (share >= settings["major_frequency"]) & (totals >= settings["major_count"])
    if not eligible.any():
        eligible = totals == totals.max()
    for threshold in np.unique(share[eligible]):
        major = source.index[eligible & (share >= threshold)].tolist()
        neighbors = one_step_neighbours(parse_karyotypes(major), settings["min_cn"], settings["max_cn"])
        labels = list(dict.fromkeys([*major, *neighbors]))
        if len(labels) <= settings["max_states"]:
            break
    return domain_from_labels(counts, labels, major)


def _weighted_median(value, count):
    order = np.argsort(value, kind="stable")
    ranks = np.cumsum(count[order])
    indices = np.searchsorted(ranks, [(ranks[-1] - 1) // 2, ranks[-1] // 2], side="right")
    return value[order[indices]].mean()


def domain_coverage(counts, domain, settings):
    values = counts.to_numpy()
    inside = counts.index.isin(domain.labels)
    majors = counts.index.isin(domain.labels[:domain.n_major])
    distance = np.full(len(values), np.nan)
    if not inside.all():
        distance[~inside] = major_distance(parse_karyotypes(counts.index[~inside]), domain.major)
    totals = values.sum(axis=1)
    shared = (values > 0).sum(axis=1) > 1
    result = []
    for name, cells in zip(counts.columns, values.T):
        seen = cells > 0
        outside = ~inside & seen
        total = cells.sum()
        result.append({
            "timepoint": name, "cells": total, "cells_in_domain": cells[inside].sum(),
            "fraction": cells[inside].sum() / total, "karyotypes": seen.sum(),
            "major_candidates": ((cells / total >= settings["major_frequency"]) & (totals >= settings["major_count"])).sum(),
            "majors": (majors & seen).sum(), "shared": (shared & seen).sum(),
            "median_distance_outside": _weighted_median(distance[outside], cells[outside]) if outside.any() else np.nan,
        })
    return pd.DataFrame(result)


def pair_domains(counts, settings):
    domain_list, coverage_list = [], []
    for i in range(counts.shape[1] - 1):
        pair = counts.iloc[:, i:i + 2]
        pair = pair.loc[pair.any(axis=1)]
        domain = build_domain(pair, settings)
        coverage = domain_coverage(pair, domain, settings)
        coverage.insert(0, "pair", i + 1)
        domain_list.append(domain)
        coverage_list.append(coverage)
    majors = list(dict.fromkeys(label for domain in domain_list for label in domain.labels[:domain.n_major]))
    labels = list(dict.fromkeys([*majors, *(label for domain in domain_list for label in domain.labels)]))
    member = np.column_stack([np.isin(labels, domain.labels) for domain in domain_list])
    pair_major = np.column_stack([np.isin(labels, domain.labels[:domain.n_major]) for domain in domain_list])
    domain = domain_from_labels(counts, labels, majors, member=member, pair_major=pair_major)
    return domain, pd.concat(coverage_list, ignore_index=True)


def design_matrix(karyotypes, centre, majors, shell_power):
    difference = karyotypes - centre
    shell = major_distance(karyotypes, majors)
    values = np.column_stack([np.ones(len(karyotypes)), difference, difference**2, shell**shell_power])
    group = np.r_[0, np.repeat(1, difference.shape[1]), np.repeat(2, difference.shape[1]), 4]
    return values, group, shell


def initial_prior(counts, initial_frequency, depth=None, first_only=False):
    if depth is None:
        depth = counts.sum(axis=0)
    first = counts[:, 0]
    mean = np.full(len(first), np.log(depth[0]) + initial_frequency[0])
    if not first_only:
        for j in range(counts.shape[1] - 1, 0, -1):
            seen = counts[:, j] > 0
            mean[seen] = np.log(counts[seen, j] * depth[0] / depth[j])
    mean[first > 0] = np.log(first[first > 0])
    return mean, np.where(first > 0, 10, initial_frequency[1])


@dataclass
class Problem:
    model: str
    domain: Domain
    kernel: DivisionKernel
    centre: np.ndarray
    shell: np.ndarray
    used: np.ndarray
    tmb_data: dict
    parameters: dict
    random: dict
    free: dict
    scale_names: list[str]
    bounds: tuple[np.ndarray, np.ndarray]

    def objective(self, settings, data=None):
        return LaplaceObjective(
            self.model, self.tmb_data if data is None else data,
            self.parameters, self.random, self.free, settings["inner_control"],
        )


def model_name(counts, settings):
    if counts.shape[1] > 2:
        series = settings.get("series")
        return "sparse_joint_series" if series in ("anchored", None) else "sparse_joint_" + series
    if settings.get("template") == "series":
        return "sparse_joint_series"
    return "sparse_joint"


def domain_and_coverage(counts, settings):
    name = model_name(counts, settings)
    if name in ("sparse_joint_pairwise", "sparse_joint_state_space"):
        return pair_domains(counts, settings)
    columns = [0, -1] if counts.shape[1] > 2 and settings.get("intermediate") == "majors" else None
    domain = build_domain(counts, settings, columns)
    return domain, domain_coverage(counts, domain, settings)


def build_problem(data, settings, domain=None):
    model = model_name(data.counts, settings)
    if domain is None:
        domain, _ = domain_and_coverage(data.counts, settings)
    kernel = sparse_division_kernel(domain.karyotypes, data.beta)
    y = domain.y.astype(float)
    n, n_time = y.shape
    times = karst_data_times(data)
    interval, depth = np.diff(times), y.sum(axis=0)
    used = y.copy()
    weights = np.ones((n, n_time - 1))
    if model == "sparse_joint_series" and n_time > 2 and settings.get("intermediate") == "majors":
        used[domain.n_major:, 1:-1] = 0
        weights[domain.n_major:, :-1] = 0
        weights[:domain.n_major, :-1] = settings.get("intermediate_weight") or 1
    total = used.sum(axis=1)
    centre = (domain.karyotypes * total[:, None]).sum(axis=0) / total.sum()
    design, group, shell = design_matrix(domain.karyotypes, centre, domain.major, settings["shell_power"])
    intercept = np.log(settings["birth_rate"])
    names = ["dosage", "curvature", "interaction", "residual"]
    free_scale = [True, True, False, True]
    scale_prior = np.repeat(np.log(0.01), 4)
    scale_prior_sd = np.repeat(1.5, 4)
    lower = np.repeat(np.log(1e-4), 3)
    upper = np.repeat(np.log(0.5), 3)
    common = dict(
        edge_from=kernel.source, edge_to=kernel.target, edge_log_prob=kernel.log_prob,
        log_self=kernel.log_self, design=design, coef_group=group, reference=total.argmax(),
        intercept_mean=intercept, intercept_sd=settings["birth_rate_log_sd"], shell_sd=settings["shell_sd"],
        n_major=domain.n_major, balance=int(settings["initial_state"] == "balance"),
        history=settings["history"], balance_sd=settings["balance_sd"],
        balance_offset=np.zeros(n), weight=np.ones(n),
    )
    init_frequency = settings["initial_frequency"]
    mean, sd = initial_prior(used, init_frequency, depth, first_only=model == "sparse_joint_state_space")
    theta_name = "theta"
    matrices = {}
    if model == "sparse_joint":
        common.update(y0=y[:, 0], y1=y[:, 1], dt=times[-1])
        theta, theta_name = np.log(depth[1] / depth[0]), "theta1"
    elif model == "sparse_joint_series":
        common.update(y=y, tau=times[1:], sample_weight=weights)
        theta = np.log(depth[1:] / depth[0])
    elif model == "sparse_joint_pairwise":
        member = domain.member
        start, end = y[:, :-1] * member, y[:, 1:] * member
        mean, sd = np.zeros_like(start), np.ones_like(start)
        for k, active in enumerate(member.T):
            mean[active, k], sd[active, k] = initial_prior(np.column_stack([start[active, k], end[active, k]]), init_frequency)
        age = 0 if settings.get("pair_history") == "same" else times[:-1]
        common.update(
            y=y, y_start=start, y_end=end, interval=interval, member=member, pair_major=domain.pair_major,
            reference=(start + end).argmax(axis=0), balance_offset=np.zeros_like(start),
            pair_history=np.full(n_time - 1, settings["history"]) + age,
        )
        theta = np.log(end.sum(axis=0) / start.sum(axis=0))
    else:
        share = y / depth
        frequency = (share[:, :-1] + share[:, 1:]) / 2 + 1 / (depth[:-1] + depth[1:])
        common.update(y=y, interval=interval, drift_weight=1 / frequency, max_drift_var=4,
                      transition_offset=np.zeros_like(y[:, 1:]))
        matrices["zeta_later"] = np.where(y[:, 1:] > 0, np.log(np.maximum(y[:, 1:], 1)), np.log(depth[1:]) + init_frequency[0])
        theta = np.log(depth[1:] / depth[:-1])
        names.append("drift")
        drift = settings.get("drift")
        scale_prior, scale_prior_sd = np.r_[scale_prior, np.log(1e-4)], np.r_[scale_prior_sd, 2.5]
        free_scale.append(drift is None)
        if drift is None:
            lower, upper = np.r_[lower, np.log(1e-9)], np.r_[upper, 0]
    common.update(zeta_prior_mean=mean, zeta_prior_sd=sd, log_scale_prior_mean=scale_prior, log_scale_prior_sd=scale_prior_sd)
    log_scale = scale_prior.copy()
    if model == "sparse_joint_state_space" and not free_scale[-1]:
        log_scale[-1] = np.log(drift)
    parameters = dict(zeta=mean, **matrices, g=np.full(n, intercept), coef=np.r_[intercept, np.zeros(design.shape[1] - 1)])
    parameters[theta_name] = theta
    parameters["log_scale"] = log_scale
    random = {name: np.ones_like(value, dtype=bool) for name, value in parameters.items() if name != "log_scale"}
    if model == "sparse_joint_pairwise":
        random["zeta"] = domain.member
    if model == "sparse_joint":
        level_free = settings.get("balance_level") == "free" and settings["initial_state"] == "balance"
        common.update(level_free=int(level_free), level_sd=settings.get("balance_level_sd", 1))
        parameters["balance_level"] = 0.0
        if level_free:
            random["balance_level"] = np.array(True)
    return Problem(model, domain, kernel, centre, shell, used, common, parameters, random,
                   {"log_scale": free_scale}, names, (lower, upper))
