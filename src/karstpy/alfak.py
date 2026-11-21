"""ALFA-K landscape, posterior draw and cross-validation outputs."""

from dataclasses import dataclass
from pathlib import Path
import math
import pickle
import sys
import tempfile
import warnings

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.special import softmax

from ._problem import build_problem, design_matrix
from ._rds import RAttributed as _RAttributed
from ._validate import integer, scalar
from .data import karst_data_times
from .karyotypes import format_karyotypes, parse_karyotypes
from .rates import KarstRateFit
from .rng import normal_array, random_stream
from .sparse_joint import KarstFit, fit_sparse_joint, rebuild_posterior, self_probability
from .tmb import KarstInnerFailure


def alfak_neighbours(labels):
    """Return distinct viable one-step neighbours in ALFA-K order."""
    if len(labels) == 0:
        return []
    karyotypes = parse_karyotypes(labels)
    steps = np.repeat(np.eye(karyotypes.shape[1], dtype=int), 2, axis=0)
    steps[::2] *= -1
    moved = (karyotypes[:, None, :] + steps[None, :, :]).reshape(-1, karyotypes.shape[1])
    neighbors = format_karyotypes(moved[np.all(moved >= 1, axis=1)])
    original = set(labels)
    return [label for label in dict.fromkeys(neighbors) if label not in original]


@dataclass
class AlfakLandscape:
    fit: KarstFit
    level: float
    correct_efflux: bool

    def predict(self, karyotypes):
        """Return a one-column matrix of fitness on the ALFA-K scale."""
        cn = np.asarray(karyotypes)
        if cn.dtype.kind in "iuf":
            if np.any(np.abs(cn - np.rint(cn)) > 1e-8):
                warnings.warn("Non-integer copy numbers were rounded.", stacklevel=2)
            karyotypes = np.rint(cn)
        prediction = self.fit.predict(karyotypes)
        value = prediction.relative_fitness.to_numpy() + self.level
        if self.correct_efflux:
            value /= 2 * self_probability(parse_karyotypes(prediction.karyotype, unique=False), self.fit.data.beta) - 1
        return value[:, None]


@dataclass
class AlfakOutput:
    bootstrap_res: dict
    landscape: pd.DataFrame
    landscape_posterior_samples: np.ndarray
    xval: dict
    landscape_data: dict | None
    fit: KarstFit
    rebuild: dict
    bootstrap_retry_log: pd.DataFrame
    bootstrap_retry_summary: dict

    def __repr__(self):
        n, draws = self.landscape_posterior_samples.shape
        return f"AlfakOutput({n} karyotypes, {draws} posterior draws, cross-validation R2R={self.xval['R2R']:g})"


def _relative_draws(fit, posterior, labels, prediction, size, stream):
    o = posterior.objective
    latent = posterior.draw(size, stream)
    g_index, coef_index = o.latent_indices("g"), o.latent_indices("coef")
    inside = pd.Index(fit.summary.karyotype).get_indexer(labels)
    log_birth = np.zeros((len(labels), size))
    mask = inside >= 0
    log_birth[mask] = latent[g_index[inside[mask]]]
    trend = fit.trend
    observed = np.isin(labels, fit.dropped) & (trend["observed_level"] is not None)
    outside = ~mask & ~observed
    if outside.any():
        cn = parse_karyotypes(np.asarray(labels)[outside])
        design, _, _ = design_matrix(cn, trend["centre"], trend["major"], trend["shell_power"])
        log_birth[outside] = design @ latent[coef_index] + trend["scales"]["residual"] * normal_array(stream, (int(outside.sum()), size))
    if observed.any():
        log_birth[observed] = trend["observed_level"]["sd"] * normal_array(stream, (int(observed.sum()), size))
    fitness = prediction.fitness_mean.to_numpy()[:, None] * log_birth
    n_major = fit.diagnostics["n_major"]
    weights = fit.summary.loc[:n_major - 1, fit.summary.columns.str.startswith("count_")].sum(axis=1).to_numpy(dtype=float)
    weights /= weights.sum()
    reference = (weights * fit.summary.fitness_mean.to_numpy()[:n_major]) @ latent[g_index[:n_major]]
    extra = np.sqrt(np.maximum(0, prediction.relative_sd**2 - prediction.relative_sd_model**2)).to_numpy()
    return fitness - reference + extra[:, None] * normal_array(stream, (len(labels), size)), latent


def _mirrored(deviation, size):
    result = np.zeros((len(deviation), size))
    pairs = size // 2
    result[:, :2 * pairs:2] = deviation[:, :pairs]
    result[:, 1:2 * pairs:2] = -deviation[:, :pairs]
    return result


def _efflux_constant(relative, divisor, log_first, first_time, growth):
    def weights(k):
        return softmax(log_first - first_time * (relative - k) / divisor)

    x0 = weights(0)
    k = ((x0 * relative / divisor).sum() - growth) / (x0 / divisor).sum()
    if first_time == 0:
        return k

    def gap(offset):
        return weights(offset) @ ((relative - offset) / divisor) - growth

    width = 1.0
    for _ in range(50):
        a, b = k - width, k + width
        if np.sign(gap(a)) != np.sign(gap(b)):
            return brentq(gap, a, b, xtol=1e-12)
        width *= 2
    raise ValueError("Could not solve the efflux correction level.")


def _cross_validate(fit, problem, frequent, means, level, divisor):
    summary, data = fit.summary, problem.tmb_data
    label_index = pd.Index(summary.karyotype)
    n_major, n = problem.domain.n_major, len(summary)
    major_count = problem.domain.y.sum(axis=1)[:n_major].astype(float)
    depth = problem.domain.y.sum(axis=0)
    init_mean, init_sd = fit.settings["initial_frequency"]
    model_counts = problem.used.sum(axis=1)
    retained = set()
    results, tallies = [], dict(rebuilt=0, restarted=0, failed=0)
    for label in frequent:
        members = [label, *alfak_neighbours([label])]
        tested = [x for x in members if x not in retained and x in label_index]
        retained.update(members)
        if not tested:
            continue
        state = label_index.get_indexer(tested)
        major = np.setdiff1d(np.arange(n_major), state)
        if not major.size:
            continue
        weight = np.ones(n)
        weight[state] = 0
        prior, prior_sd = data["zeta_prior_mean"].copy(), data["zeta_prior_sd"].copy()
        prior[state] = np.log(data["y_start"].sum(axis=0)) + init_mean if prior.ndim == 2 else math.log(depth[0]) + init_mean
        prior_sd[state] = init_sd
        cn = problem.domain.karyotypes
        total = model_counts * weight
        centre = (cn * total[:, None]).sum(axis=0) / total.sum()
        design, _, _ = design_matrix(cn, centre, cn[major], fit.settings["shell_power"])
        changed = dict(data, weight=weight, zeta_prior_mean=prior, zeta_prior_sd=prior_sd, design=design)
        if "drift_weight" in data:
            changed["drift_weight"] = data["drift_weight"].copy()
            changed["drift_weight"][state] = depth[:-1] + depth[1:]
        objective = problem.objective(fit.settings, changed)
        stored = objective.pack(fit.posterior)
        fixed = stored[objective.fixed]
        try:
            try:
                found = objective.find_mode(fixed, stored[objective.random])
            except KarstInnerFailure:
                found = objective.find_mode(fixed, objective.initial[objective.random])
                tallies["restarted"] += 1
        except KarstInnerFailure:
            tallies["failed"] += 1
            continue
        tallies["rebuilt"] += 1
        held_out = (2 * np.exp(problem.kernel.log_self) - 1) * np.exp(objective.unpack(found.point)["g"])
        shift = major_count[major] @ (held_out[major] - summary.fitness_mean.to_numpy()[major]) / major_count[major].sum()
        relative = held_out[state] - shift - fit.trend["reference"]["fitness"]
        predicted = (relative + level) / divisor.loc[tested].to_numpy()
        results.extend(zip(tested, means.loc[tested].to_numpy(), predicted))
    frame = pd.DataFrame(results, columns=["k", "observation", "prediction"]).dropna()
    obs, pred = (frame[name].to_numpy(dtype=float) for name in ("observation", "prediction"))
    r2r = np.nan
    if len(obs) >= 2:
        observed = obs - obs.mean()
        predicted = pred - pred.mean()
        if np.any(observed):
            r2r = 1 - np.sum((predicted - observed)**2) / np.sum(observed**2)
    return dict(R2R=r2r, predictions=pred, observations=obs, xval_data=frame,
                tmp=np.column_stack([obs, pred])), tallies


def as_alfak(fit, minobs=20, nboot=45, n0=1e5, nb=1e7, correct_efflux=False,
             xval=True, landscape_data_output=False, first_time=0, *, seed=1):
    """Summarize the stored fit in ALFA-K's landscape and bootstrap formats."""
    if isinstance(fit, KarstRateFit):
        fit = fit.fits[fit.point]
    if not isinstance(fit, KarstFit):
        raise ValueError("`fit` must come from fit_sparse_joint or fit_sparse_joint_rates.")
    minobs = scalar(minobs, "minobs", 0)
    nboot = integer(nboot, "nboot", 1)
    first_time = scalar(first_time, "first_time")
    if (n0 is None) != (nb is None):
        raise ValueError("Give both `n0` and `nb`, or neither.")
    for name, flag in (("correct_efflux", correct_efflux), ("xval", xval), ("landscape_data_output", landscape_data_output)):
        if not isinstance(flag, (bool, np.bool_)):
            raise ValueError(f"`{name}` must be Boolean.")
    original = fit.input if fit.input is not None else fit.data
    counts = original.counts
    frequent = counts.index[counts.sum(axis=1) >= minobs].tolist()
    if not frequent:
        raise ValueError(f"No frequent karyotypes detected for minobs = {minobs:g}.")
    neighbours = alfak_neighbours(frequent)
    labels = list(dict.fromkeys([*frequent, *neighbours, *alfak_neighbours(neighbours)]))
    prediction = fit.predict(labels)
    divisor = 2 * self_probability(parse_karyotypes(labels), fit.data.beta) - 1 if correct_efflux else np.ones(len(labels))
    problem = build_problem(fit.data, fit.settings)
    posterior = rebuild_posterior(fit, problem)
    stream = random_stream(seed)
    deviation, draws = _relative_draws(fit, posterior, labels, prediction, max(nboot // 2, 1), stream)
    mirror = _mirrored(deviation, nboot)
    if nboot > 1:
        spread = np.sqrt((mirror**2).sum(axis=1) / (nboot - 1))
        ratio = np.divide(prediction.relative_sd.to_numpy(), spread, out=np.zeros_like(spread), where=spread > 0)
        mirror *= ratio[:, None]
    relative = prediction.relative_fitness.to_numpy()[:, None] + mirror
    n_frequent = len(frequent)
    zeta = posterior.mode["zeta"]
    zeta_index = posterior.objective.latent_indices("zeta")
    if problem.model == "sparse_joint_pairwise":
        member = problem.domain.member
        first_locations = np.full(member.shape, -1, dtype=int, order="F")
        first_locations.ravel(order="F")[member.ravel(order="F")] = zeta_index
        zeta_index = first_locations[:, 0]
        zeta = zeta[:, 0]
    inside = pd.Index(fit.summary.karyotype).get_indexer(frequent)
    rows = (inside >= 0) & (zeta_index[np.maximum(inside, 0)] >= 0)
    log_first = np.log(np.maximum(counts.loc[frequent].to_numpy(dtype=float)[:, 0], np.finfo(float).tiny))
    log_first[rows] = zeta[inside[rows]]
    log_count = np.repeat(log_first[:, None], nboot, axis=1)
    log_count[rows] += _mirrored(draws[zeta_index[inside[rows]]], nboot)
    if n0 is None:
        level = fit.trend["reference"]["fitness"]
    else:
        n0 = scalar(n0, "n0", 0, lower_open=True)
        nb = scalar(nb, "nb", 0, lower_open=True)
        growth = math.log(nb / n0) / np.diff(karst_data_times(original))[0]
        level = growth
        if correct_efflux:
            level = -_efflux_constant(prediction.relative_fitness.to_numpy()[:n_frequent],
                                      divisor[:n_frequent], log_first, first_time, growth)
    samples = (relative + level) / divisor[:, None]
    landscape = pd.DataFrame({
        "k": labels, "mean": samples.mean(axis=1), "median": np.median(samples, axis=1),
        "sd": samples.std(axis=1, ddof=1) if nboot > 1 else np.nan,
        "fq": np.isin(labels, frequent), "nn": np.isin(labels, neighbours),
    })
    final_fitness = samples[:n_frequent].T
    frequency = softmax(log_count - first_time * final_fitness.T, axis=0).T
    nn_fitness = samples[n_frequent:n_frequent + len(neighbours)].T
    bootstrap = dict(initial_fitness=final_fitness, final_fitness=final_fitness,
                     initial_frequencies=frequency, final_frequencies=frequency, nn_fitness=nn_fitness)
    retry_log = pd.DataFrame({
        "target_replicate": np.empty(0, dtype=np.int32), "attempt": np.empty(0, dtype=np.int32),
        "error_message": np.empty(0, dtype=object), "warning_message": np.empty(0, dtype=object),
    })
    retry_summary = dict(requested_nboot=nboot, successful_replicates=nboot, total_attempts=nboot,
                         failed_attempts=0, max_attempts=nboot * 20)
    cvs, tallies = _cross_validate(
        fit, problem, frequent if xval else [], landscape.set_index("k")["mean"], level, pd.Series(divisor, index=labels),
    )
    gradient = np.max(np.abs(posterior.objective.gradient(posterior.point)[posterior.objective.random]))
    model = AlfakLandscape(fit, level, correct_efflux)
    landscape_data = {"mean": model, "median": model} if landscape_data_output else None
    return AlfakOutput(bootstrap, landscape, samples, cvs, landscape_data, fit,
                       dict(gradient=gradient, xval_folds=tallies), retry_log, retry_summary)


def _r_constructor(value, converter):
    from rdata.conversion.to_r import build_r_object
    from rdata.parser import RObjectType

    attributes = value.attributes.copy()
    if isinstance(value.value, dict):
        items = value.value
        attributes["names"] = np.array(list(items), dtype=str)
        r_type = RObjectType.VEC
        contents = [converter.convert_to_r_object(x) for x in items.values()]
    else:
        array = np.asarray(value.value)
        if array.ndim > 1:
            attributes["dim"] = np.array(array.shape, dtype=np.int32)
        converted = converter.convert_to_r_object(array.ravel(order="F"))
        r_type, contents = converted.info.type, converted.value
    return build_r_object(r_type, value=contents, is_object="class" in attributes,
                          attributes=converter.convert_to_r_attributes(attributes))


def _r_named_array(values, names):
    return _RAttributed(values, {"names": np.asarray(names, dtype=str)})


def _r_matrix(values, rows=None, columns=None):
    if rows is None and columns is None:
        return np.asarray(values)
    attributes = {"dimnames": [None if x is None else np.asarray(x, dtype=str) for x in (rows, columns)]}
    return _RAttributed(values, attributes)


def _r_frame(frame):
    frame = frame.copy()
    frame.index = pd.RangeIndex(1, len(frame) + 1)
    return frame


def _r_data(data):
    values = dict(counts=_r_matrix(data.counts.to_numpy(dtype=np.int32), data.counts.index, data.counts.columns),
                  labels=np.asarray(data.labels), dt=float(data.dt), beta=float(data.beta), metadata=data.metadata)
    if data.times is not None:
        values["times"] = data.times
    return _RAttributed(values, {"class": "karst_data"})


def _r_fit(fit):
    trend = fit.trend.copy()
    for name in ("scales", "calibration", "observed_level"):
        value = trend[name]
        if value is not None:
            trend[name] = _r_named_array(np.asarray(list(value.values())), list(value))
    diagnostics = {
        name: (_r_frame(value) if isinstance(value, pd.DataFrame) else value)
        for name, value in fit.diagnostics.items() if name != "optimisation"
    }
    summary = _r_frame(fit.summary)
    for name in summary.columns[summary.columns.str.startswith("count_")]:
        summary[name] = summary[name].astype(float)
    settings = fit.settings.copy()
    settings["initial_frequency"] = _r_named_array(np.asarray(settings["initial_frequency"]), ["mean", "sd"])
    settings["calibration"] = trend["calibration"]
    values = dict(data=_r_data(fit.data), summary=summary, dropped=np.asarray(fit.dropped, dtype=str),
                  settings=settings, posterior=fit.posterior, trend=trend, diagnostics=diagnostics)
    if fit.input is not None:
        values["input"] = _r_data(fit.input)
    return _RAttributed(values, {"class": "karst_sparse_joint_fit"})


def _r_values(output):
    landscape = output.landscape
    frequent = landscape.k[landscape.fq].tolist()
    neighbours = landscape.k[landscape.nn].tolist()
    bootstrap = {
        name: _r_matrix(values, columns=neighbours if name == "nn_fitness" else frequent)
        for name, values in output.bootstrap_res.items()
    }
    xval = output.xval.copy()
    keys = xval["xval_data"].k
    for name in ("predictions", "observations"):
        if len(keys):
            xval[name] = _r_named_array(xval[name], keys)
    xval["xval_data"] = _r_frame(xval["xval_data"].astype({"k": str, "observation": float, "prediction": float}))
    xval["tmp"] = _r_matrix(xval["tmp"], keys if len(keys) else None, ["test_f", "est_f"])
    values = {
        "bootstrap_res.Rds": _RAttributed(bootstrap, dict(
            bootstrap_retry_log=_r_frame(output.bootstrap_retry_log),
            bootstrap_retry_summary=output.bootstrap_retry_summary,
        )),
        "landscape.Rds": _r_frame(landscape),
        "landscape_posterior_samples.Rds": output.landscape_posterior_samples,
        "xval.Rds": xval,
    }
    if output.landscape_data is not None:
        fit = _r_fit(output.fit)
        values["landscape_data.Rds"] = {
            name: _RAttributed(dict(fit=fit, level=float(model.level), correct_efflux=model.correct_efflux),
                              {"class": "karst_alfak_landscape"})
            for name, model in output.landscape_data.items()
        }
    return values


def write_alfak_outputs(x, outdir):
    """Write the ALFA-K RDS files and a Python pickle of the fitted model."""
    if not isinstance(x, AlfakOutput):
        raise ValueError("`x` must come from as_alfak.")
    if sys.version_info < (3, 11):
        from ._rds import write_rds

        options = {}
    else:
        try:
            from rdata import write_rds
            from rdata.conversion.to_r import DEFAULT_CLASS_MAP
        except ImportError as error:
            raise ImportError("Install karstpy[alfak] to write RDS files.") from error
        options = {"constructor_dict": {
            **DEFAULT_CLASS_MAP,
            _RAttributed: _r_constructor,
            np.generic: lambda scalar, converter: converter.convert_to_r_object(np.asarray(scalar)),
        }}
    dest = Path(outdir)
    dest.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=dest) as staging:
        staging = Path(staging)
        values = _r_values(x)
        for name, value in values.items():
            write_rds(staging / name, value, **options)
        with (staging / "karst_fit.pkl").open("wb") as sink:
            pickle.dump(x.fit, sink, protocol=pickle.HIGHEST_PROTOCOL)
        return [(staging / name).replace(dest / name) for name in [*values, "karst_fit.pkl"]]


def alfak_compat(yi, outdir, passage_times=None, minobs=20, nboot=45, n0=1e5, nb=1e7, pm=5e-5,
                 allow_noninteger_counts=False, correct_efflux=False, landscape_data_output=False, *, seed=1, **kwargs):
    """Fit ALFA-K inputs, write its output files and return cross-validation R2R."""
    if not isinstance(yi, dict) or "x" not in yi or not isinstance(yi["x"], pd.DataFrame):
        raise ValueError("`yi` needs a labelled count DataFrame under key `x`.")
    counts = yi["x"]
    if counts.shape[1] < 2:
        raise ValueError("`yi['x']` needs two or more time points.")
    values = counts.to_numpy(dtype=float)
    if np.any(abs(values - np.rint(values)) > np.finfo(float).eps**0.5):
        if not allow_noninteger_counts:
            raise ValueError("Use `allow_noninteger_counts=True` to round non-integer counts.")
        warnings.warn("Non-integer counts in `yi['x']` were rounded.", stacklevel=2)
        counts = counts.round()
    if passage_times is None:
        interval = scalar(yi.get("dt"), "yi['dt']")
        try:
            passage_times = counts.columns.to_numpy(dtype=float) * interval
        except (ValueError, TypeError) as error:
            raise ValueError("Give numeric column names or explicit `passage_times`.") from error
    fit = fit_sparse_joint(counts, times=passage_times, beta=pm, **kwargs)
    out = as_alfak(fit, minobs, nboot, n0, nb, correct_efflux,
                   landscape_data_output=landscape_data_output, first_time=passage_times[0], seed=seed)
    write_alfak_outputs(out, outdir)
    return out.xval["R2R"]
