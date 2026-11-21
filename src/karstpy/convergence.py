"""Scale optimisation and checks of the resulting inner modes."""

import warnings

import numpy as np
from scipy.optimize import minimize
from scipy.sparse.linalg import eigsh, ArpackNoConvergence

from .options import get_option
from .tmb import KarstInnerFailure


class KarstConvergenceWarning(UserWarning):
    """The outer optimum did not pass the numerical convergence checks."""


def _smallest_eigenpair(matrix):
    if matrix.shape[0] <= 250:
        values, vectors = np.linalg.eigh(matrix.toarray())
        return values[0], vectors[:, 0]
    try:
        values, vectors = eigsh(matrix, k=1, sigma=0, which="LM", v0=np.ones(matrix.shape[0]), tol=1e-7)
        return values[0], vectors[:, 0]
    except (ArpackNoConvergence, RuntimeError):
        return np.nan, None


def _flat_direction(vector, labels):
    if vector is None or labels is None:
        return None, None
    energy = vector**2
    order = np.argsort(-energy, kind="stable")
    result = labels.iloc[order].copy()
    result["share"] = energy[order] / energy.sum()
    composition = result.groupby("variable", sort=False, dropna=False)["share"].sum().sort_values(ascending=False)
    return result.iloc[:5].reset_index(drop=True), composition


def _check(objective, point, value, lower, upper, labels):
    saved = objective.state()
    try:
        refit = objective.fn(point)
        mode = objective.last_par
        eigenvalue, direction = _smallest_eigenpair(objective.hessian(mode)) if np.isfinite(refit) else (np.nan, None)
        joint_fit = objective.value(mode)
        objective.random_start = objective.initial[objective.random]
        fixed = objective.fn(point)
        joint_fixed = objective.value(objective.last_par) if np.isfinite(fixed) else np.nan
        objective.random_start = mode[objective.random]
        derivative = np.full(point.size, np.nan)
        for i in range(point.size):
            plus, minus = point.copy(), point.copy()
            plus[i], minus[i] = min(upper[i], point[i] + 1e-3), max(lower[i], point[i] - 1e-3)
            derivative[i] = (objective.fn(plus) - objective.fn(minus)) / (plus[i] - minus[i])
        blocked = ((point <= lower + 1e-8) & (derivative > 0)) | ((point >= upper - 1e-8) & (derivative < 0))
        gradient = np.max(abs(np.where(blocked, 0, derivative)))
        reasons = []
        if not np.isfinite(value) or not np.isfinite(refit):
            reasons.append("objective not finite")
        elif abs(refit - value) > 1e-3:
            reasons.append("objective not reproduced at the fit's mode")
        if not eigenvalue >= 2e-3:
            reasons.append("near-singular inner Hessian")
        if joint_fixed < joint_fit - 1e-3:
            reasons.append("a better inner mode exists")
        if not np.isfinite(gradient):
            reasons.append("objective not finite next to the point")
        elif gradient > 0.05:
            reasons.append("gradient not small")
        flat, composition = _flat_direction(direction, labels) if reasons else (None, None)
        return dict(reasons=reasons, objective_difference=fixed - value, joint_difference=joint_fixed - joint_fit,
                    min_eigen=eigenvalue, projected_gradient=gradient, flat_direction=flat, flat_composition=composition)
    finally:
        objective.restore(saved)


def _compass(objective, start, lower, upper):
    point = np.clip(start, lower, upper)
    objective.random_start = objective.initial[objective.random]
    best = objective.fn(point)
    best_mode = objective.last_par.copy()
    evaluations = 1
    steps = np.eye(point.size)
    step = 0.5
    while step >= 0.005 and evaluations < 200:
        trial_improved = False
        for direction in steps:
            for sign in (-1, 1):
                candidate = np.clip(point + sign * step * direction, lower, upper)
                if np.array_equal(candidate, point):
                    continue
                value = objective.fn(candidate)
                evaluations += 1
                if np.isfinite(value) and (not np.isfinite(best) or value < best - 1e-6):
                    spectrum, _ = _smallest_eigenpair(objective.hessian(objective.last_par))
                    if spectrum >= 2e-3:
                        point, best, best_mode = candidate, value, objective.last_par.copy()
                        objective.random_start = best_mode[objective.random]
                        trial_improved = True
        if not trial_improved:
            step /= 2
    return point, best, best_mode, evaluations


def optimise_scales(objective, bounds, control, labels=None):
    """Optimize the Laplace objective with the approved bounded fallback solver."""
    lower, upper = bounds

    def evaluate(point):
        value = objective.fn(point)
        if not np.isfinite(value):
            return np.inf, np.zeros_like(point)
        try:
            gradient = objective.gr(point)
        except KarstInnerFailure:
            return np.inf, np.zeros_like(point)
        return value, gradient

    first = minimize(
        evaluate, objective.par, jac=True, method="L-BFGS-B", bounds=list(zip(lower, upper)),
        options={"maxiter": int(control.get("iter.max", 100)), "maxfun": int(control.get("eval.max", 200)),
                 "ftol": float(control.get("rel.tol", 1e-10)), "gtol": 1e-7, "maxls": 30},
    )
    point, value, stage = first.x, first.fun, "first"
    if not np.isfinite(value):
        point, stage = objective.par, "none"
        value = objective.fn(point)
    check = _check(objective, point, value, lower, upper, labels)
    if stage == "none":
        check["reasons"].insert(0, "scale optimisation failed")
    restarts, evaluations = 0, first.nfev
    if check["reasons"] and get_option("convergence_restarts"):
        saved = objective.state()
        candidates = [(point, value, objective.last_par_best, stage, check)]
        for name, start in (("end_of_first_run", point), ("prior_centred", objective.par)):
            new_point, new_value, new_mode, calls = _compass(objective, start, lower, upper)
            restarts += 1
            evaluations += calls
            if np.isfinite(new_value):
                objective.random_start = new_mode[objective.random]
                new_check = _check(objective, new_point, new_value, lower, upper, labels)
                candidates.append((new_point, new_value, new_mode, name, new_check))
                if not new_check["reasons"]:
                    break
        objective.restore(saved)
        point, value, mode, stage, check = min(
            candidates, key=lambda entry: (bool(entry[4]["reasons"]), entry[1] if np.isfinite(entry[1]) else np.inf)
        )
        objective.last_par_best, objective.value_best = mode, value
    passed = not check["reasons"]
    check.update(status="converged" if passed else "not converged", stage=stage, restarts=restarts,
                 evaluations=evaluations, optimizer="L-BFGS-B", optimizer_status=int(first.status),
                 optimizer_message=str(first.message), nlminb_convergence=None, nlminb_message=None)
    message = "converged" if passed else "not converged: " + "; ".join(check["reasons"])
    return point, value, int(not passed), message, check


def warn_if_unconverged(check):
    if check["reasons"]:
        message = "Scale optimisation did not converge: " + "; ".join(check["reasons"])
        flat = check["flat_direction"]
        if flat is not None and len(flat):
            leading = flat.iloc[0]
            message += f". The flattest direction is led by {leading.variable}"
            if isinstance(leading.karyotype, str):
                message += f" of karyotype {leading.karyotype}"
        warnings.warn(message + ".", KarstConvergenceWarning, stacklevel=3)
