"""Sparse Laplace approximation for the four native count likelihoods."""

from dataclasses import dataclass
import math

import numpy as np
from scipy.sparse import coo_matrix

from . import _core


class KarstInnerFailure(RuntimeError):
    """The latent posterior mode could not be found."""


@dataclass
class InnerMode:
    point: np.ndarray
    joint: float
    hessian: object
    hessian_values: np.ndarray
    laplace: float


class LaplaceObjective:
    """Parameter packing, sparse derivatives, and inner-solver state."""

    def __init__(self, model, data, parameters, random, fixed, inner_control=None):
        packed = {name: np.atleast_1d(value) for name, value in parameters.items()}
        constants = {name: np.atleast_1d(value) for name, value in data.items()}
        self.tape = _core.ModelTape(model, constants, packed)
        self.initial = np.asarray(self.tape.start)
        self.slices = {}
        self.shapes = {}
        offset = 0
        random_index, fixed_index = [], []
        for name, value in parameters.items():
            array = np.asarray(value)
            indices = np.arange(offset, offset + array.size)
            self.slices[name] = indices
            self.shapes[name] = array.shape
            if name in random:
                random_index.extend(indices[np.ravel(random[name], order="F")])
            if name in fixed:
                fixed_index.extend(indices[np.asarray(fixed[name], dtype=bool)])
            offset += array.size
        self.random = np.asarray(random_index, dtype=int)
        self.fixed = np.asarray(fixed_index, dtype=int)
        self.par = self.initial[self.fixed].copy()
        self.row, self.column = (np.asarray(index, dtype=np.int32) for index in self.tape.pattern)
        reverse_index = np.full(offset, -1, dtype=np.int32)
        reverse_index[self.random] = np.arange(self.random.size)
        self.random_entries = (reverse_index[self.row] >= 0) & (reverse_index[self.column] >= 0)
        self.hrows = reverse_index[self.row[self.random_entries]]
        self.hcols = reverse_index[self.column[self.random_entries]]
        self.factor = _core.Cholesky(self._matrix(np.ones(self.row.size)))
        self.last_par = self.last_par1 = self.last_par2 = self.last_par_ok = self.initial.copy()
        self.last_par_best = self.initial.copy()
        self.value_best = np.inf
        self.random_start = None
        self.inner_control = dict(inner_control or {})

    def pack(self, values):
        point = self.initial.copy()
        for name, value in values.items():
            if name in self.slices:
                point[self.slices[name]] = np.ravel(value, order="F")
        return point

    def unpack(self, point):
        return {name: point[index].reshape(self.shapes[name], order="F").copy()
                for name, index in self.slices.items()}

    def latent_indices(self, name):
        return np.flatnonzero(np.isin(self.random, self.slices[name]))

    def value(self, point):
        return self.tape.value(point)

    def gradient(self, point):
        return self.tape.gradient(point)

    def report(self, point=None):
        return self.tape.report(self.last_par if point is None else point)

    def _matrix(self, values):
        n = self.random.size
        return coo_matrix((values[self.random_entries], (self.hrows, self.hcols)), shape=(n, n)).tocsc()

    def hessian(self, point, random=True):
        values = self.tape.hessian_values(point)
        if random:
            return self._matrix(values)
        n = self.initial.size
        return coo_matrix((values, (self.row, self.column)), shape=(n, n)).tocsc()

    def state(self):
        return (self.last_par, self.last_par1, self.last_par2, self.last_par_ok,
                self.last_par_best, self.value_best, self.random_start)

    def restore(self, state):
        (self.last_par, self.last_par1, self.last_par2, self.last_par_ok,
         self.last_par_best, self.value_best, self.random_start) = state

    def find_mode(self, fixed, start=None):
        if start is None:
            start = self.random_start
        if start is None:
            start = self.last_par_best[self.random]
        x = self.initial.copy()
        x[self.fixed], x[self.random] = fixed, start
        config = self.inner_control
        tol = config.get("tol", 1e-10)
        grad_tol, step_tol = config.get("grad.tol", tol), config.get("step.tol", tol)
        u, u0, power = config.get("ustep", 1.0), config.get("u0", 1e-4), config.get("power", 0.5)
        joint = self.value(x)
        history, failed_checks = [], 0
        for iteration in range(int(config.get("maxit", 200))):
            gradient = self.gradient(x)[self.random]
            magnitude = np.max(np.abs(gradient))
            if not np.isfinite(magnitude) or magnitude > config.get("mgcmax", 1e60):
                raise KarstInnerFailure("Inner gradient is non-finite or too steep.")
            if magnitude < grad_tol:
                break
            hessian = self.hessian(x)
            old = x[self.random].copy()
            if config.get("smartsearch", True):
                min_diag = hessian.diagonal().min()
                u = u0 + (1 - u0) * u**power
                trial = x.copy()
                for _ in range(2000):
                    shift = 1 / u - 1
                    if min_diag < 0 and shift <= -min_diag:
                        u = min(u, 1 / (1 - min_diag))
                    elif self.factor.factorize(hessian, shift):
                        trial[self.random] = old - self.factor.solve(gradient).reshape(-1)
                        value = self.value(trial)
                        if np.isfinite(value) and value <= joint + math.sqrt(np.finfo(float).eps):
                            x, joint = trial, value
                            break
                    u = 1 - (u0 + (1 - u0) * (1 - u)**power) if u > 1e-10 else (1 - u0) * power * u
                    if u <= np.finfo(float).tiny:
                        raise KarstInnerFailure("No admissible inner Newton step.")
                else:
                    raise KarstInnerFailure("No admissible inner Newton step.")
            else:
                if not self.factor.factorize(hessian):
                    raise KarstInnerFailure("Inner Hessian is not positive definite.")
                x[self.random] = old - config.get("alpha", 1) * self.factor.solve(gradient).reshape(-1)
                joint = self.value(x)
            history.append(joint)
            if iteration >= 10 and history[-10] - joint < config.get("tol10", 1e-3):
                if self.factor.factorize(self.hessian(x)):
                    break
                failed_checks += 1
            if np.linalg.norm(x[self.random] - old) < step_tol:
                break
            if failed_checks > 5:
                raise KarstInnerFailure("Too many inner Newton attempts.")
        values = self.tape.hessian_values(x)
        hessian = self._matrix(values)
        if not np.isfinite(joint) or not np.isfinite(values[self.random_entries]).all() or not self.factor.factorize(hessian):
            raise KarstInnerFailure("Inner Newton search did not find a finite positive-definite mode.")
        laplace = joint + 0.5 * (self.factor.log_det() - self.random.size * math.log(2 * math.pi))
        return InnerMode(x, joint, hessian, values, laplace)

    def fn(self, fixed):
        try:
            mode = self.find_mode(fixed)
        except KarstInnerFailure:
            return np.nan
        self.last_par = self.last_par1 = self.last_par_ok = mode.point
        if mode.laplace < self.value_best:
            self.last_par_best, self.value_best = mode.point, mode.laplace
        return mode.laplace

    def gr(self, fixed):
        mode = self.find_mode(fixed)
        self.last_par = self.last_par2 = mode.point
        weights = np.zeros(self.row.size)
        weights[self.random_entries] = 0.5 * self.factor.inverse_weights(self.hrows, self.hcols)
        gradient = self.gradient(mode.point) + self.tape.hessian_reverse(mode.point, weights)
        adjustment = self.factor.solve(gradient[self.random]).reshape(-1)
        n = self.initial.size
        joint_hessian = coo_matrix((mode.hessian_values, (self.row, self.column)), shape=(n, n)).tocsr()
        return gradient[self.fixed] - joint_hessian[self.fixed][:, self.random] @ adjustment
