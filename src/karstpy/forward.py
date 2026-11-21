"""Sparse first- and second-order division dynamics."""

from dataclasses import dataclass
from itertools import product
import math

import numpy as np
from scipy.sparse import coo_matrix
from scipy.stats import poisson

from . import _core


@dataclass
class DivisionKernel:
    """Directed one-step edges and single-daughter retention probabilities."""

    n: int
    source: np.ndarray
    target: np.ndarray
    log_prob: np.ndarray
    log_self: np.ndarray


def sparse_division_kernel(karyotypes, beta):
    states = [tuple(map(int, state)) for state in np.asarray(karyotypes)]
    index = {state: i for i, state in enumerate(states)}
    by_copy = {cn for state in states for cn in state}
    stay = {cn: math.log(_core.pij(cn, cn, beta)) for cn in by_copy}
    log_self = np.array([sum(stay[cn] for cn in state) for state in states])
    source, target, log_prob = [], [], []
    for chromosome in range(len(states[0])):
        for step in (-1, 1):
            move = {cn: _core.pij(cn, cn + step, beta) for cn in by_copy}
            for i, state in enumerate(states):
                cn = state[chromosome]
                shifted = state[:chromosome] + (cn + step,) + state[chromosome + 1:]
                j = index.get(shifted)
                if j is not None and move[cn] > 0:
                    source.append(i)
                    target.append(j)
                    log_prob.append(math.log(move[cn]) + log_self[i] - stay[cn])
    return DivisionKernel(
        len(states), np.array(source, dtype=np.int32), np.array(target, dtype=np.int32),
        np.array(log_prob), log_self,
    )


def log_phi1(a, b, dt):
    """Log integral of exp(a * (dt - s) + b * s) over the interval."""
    gap = abs(np.asarray(a) - np.asarray(b)) * dt
    ratio = np.where(gap < 1e-6, 1 - gap / 2, -np.expm1(-gap) / np.maximum(gap, np.finfo(float).tiny))
    return np.maximum(a, b) * dt + math.log(dt) + np.log(ratio)


def log_phi2(a, b, c, dt):
    """Log integral for a two-edge path through three growth rates."""
    top = np.maximum(np.maximum(a, b), c)
    low = np.minimum(np.minimum(a, b), c)
    middle = a + b + c - top - low
    u, v = (top - middle) * dt, (top - low) * dt

    def h(z):
        return np.where(z < 1e-6, 1 - z / 2 + z**2 / 6, -np.expm1(-z) / np.maximum(z, np.finfo(float).tiny))

    z = (u + v) / 2
    minus_dh = np.where(
        z < 1e-4, 0.5 - z / 3 + z**2 / 8,
        (1 - np.exp(-z) * (1 + z)) / np.maximum(z, 1e-4)**2,
    )
    ratio = np.where(v - u < 1e-6 * np.maximum(1, v), minus_dh, (h(u) - h(v)) / np.maximum(v - u, np.finfo(float).tiny))
    return top * dt + 2 * math.log(dt) + np.log(ratio)


def row_log_sum_exp(value, group, size):
    peaks = np.full(size, -np.inf)
    np.maximum.at(peaks, group, value)
    totals = np.zeros(size)
    finite = np.isfinite(peaks[group])
    np.add.at(totals, group[finite], np.exp(value[finite] - peaks[group[finite]]))
    valid = totals > 0
    peaks[valid] += np.log(totals[valid])
    return peaks


def closed_form_forward(log_n0, log_birth, kernel, dt, order=1):
    """Predict log abundance using one- or two-edge paths on the domain."""
    source, target = kernel.source, kernel.target
    log_n0, log_birth = np.asarray(log_n0), np.asarray(log_birth)
    net = (2 * np.exp(kernel.log_self) - 1) * np.exp(log_birth)
    log_flux = math.log(2) + kernel.log_prob + log_birth[source]
    terms = [log_n0 + net * dt, log_flux + log_n0[source] + log_phi1(net[target], net[source], dt)]
    groups = [np.arange(kernel.n), target]
    if order >= 2 and source.size:
        incoming = [[] for _ in range(kernel.n)]
        outgoing = [[] for _ in range(kernel.n)]
        for edge, (i, j) in enumerate(zip(source, target)):
            outgoing[i].append(edge)
            incoming[j].append(edge)
        paths = [pair for enters, leaves in zip(incoming, outgoing) for pair in product(leaves, enters)]
        if paths:
            second, first = np.asarray(paths).T
            i, k, j = source[first], target[first], target[second]
            terms.append(log_flux[first] + log_flux[second] + log_n0[i] + log_phi2(net[j], net[k], net[i], dt))
            groups.append(j)
    return row_log_sum_exp(np.concatenate(terms), np.concatenate(groups), kernel.n)


def exact_forward_reference(log_n0, log_birth, kernel, dt, tolerance=1e-13):
    """Reference matrix exponential by uniformisation on the same domain."""
    source, target = kernel.source, kernel.target
    birth = np.exp(log_birth)
    net = (2 * np.exp(kernel.log_self) - 1) * birth
    shift = net.max()
    rate = (shift - net).max() + np.finfo(float).eps
    rows = np.arange(kernel.n)
    flux = 2 * np.exp(kernel.log_prob) * birth[source]
    operator = coo_matrix(
        (np.r_[1 + (net - shift) / rate, flux / rate], (np.r_[rows, target], np.r_[rows, source])),
        shape=(kernel.n, kernel.n),
    ).tocsr()
    offset = np.max(log_n0)
    term = np.exp(log_n0 - offset)
    mean = rate * dt
    steps = math.ceil(mean + 12 * math.sqrt(mean) + 30)
    weight = poisson.pmf(np.arange(steps + 1), mean)
    total = weight[0] * term
    for jump in range(1, steps + 1):
        term = operator @ term
        add = weight[jump] * term
        total += add
        if jump > mean and add.sum() <= tolerance * total.sum():
            break
    return np.log(total) + offset + shift * dt
