"""Simulate sparse karyotype counts on a lazy fitness landscape."""

import math

import numpy as np
import pandas as pd

from ._validate import integer, scalar
from .karyotypes import format_karyotypes, parse_karyotypes
from .rng import random_stream
from .simulate_l1_gp import GRFLandscape, predict_landscape_fitness


def _initial_karyotypes(founder, depth, lower, upper):
    nodes = {tuple(founder): None}
    for _ in range(depth):
        for state in tuple(nodes):
            for chromosome, copy in enumerate(state):
                for step in (-1, 1):
                    if lower <= copy + step <= upper:
                        nodes[state[:chromosome] + (copy + step,) + state[chromosome + 1:]] = None
    return np.asarray(list(nodes))


def _multinomial(stream, weights, size, overdispersion=0):
    # Match the repeated normalization in the R caller and in rmultinom.
    prob = weights / sum(map(float, weights))
    prob = prob / sum(map(float, prob))
    if overdispersion > 0:
        prob = np.array([stream.gamma(max(1e-10, overdispersion * p)) for p in prob])
        total = sum(map(float, prob))
        if not math.isfinite(total) or total <= 0:
            raise ValueError("`overdispersion` is too small: every Dirichlet proportion underflowed to zero.")
        prob /= total
    return stream.multinomial(size, (prob / sum(map(float, prob))).tolist())


def _draw_counts(stream, population, size, overdispersion=0):
    draw = _multinomial(stream, np.asarray(list(population.values()), dtype=float), size, overdispersion)
    return dict(zip(population, draw))


def _daughters(stream, karyotype, beta):
    copies = np.repeat(np.arange(len(karyotype)), karyotype)
    altered = int(stream.binomial(len(copies), beta))
    while altered == 0:
        altered = int(stream.binomial(len(copies), beta))
    picked = copies[stream.sample(len(copies), altered, False)]
    signs = 2 * np.array(stream.sample(2, altered, True)) - 1
    delta = np.zeros(len(karyotype), dtype=int)
    np.add.at(delta, picked, signs)
    return karyotype + delta, karyotype - delta


def simulate_sparse_counts(
    landscape, beta=0.00005, dt=1, n0=500, n1=500,
    detection_threshold=1, dropout_prob=0, seed=1, init_concentration=1.5,
    overdispersion=0, initial_population=None, initial_shell_depth=1,
    time_step=0.1, carrying_capacity=np.inf, passage_fraction=0.2, max_unique=50000,
):
    """Simulate counts at two time points and return the corresponding birth rates."""
    if not isinstance(landscape, GRFLandscape):
        raise ValueError("`landscape` must be a GRFLandscape object.")
    beta = scalar(beta, "beta", 0, 1)
    dt = scalar(dt, "dt", 0, lower_open=True)
    n0, n1 = integer(n0, "n0", 1), integer(n1, "n1", 1)
    detection_threshold = integer(detection_threshold, "detection_threshold", 0)
    dropout_prob = scalar(dropout_prob, "dropout_prob", 0, 1)
    init_concentration = scalar(init_concentration, "init_concentration", 0)
    overdispersion = scalar(overdispersion, "overdispersion", 0)
    if initial_population is None:
        initial_population = max(100, 5 * max(n0, n1))
    initial_population = integer(initial_population, "initial_population", 1)
    initial_shell_depth = integer(initial_shell_depth, "initial_shell_depth", 0)
    time_step = scalar(time_step, "time_step", 0, lower_open=True)
    if carrying_capacity != np.inf:
        carrying_capacity = scalar(carrying_capacity, "carrying_capacity", 1)
    passage_fraction = scalar(passage_fraction, "passage_fraction", 0, 1, lower_open=True, upper_open=True)
    max_unique = integer(max_unique, "max_unique", 1)
    stream = random_stream(seed)
    lower, upper = landscape.min_cn, landscape.max_cn
    init = _initial_karyotypes(landscape.founder, initial_shell_depth, lower, upper)
    weights = np.array([math.exp(-init_concentration * float(d)) for d in ((init - landscape.founder)**2).sum(axis=1)])
    draw = _multinomial(stream, weights, initial_population)
    population = {label: count for label, count in zip(format_karyotypes(init), draw) if count > 0}
    visited = set(population)
    y0 = _draw_counts(stream, population, n0, overdispersion)
    steps = math.ceil(dt / time_step)
    step = dt / steps
    for _ in range(steps):
        labels = list(population)
        karyotypes = parse_karyotypes(labels)
        fitness = predict_landscape_fitness(landscape, karyotypes)
        if np.any(~np.isfinite(fitness) | (fitness < 0)):
            raise ValueError("The division simulator requires finite non-negative birth rates.")
        next_counts = population.copy()
        for label, x, rate in zip(labels, karyotypes, fitness):
            count = population[label]
            faithful_prob = math.pow(1 - beta, int(x.sum()))
            net = 2 * faithful_prob - 1
            mean = math.expm1(float(rate) * net * step) / net if abs(net) > 1e-8 else rate * step
            divisions = min(count, int(stream.poisson(count * mean)))
            if divisions == 0:
                continue
            faithful = int(stream.binomial(divisions, faithful_prob))
            next_counts[label] += 2 * faithful - divisions
            for _ in range(divisions - faithful):
                for daughter in _daughters(stream, x, beta):
                    if np.all((lower <= daughter) & (daughter <= upper)):
                        child = format_karyotypes(daughter)[0]
                        next_counts[child] = next_counts.get(child, 0) + 1
                        visited.add(child)
        population = {label: next_counts[label] for label in sorted(next_counts) if next_counts[label] > 0}
        if not population:
            raise ValueError("All karyotypes went extinct during simulation.")
        if len(population) > max_unique:
            raise ValueError("Simulation exceeded `max_unique`.")
        if sum(population.values()) > carrying_capacity:
            retained = max(1, round(carrying_capacity * passage_fraction))
            population = {label: count for label, count in _draw_counts(stream, population, retained).items() if count > 0}
    y1 = _draw_counts(stream, population, n1, overdispersion)
    labels = list(dict.fromkeys([*y0, *y1]))
    counts = np.array([(y0.get(label, 0), y1.get(label, 0)) for label in labels])
    totals = counts.sum(axis=1)
    keep = (totals > 0) & (totals >= detection_threshold)
    if dropout_prob > 0:
        keep &= np.array([stream.uniform() >= dropout_prob for _ in keep])
    if not keep.any():
        keep[totals.argmax()] = True
    observed = np.asarray(labels)[keep].tolist()
    frame = pd.DataFrame(counts[keep], columns=["t0", "t1"], index=observed)
    p = totals[keep] / (n0 + n1)
    return {
        "counts": frame,
        "observed_labels": observed,
        "truth_observed": predict_landscape_fitness(landscape, observed),
        "landscape": landscape,
        "sparsity": {
            "observed_nodes": len(observed), "visited_nodes": len(visited),
            "live_nodes": len(population), "fraction_observed": len(observed) / len(visited),
            "effective_support_entropy": -sum(x * math.log(x) for x in p),
            "sampling_depth": n0 + n1, "detection_threshold": detection_threshold,
            "dropout_prob": dropout_prob,
        },
    }
