# coding=utf-8
"""
.. moduleauthor:: Shuhua Gao

This :mod:`copt` module provides advanced routines for numerical constant optimization in gene expression programming.
It is known that both genetic programming (GP) and GEP can usually find a near-optimal solution quite quickly. However,
they both have difficulty in optimizing the numerical constants integrated in the model, especially when the constants
are continuous real numbers. To optimize the constants efficiently, a local optimizer is often incorporated, i.e., GEP
mainly aims to formulate the *true* structure of the model, while the embedded local optimizer is dedicated to improving
the constants existent in the individuals.

To reduce computational cost, such local optimization is applied every *k* generations and often only applied to
selected individuals. Accordingly, two arguments *opt_period* and *opt_selector* are required by the algorithms.

Such kind of algorithms combining evolutionary algorithms and local improvement procedures are often referred to as
`memetic algorithm <https://en.wikipedia.org/wiki/Memetic_algorithm>`_ in literature.

A common interface for a local optimizer is ``optimize(population, selector, generation, **kwargs) -> population``.
"""


def gep_opt(population, toolbox, n_generations=100, n_elites=1,
            stats=None, hall_of_fame=None, verbose=__debug__,
            optimizer=None, opt_period=None, opt_selector=None, **kwargs):
    pass


def _hill_climbing(population, selector, generation, **kwargs):
    pass


def gep_hill_climbing(population, toolbox, n_generations=100, n_elites=1,
            stats=None, hall_of_fame=None, verbose=__debug__,
            opt_period=30, opt_selector='best 10', max_step=1, rnc_gen=None):
    pass