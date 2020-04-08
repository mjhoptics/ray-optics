#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2020 Michael J. Hayford
"""Various generators and utilities for producing 2d distributions

.. Created on Tue Mar 24 21:14:31 2020

.. codeauthor: Michael J. Hayford
"""

import math
import numpy as np


def grid_ray_generator(grid_rng):
    """Generator function to produce a 2d square regular grid.

    arguments:
        grid_rng: start, stop, num
        start: 2d numpy array of lower left grid coords
        stop: 2d numpy array of upper right grid coords
        num: the number of samples along each axis

    A sample input might be:
        grid_start = np.array([-1., -1.])
        grid_stop = np.array([1., 1.])
        grid_rng = grid_start, grid_stop, num_rays

    """
    start, stop, num = grid_rng
    sample_pt = np.array(start)
    step = np.array((stop - start)/(num - 1))
    for i in range(num):
        for j in range(num):
            yield np.array(sample_pt)
            sample_pt[1] += step[1]

        sample_pt[0] += step[0]
        sample_pt[1] = start[1]


def csd_grid_ray_generator(grid_rng):
    start = np.array(grid_rng[0])
    stop = grid_rng[1]
    num = grid_rng[2]
    step = np.array((stop - start)/(num - 1))
    for i in range(num):
        for j in range(num):
            xy = concentric_sample_disk(start, offset=False)
            yield xy
            start[1] += step[1]

        start[0] += step[0]
        start[1] = grid_rng[0][1]


def polar_grid_ray_generator(grid_rng):
    start = np.array(grid_rng[0])
    stop = grid_rng[1]
    num = grid_rng[2]
    step = np.array((stop - start)/(num - 1))
    for i in range(num):
        for j in range(num):
            yield np.array(start)
            start[1] += step[1]

        start[0] += step[0]
        start[1] = grid_rng[0][1]


# Using the above nested radical formula for g=phi_d
# or you could just hard-code it.
# phi(1) = 1.61803398874989484820458683436563
# phi(2) = 1.32471795724474602596090885447809
def phi(d):
    x = 2.0000
    for i in range(10):
        x = pow(1+x, 1/(d+1))
    return x


def R_2_quasi_random_generator(n):
    """A 2d sequence based on a R**2 quasi-random sequence

    See `The Unreasonable Effectiveness of Quasirandom Sequences 
    <http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/ >`
    """
    d = 2
    g = phi(d)
    alpha = np.zeros(d)
    for j in range(d):
        alpha[j] = pow(1/g, j+1) % 1

    # This number can be any real number.
    # Common default setting is typically seed=
    # But seed = 0.5 is generally better.
    seed = 0.5

    z = np.zeros((n, d))
    for i in range(n):
        z[i] = (seed + alpha*(i+1)) % 1
        yield z[i]


def concentric_sample_disk(u, offset=True):
    """Map a 2d unit square sample to the unit disk."""
    if offset:
        uOffset = 2*u - np.array([1, 1])
    else:
        uOffset = u

    if uOffset[0] == 0 and uOffset[1] == 0:
        return np.array([0, 0])

    if abs(uOffset[0]) > abs(uOffset[1]):
        r = uOffset[0]
        theta = np.pi/4 * (uOffset[1]/uOffset[0])
    else:
        r = uOffset[1]
        theta = np.pi/2 - np.pi/4 * (uOffset[0]/uOffset[1])

    return r*np.array([math.cos(theta), math.sin(theta)])


def create_generator(sampler, *sampler_args, mapper=None, **kwargs):
    def gen():
        for xy in sampler(*sampler_args):
            if mapper:
                yield mapper(xy, **kwargs)
            else:
                yield xy
    return gen()
