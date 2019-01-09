#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
"""
.. Created on Wed Apr 18 11:04:53 2018

.. codeauthor: Michael J. Hayford
"""

import numpy as np


def line(p1, p2):
    A = (p1[1] - p2[1])
    B = (p2[0] - p1[0])
    C = (p1[0]*p2[1] - p2[0]*p1[1])
    return A, B, -C


def intersection(L1, L2):
    D = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    Dy = L1[0] * L2[2] - L1[2] * L2[0]
    if D != 0:
        x = Dx / D
        y = Dy / D
        return x, y
    else:
        return False


def get_intersect(a1, a2, b1, b2):
    """
    Returns the point of intersection of the lines passing through a2,a1 and b2,b1.
    a1: [x, y] a point on the first line
    a2: [x, y] another point on the first line
    b1: [x, y] a point on the second line
    b2: [x, y] another point on the second line
    """
    s = np.vstack([a1, a2, b1, b2])      # s for stacked
    h = np.hstack((s, np.ones((4, 1))))  # h for homogeneous
    l1 = np.cross(h[0], h[1])            # get first line
    l2 = np.cross(h[2], h[3])            # get second line
    x, y, z = np.cross(l1, l2)           # point of intersection
    if z == 0:                           # lines are parallel
        return (float('inf'), float('inf'))
    return (x/z, y/z)


def do_intersect(a1, a2, b1, b2, soln, delta):
    xy = get_intersect(a1, a2, b1, b2)
    d = xy - a1
    if abs(d[0]/delta[0] - d[1]/delta[1]) > 1e-14:
        print("t", d[0]/delta[0], d[1]/delta[1])
    if delta[0] != 0.0:
        t = d[0]/delta[0]
    else:
        t = d[1]/delta[1]
    if soln[0] < t and t < 1.0:
        soln = t, xy
    return soln


def intersect_with_3lines(pt, wht, bg, gr, rb):
#    print("wht_dist", wht_dist)
    soln = 0.0, pt
    delta = wht - pt
#    print("bg")
    soln = do_intersect(pt, wht, *bg, soln, delta)
#    print("gr")
    soln = do_intersect(pt, wht, *gr, soln, delta)
#    print("rb")
    soln = do_intersect(pt, wht, *rb, soln, delta)
    return soln[1]
