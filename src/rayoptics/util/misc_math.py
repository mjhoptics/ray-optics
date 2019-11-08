#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" miscellaneous functions for working with numpy vectors and floats

.. Created on Wed May 23 15:27:06 2018

.. codeauthor: Michael J. Hayford
"""
import numpy as np
from numpy.linalg import norm
from math import sqrt


def normalize(v):
    """ return normalized version of input vector v """
    length = norm(v)
    if length == 0.0:
        return v
    else:
        return v/length


def distance_sqr_2d(pt0, pt1):
    """ return distance squared between 2d points pt0 and pt1 """
    return (pt0[0] - pt1[0])**2 + (pt0[1] - pt1[1])**2


def perpendicular_distance_2d(pt, pt1, pt2):
    """ return perpendicular distance of pt from the line between pt1 and pt2
    """
    return (((pt2[0] - pt1[0])*(pt1[1] - pt[1])
            - (pt1[0] - pt[0])*(pt2[1] - pt1[1]))
            / sqrt(distance_sqr_2d(pt2, pt1)))


def perpendicular_to_radial(pt, pt2):
    """ return perpendicular distance of pt from the line between the origin
    and pt2
    """
    return (pt[0]*pt2[1] - pt2[0]*pt[1]) / sqrt(pt2[0]**2 + pt2[1]**2)


def perpendicular_to_line(pt, pt1, pt2):
    """ return perpendicular distance of pt from the line between pt1 and pt2
    """
    d = perpendicular_distance_2d(pt, pt1, pt2)
    return (((pt2[0] - pt1[0])*(pt1[1] - pt[1])
            - (pt1[0] - pt[0])*(pt2[1] - pt1[1]))
            / sqrt(distance_sqr_2d(pt2, pt1)))


def perpendicular_from_origin(pt1, pt2):
    """ return perpendicular distance of the origin from the line between
    pt1 and pt2
    """
    return (((pt2[0] - pt1[0])*(pt1[1]) - (pt1[0])*(pt2[1] - pt1[1]))
            / sqrt(distance_sqr_2d(pt2, pt1)))


def projected_point_on_line(pt, pt1, pt2):
    e1 = pt2 - pt1
    e2 = pt - pt1
    # get dot product of e1, e2
    dot_prod = np.dot(e1, e2)
    # get squared length of e1
    len_sqr_e1 = e1[0]**2 + e1[1]**2
    p = np.array([(pt1[0] + (dot_prod * e1[0]) / len_sqr_e1),
                 (pt1[1] + (dot_prod * e1[1]) / len_sqr_e1)])
    return p


def projected_point_on_radial_line(pt, radial_pt):
    # get dot product of radial_pt, pt
    dot_prod = np.dot(radial_pt, pt)
    # get squared length of radial_pt
    len_sqr = radial_pt[0]**2 + radial_pt[1]**2
    p = np.array([(dot_prod * radial_pt[0]) / len_sqr,
                 (dot_prod * radial_pt[1]) / len_sqr])
    return p


def projected_point_on_radial_line_full(pt, radial_pt):
    x_prod = pt[0]*radial_pt[1] - radial_pt[0]*pt[1]
    # get dot product of radial_pt, pt
    dot_prod = np.dot(radial_pt, pt)
    # get squared length of radial_pt
    len_sqr = radial_pt[0]**2 + radial_pt[1]**2
    p = np.array([(dot_prod * radial_pt[0]) / len_sqr,
                 (dot_prod * radial_pt[1]) / len_sqr])
    dist = x_prod/sqrt(len_sqr)
    return p, dist


def euler2opt(e):
    """ convert right-handed euler angles to optical design convention,
        i.e. alpha and beta are left-handed
    """
    return np.array([-e[0], -e[1], e[2]])


def isanumber(a):
    """ returns true if input a can be converted to floating point number """
    try:
        float(a)
        bool_a = True
    except ValueError:
        bool_a = False

    return bool_a


def transpose(mat):
    """ transposes a m x n input list and returns the result """
    mat_t = []
    for j in range(len(mat[0])):
        mat_t.append([mat[i][j] for i in range(len(mat))])
    return mat_t
