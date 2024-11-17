#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" miscellaneous functions for working with numpy vectors and floats

.. Created on Wed May 23 15:27:06 2018

.. codeauthor: Michael J. Hayford
"""
import numpy as np
from numpy.linalg import norm
from math import sqrt, pi, acos
import transforms3d as t3d


def infinity_guard(x: float, big: float = 1e12) -> float:
    """ Replace IEEE inf with a signed big number. """
    return (-big if np.isneginf(x) else big) if np.isinf(x) else x


def is_kinda_big(x: float, kinda_big: float = 1e8) -> bool:
    """ Test for IEEE inf as well as any \|x| > kinda_big  """
    if np.isinf(x):
        return True
    elif np.abs(x) > kinda_big:
        return True
    else:
        return False


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


def rot_v1_into_v2(v1, v2):
    """ rotate v1 into v2 using equivalent angle rotation. 
    
    Compute a rotation matrix from v1 to v2. Take the cross product of the input vectors to get the rotation axis. The eqivalent angle rotation is equation 2.80 from Introduction to Robotics, 2nd ed, by John J Craig.
    """
    rot_axis = -np.cross(v1, v2)
    s = np.linalg.norm(rot_axis)
    c = cosine_ang = np.dot(v1, v2)
    v = 1 - cosine_ang
    ax = normalize(rot_axis)
    rot_mat = np.array(
        [[ax[0]*ax[0]*v + c, ax[0]*ax[1]*v - ax[2]*s, ax[0]*ax[2]*v + ax[1]*s],
         [ax[0]*ax[1]*v + ax[2]*s, ax[1]*ax[1]*v + c, ax[1]*ax[2]*v + ax[0]*s],
         [ax[0]*ax[2]*v + ax[1]*s, ax[1]*ax[2]*v + ax[0]*s, ax[2]*ax[2]*v + c]]
         )
    return rot_mat


def euler2opt(e):
    """ convert right-handed euler angles to optical design convention,
        i.e. alpha and beta are left-handed
    """
    return np.array([-e[0], -e[1], e[2]])


def euler2rot3d(euler):
    """ convert euler angle vector to a rotation matrix. """
    rot_mat = t3d.euler.euler2mat(*np.deg2rad(euler2opt(euler)))
    return rot_mat


def isanumber(a):
    """ returns true if input a can be converted to floating point number """
    try:
        float(a)
        bool_a = True
    except ValueError:
        bool_a = False
    except TypeError:
        bool_a = False

    return bool_a


def transpose(mat):
    """ transposes a m x n input list and returns the result """
    mat_t = []
    for j in range(len(mat[0])):
        mat_t.append([mat[i][j] for i in range(len(mat))])
    return mat_t


def circle_intersection_area(ra, rb, d):
    """ return the area of the intersection of 2 circles

    Args:
        ra: radius of first circle
        rb: radius of second circle
        d: separation of the circles' centers of curvature

    Returns:
        area of the circle intersection

    `Weisstein, Eric W. "Circle-Circle Intersection." From MathWorld--A Wolfram Web
    Resource. <http://mathworld.wolfram.com/Circle-CircleIntersection.html>`_
    """
    if ra < rb:  # sort into ascending order for convenience
        r, R = ra, rb
    else:
        r, R = rb, ra

    r2 = r**2
    if R >= r + d:  # smaller circle contained inside the larger one
        return pi*r2

    if d > r + R:  # circles are completely separated - no overlap
        return 0

    R2 = R**2
    d2 = d**2

    # calculate area via eq 14 in the referenced link
    p1 = r2*acos((d2 + r2 - R2)/(2*d*r))
    p2 = R2*acos((d2 + R2 - r2)/(2*d*R))
    p3 = sqrt((-d + r + R)*(d + r - R)*(d - r + R)*(d + r + R))/2
    area = p1 + p2 - p3
    return area


def compute_tangent_point_to_circle(CofC, r, pt):
    """ return the area of the intersection of 2 circles

    Args:
        CofC: center of curvature of circle (2d numpy array)
        r: radius of circle
        pt: 2d numpy array of point outside of circle

    Returns:
        the 2 tangent points for lines from pt to circle

    `gboffi. <https://math.stackexchange.com/users/467357/gboffi>`_ 
    `"How to find the equation of a line, tangent to a circle, that passes 
    through a given external point." StackExchange (version: 2019-05-30) 
    <https://math.stackexchange.com/a/3190374>`_
    """
    dxdy = pt - CofC
    dxdyr = np.array([-dxdy[1], dxdy[0]])
    d = sqrt(dxdy[0]**2 + dxdy[1]**2)
    if d > r:
        rho = r/d
        ad = rho**2
        bd = rho * sqrt(1 - rho**2)
        T1 = CofC + ad*dxdy + bd*dxdyr
        T2 = CofC + ad*dxdy - bd*dxdyr
        return T1, T2
