#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" utility functions for gui functions

.. Created on Wed Oct 16 15:27:40 2019

.. codeauthor: Michael J. Hayford
"""


from collections import namedtuple
import numpy as np

GUIHandle = namedtuple('GUIHandle', ['poly', 'bbox'])
""" tuple grouping together graphics entity and bounding box

    Attributes:
        poly: poly entity for underlying graphics system (e.g. mpl)
        bbox: bounding box for poly
"""


def transform_ray_seg(poly, r, tfrm):
    rot, trns = tfrm
    p = rot.dot(r.p) + trns
    poly.append([p[2], p[1]])


def bbox_from_poly(poly):
    minx, miny = np.min(poly, axis=0)
    maxx, maxy = np.max(poly, axis=0)
    return np.array([[minx, miny], [maxx, maxy]])


def scale_bounds(bbox, oversize_factor):
    inc_x = oversize_factor*(bbox[1][0] - bbox[0][0])
    inc_y = oversize_factor*(bbox[1][1] - bbox[0][1])
    incr = max(inc_x, inc_y)
    return np.array([[bbox[0][0]-incr, bbox[0][1]-incr],
                     [bbox[1][0]+incr, bbox[1][1]+incr]])


def transform_poly(tfrm, poly):
    coord_flip = np.array([[0., 1.], [1., 0.]])

    poly = np.matmul(coord_flip, poly.T)
    poly = np.matmul(tfrm[0][1:, 1:], poly).T

    t = np.array([tfrm[1][1], tfrm[1][2]])
    poly += t

    # flip coordinates back to 2D plot coordinates, +y points up
    poly = np.matmul(poly, coord_flip)
    bbox = bbox_from_poly(poly)
    return poly, bbox


def inv_transform_poly(tfrm, poly):
    coord_flip = np.array([[0., 1.], [1., 0.]])
    try:
        poly = np.matmul(coord_flip, poly.T)
    except TypeError:
        print(poly)

    t = np.array([tfrm[1][1], tfrm[1][2]])
    poly -= t

    poly = np.matmul(tfrm[0][1:, 1:], poly).T

    # flip coordinates back to 2D plot coordinates, +y points up
    poly = np.matmul(poly, coord_flip)
    return poly
