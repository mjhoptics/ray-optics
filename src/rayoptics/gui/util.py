#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" utility functions for gui functions

.. Created on Wed Oct 16 15:27:40 2019

.. codeauthor: Michael J. Hayford
"""

from collections import namedtuple
import numpy as np

import opticalglass.glasspolygons as gp  # type: ignore

GUIHandle = namedtuple('GUIHandle', ['poly', 'bbox'])
GUIHandle.poly.__doc__ = "poly entity for underlying graphics system (e.g. mpl)"
GUIHandle.bbox.__doc__ = "bounding box for poly"


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
    if len(np.array(poly).shape) > 1:
        minx, miny = np.min(poly, axis=0)
        maxx, maxy = np.max(poly, axis=0)
        bbox = np.array([[minx, miny], [maxx, maxy]])
    else:
        x = poly[0]
        y = poly[1]
        bbox = np.array([[x, y], [x, y]])
    return bbox


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
    return poly


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


def fit_data_range(x_data, margin=0.05, range_trunc=0.25, **kwargs):
    x_min = min(0., min(x_data))
    x_max = max(0., max(x_data))
    x_range = x_max - x_min
    if x_range != 0.0 and len(x_data) > 2:
        x1_min = min(0., min(x_data[1:]))
        x1_max = max(0., max(x_data[1:]))
        x1_range = x1_max - x1_min
        if abs(x1_range/x_range) < range_trunc:
            x_min = x1_min
            x_max = x1_max
            x_range = x1_range

    if x_range > 0.:
        x_margin = margin*x_range
    else:
        x_margin = 0.01
    return x_min-x_margin, x_max+x_margin


def calc_render_color_for_material(matl):
    """ get element color based on V-number of glass"""
    try:
        gc = float(matl.glass_code())
    except AttributeError:
        return (255, 255, 255, 64)  # white
    else:
        # set element color based on V-number
        indx = round(1.0 + (int(gc)/1000), 3)
        vnbr = round(100.0*(gc - int(gc)), 3)
        dsg, rgb = gp.find_glass_designation(indx, vnbr)
        if rgb is None:
            return [228, 237, 243, 64]  # ED designation
        return rgb
