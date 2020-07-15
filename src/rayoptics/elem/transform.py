#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Useful transforms for processing sequential models

.. Created on Fri Feb  9 10:09:58 2018

.. codeauthor: Michael J. Hayford
"""

import numpy as np


def forward_transform(s1, zdist, s2):
    """ generate transform rotation and translation from
        s1 coords to s2 coords """

    # calculate origin of s2 wrt to s1
    t_orig = np.array([0., 0., zdist])
    r_after_s1 = r_before_s2 = None
    if s1.decenter:
        # get transformation info after s1
        r_after_s1, t_after_s1 = s1.decenter.tform_after_surf()
        t_orig += t_after_s1

    if s2.decenter:
        # get transformation info before s2
        r_before_s2, t_before_s2 = s2.decenter.tform_before_surf()
        t_orig += t_before_s2

    r_cascade = np.identity(3)
    if r_after_s1 is not None:
        # rotate the origin of s2 around s1 "after" transformation
        t_orig = r_after_s1.dot(t_orig)
        r_cascade = r_after_s1
        if r_before_s2 is not None:
            r_cascade = r_after_s1.dot(r_before_s2)
    elif r_before_s2 is not None:
        r_cascade = r_before_s2

    return r_cascade, t_orig


def reverse_transform(s1, zdist, s2):
    """ generate transform rotation and translation from
        s2 coords to s1 coords """

    # calculate origin of s2 wrt to s1
    t_orig = np.array([0., 0., zdist])
    r_after_s1 = r_before_s2 = None
    if s1.decenter:
        # get transformation info after s1
        r_after_s1, t_after_s1 = s1.decenter.tform_after_surf()
        t_orig += t_after_s1

    if s2.decenter:
        # get transformation info before s2
        r_before_s2, t_before_s2 = s2.decenter.tform_before_surf()
        t_orig += t_before_s2

    # going in reverse direction so negate translation
    t_orig = -t_orig

    r_cascade = np.identity(3)
    if r_before_s2 is not None:
        # rotate the origin of s1 around s2 "before" transformation
        r_cascade = r_before_s2.transpose()
        t_orig = r_cascade.dot(t_orig)
        if r_after_s1 is not None:
            r_cascade = r_cascade.dot(r_after_s1.transpose())
    elif r_after_s1 is not None:
        r_cascade = r_after_s1.transpose()

    return r_cascade, t_orig


def cascade_transform(r_prev, t_prev, r_seg, t_seg):
    """ take the seg transform and cascade it with the prev transform """
    return r_prev.dot(r_seg), r_prev.dot(t_seg) + t_prev


def transfer_coords(r_seg, t_seg, pt_s1, dir_s1):
    """ take p and d in s1 coords of seg and transfer them to s2 coords """
    rt = r_seg.transpose()
    return rt.dot(pt_s1 - t_seg), rt.dot(dir_s1)


def transform_before_surface(interface, ray_seg):
    """Transform ray_seg from interface to previous seg.

    Args:
        interface: the :class:'~.Interface' for the path sequence
        ray_seg: ray segment exiting from **interface**

    Returns:
        (**b4_pt**, **b4_dir**)

        - **b4_pt** - ray intersection pt wrt following seg
        - **b4_dir** - ray direction cosine wrt following seg
    """
    if interface.decenter:
        # get transformation info after surf
        r, t = interface.decenter.tform_before_surf()
        if r is None:
            b4_pt, b4_dir = (ray_seg[0] - t), ray_seg[1]
        else:
            rt = r.transpose()
            b4_pt, b4_dir = rt.dot(ray_seg[0] - t), rt.dot(ray_seg[1])
    else:
        b4_pt, b4_dir = ray_seg[0], ray_seg[1]

    return b4_pt, b4_dir


def transform_after_surface(interface, ray_seg):
    """Transform ray_seg from interface to following seg.

    Args:
        interface: the :class:'~.Interface' for the path sequence
        ray_seg: ray segment exiting from **interface**

    Returns:
        (**b4_pt**, **b4_dir**)

        - **b4_pt** - ray intersection pt wrt following seg
        - **b4_dir** - ray direction cosine wrt following seg
    """
    if interface.decenter:
        # get transformation info after surf
        r, t = interface.decenter.tform_after_surf()
        if r is None:
            b4_pt, b4_dir = (ray_seg[0] - t), ray_seg[1]
        else:
            rt = r.transpose()
            b4_pt, b4_dir = rt.dot(ray_seg[0] - t), rt.dot(ray_seg[1])
    else:
        b4_pt, b4_dir = ray_seg[0], ray_seg[1]

    return b4_pt, b4_dir
