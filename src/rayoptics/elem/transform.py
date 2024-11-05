#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Useful transforms for processing sequential models

.. Created on Fri Feb  9 10:09:58 2018

.. codeauthor: Michael J. Hayford
"""

import numpy as np
import itertools


def compute_global_coords(seq_model, glo=1, origin=None):
    """ Return global surface coordinates (rot, t) wrt surface `glo`. 
    
    Args:
        seq_model: sequential model
        glo: global reference surface index
        origin: None | (r_origin, t_origin)

    The tuple (r_origin, t_origin) is the transform from the desired 
    global origin to the global surface `glo`.
    """
    def accumulate_transforms(seq, b4_seg, transform_calc, 
                              tfrm_prev, tfrm_dir: int):
        b4_ifc, b4_gap, b4_z_dir = b4_seg
        r_prev, t_prev = tfrm_prev
        for (ifc, gap, z_dir) in seq:
            zdist = tfrm_dir * b4_gap.thi
            r, t = transform_calc(b4_ifc, zdist, ifc)

            t_new = np.matmul(r_prev, t) + t_prev
            r_new = np.matmul(r_prev, r)

            tfrms.append((r_new, t_new))
            r_prev, t_prev = r_new, t_new
            b4_ifc, b4_gap, b4_z_dir = ifc, gap, z_dir

    # Initialize origin of global coordinate system.
    tfrms = []
    if origin is None:
        r_origin, t_origin = np.identity(3), np.array([0., 0., 0.])
    else:
        r_origin, t_origin = origin
    tfrm_origin = r_origin, t_origin
    tfrms.append(tfrm_origin)

    # Compute transforms from global surface to object surface
    if glo > 0:
        # iterate in reverse over the segments before the
        #  global reference surface
        step = -1
        seq = itertools.zip_longest(seq_model.ifcs[glo::step],
                                    seq_model.gaps[glo-1::step],
                                    seq_model.z_dir[glo-1::step])
        b4_seg = next(seq)
        # loop of remaining surfaces in path
        accumulate_transforms(seq, b4_seg, reverse_transform, 
                              tfrm_origin, -1)
        tfrms.reverse()

    # Compute transforms from global surface to image surface
    seq = itertools.zip_longest(seq_model.ifcs[glo:], 
                                seq_model.gaps[glo:], 
                                seq_model.z_dir[glo:])
    b4_seg = next(seq)
    accumulate_transforms(seq, b4_seg, forward_transform, 
                          tfrm_origin, +1)

    return tfrms


def compute_local_transforms(seq_model, seq, step):
    """ Return forward surface coordinates (r.T, t) for each interface. """
    def local_transform(seq, transform_calc, tfrm_dir: int):
        b4_ifc, b4_gap = next(seq)
        for (ifc, gap) in seq:
            zdist = tfrm_dir * b4_gap.thi
            r, t = transform_calc(b4_ifc, zdist, ifc)
            rt = r.transpose()
            tfrms.append((rt, t))
            b4_ifc, b4_gap = ifc, gap
        return tfrms

    num_ifcs = seq_model.get_num_surfaces()
    tfrms = []
    if step == -1:
        if seq is None:
            seq = itertools.zip_longest(seq_model.ifcs[num_ifcs::step], 
                                        seq_model.gaps[num_ifcs-1::step])
        tfrms = local_transform(seq, reverse_transform, -1)
    elif step == 1:
        if seq is None:
            seq = itertools.zip_longest(seq_model.ifcs[::step],
                                        seq_model.gaps[::step])
        tfrms = local_transform(seq, forward_transform, 1)
    else:
        tfrms = []

    tfrms.append((np.identity(3), np.array([0., 0., 0.])))
    return tfrms


def forward_transform(s1, zdist, s2):
    """ generate transform rotation and translation from
        s1 coords to s2 coords """

    t_orig = np.array([0., 0., zdist])
    r_after_s1 = r_before_s2 = None
    if s1.decenter:
        r_after_s1, t_after_s1 = s1.decenter.tform_after_surf()
        t_orig += t_after_s1

    if s2.decenter:
        r_before_s2, t_before_s2 = s2.decenter.tform_before_surf()
        t_orig += t_before_s2

    r_cascade = np.identity(3)
    if r_after_s1 is not None:
        t_orig = np.matmul(r_after_s1, t_orig)
        r_cascade = r_after_s1
        if r_before_s2 is not None:
            r_cascade = np.matmul(r_after_s1, r_before_s2)
    elif r_before_s2 is not None:
        r_cascade = r_before_s2

    return r_cascade, t_orig


def reverse_transform(s2, zdist, s1):
    """ generate transform rotation and translation from
        s2 coords to s1 coords, applying transforms in the reverse order """
    t_orig = np.array([0., 0., zdist])
    r_before_s2 = r_after_s1 = None
    if s2.decenter:
        r_before_s2, t_before_s2 = s2.decenter.tform_before_surf()
        t_orig += t_before_s2

    if s1.decenter:
        r_after_s1, t_after_s1 = s1.decenter.tform_after_surf()
        t_orig += t_after_s1

    r_cascade = np.identity(3)
    if r_before_s2 is not None:
        r_cascade = r_before_s2.transpose()
        t_orig = np.matmul(r_cascade, t_orig)
        if r_after_s1 is not None:
            r_cascade = np.matmul(r_cascade, r_after_s1.transpose())
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
        interface: the :class:`~.interface.Interface` for the path sequence
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
        interface: the :class:`~.interface.Interface` for the path sequence
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
