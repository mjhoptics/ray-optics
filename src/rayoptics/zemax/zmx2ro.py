#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2020 Michael J. Hayford
"""Post processing functions for Zemax import

.. Created on Mon Aug 10 18:17:55 2020

.. codeauthor: Michael J. Hayford
"""
import numpy as np
from opticalglass import opticalmedium as om
import rayoptics.elem.profiles as profiles
import rayoptics.elem.surface as surface


def apply_fct_to_sm(opt_model, fct, start=None, stop=None, step=None):
    """Iterate in reverse over seq_model.ifcs. Override if needed."""
    sm = opt_model.seq_model
    start = len(sm.ifcs)-1 if start is None else start
    stop = 0 if stop is None else stop
    step = -1 if step is None else step
    num_changes = 0
    for cur in range(start, stop, step):
        if fct(opt_model, cur):
            num_changes += 1
    return num_changes


def convert_to_bend(opt_model, cur):
    """Scan the zemax import for tilted mirrors and convert to BEND types."""
    sm = opt_model.seq_model
    ifc = sm.ifcs[cur]
    if ifc.interact_mode == 'reflect':
        ifc_p = sm.ifcs[cur-1]
        ifc_f = sm.ifcs[cur+1]
        if (ifc_p.z_type == 'COORDBRK' and ifc_f.z_type == 'COORDBRK'):
            if np.array_equal(ifc_f.decenter.euler, ifc_p.decenter.euler):
                ifc.decenter = ifc_p.decenter
                ifc.decenter.dtype = 'bend'
                sm.remove(cur+1, prev=True)
                sm.remove(cur-1)
                return True
    return False


def convert_to_dar(opt_model, cur):
    """Scan the zemax import for tilted surfs and convert to DAR types."""
    sm = opt_model.seq_model
    if cur < len(sm.ifcs)-1:
        ifc = sm.ifcs[cur]
        ifc_p = sm.ifcs[cur-1]
        ifc_f = sm.ifcs[cur+1]
        if (ifc_p.z_type == 'COORDBRK' and ifc_f.z_type == 'COORDBRK'):
            acum_dec = ifc_f.decenter.dec + ifc_p.decenter.dec
            acum_euler = ifc_f.decenter.euler + ifc_p.decenter.euler
            if np.all(acum_dec == 0) and np.all(acum_euler == 0):
                ifc.decenter = ifc_p.decenter
                ifc.decenter.dtype = 'dec and return'
                sm.remove(cur+1, prev=True)
                sm.remove(cur-1)
                return True
    return False


def collapse_coordbrk(opt_model, cur):
    """Attempt to apply the cur COORDBRK to an adjacent real interface."""
    sm = opt_model.seq_model
    ifc_cb = sm.ifcs[cur]
    if ifc_cb.z_type == 'COORDBRK':
        if ifc_cb.decenter.dtype == 'reverse':
            ifc = sm.ifcs[cur-1]
            prev = True
        else:
            ifc = sm.ifcs[cur+1]
            prev = False

        if ifc.decenter is not None:
            return False
        else:
            ifc.decenter = ifc_cb.decenter
            sm.remove(cur, prev=prev)
            return True

    return False


def remove_null_sg(opt_model, cur):
    """Remove sg with planar profile and an adjacent zero thickness air gap."""
    sm = opt_model.seq_model
    ifc = sm.ifcs[cur]
    if is_null_ifc(ifc):
        prev = None
        cur_gap = False if len(sm.gaps)-1 < cur else True
        prev_gap = True if 0 < cur else False
        if cur_gap and is_null_gap(sm.gaps[cur]):
            prev = False
        elif prev_gap and is_null_gap(sm.gaps[cur-1]):
            prev = True
        if prev is not None:
            sm.remove(cur, prev=prev)
            return True

    return False


def is_null_ifc(ifc):
    if isinstance(ifc, surface.Surface):
        if isinstance(ifc.profile, profiles.Spherical):
            if (
                    ifc.profile.cv == 0 and
                    ifc.decenter is None and
                    ifc.interact_mode == 'transmit'
                    ):
                return True
    return False


def is_null_gap(gap):
    if gap.thi == 0 and isinstance(gap.medium, om.Air):
        return True
    else:
        return False
