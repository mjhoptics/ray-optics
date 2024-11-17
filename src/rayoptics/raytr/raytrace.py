#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Functions to support ray tracing a sequential optical model

.. Created on Thu Jan 25 11:01:04 2018

.. codeauthor: Michael J. Hayford
"""

import numpy as np
from numpy.linalg import norm
from math import sqrt, copysign

import rayoptics.optical.model_constants as mc
from .traceerror import (TraceMissedSurfaceError, TraceTIRError,
                         TraceRayBlockedError, TraceEvanescentRayError)

def bend(d_in, normal, n_in, n_out):
    """ refract incoming direction, d_in, about normal """
    try:
        normal_len = norm(normal)
        cosI = np.dot(d_in, normal)/normal_len
        sinI_sqr = 1.0 - cosI*cosI
        n_cosIp = copysign(sqrt(n_out*n_out - n_in*n_in*sinI_sqr), cosI)
        alpha = n_cosIp - n_in*cosI
        d_out = (n_in*d_in + alpha*normal)/n_out
        return d_out
    except ValueError:
        raise TraceTIRError(d_in, normal, n_in, n_out)


def reflect(d_in, normal):
    """ reflect incoming direction, d_in, about normal """
    normal_len = norm(normal)
    cosI = np.dot(d_in, normal)/normal_len
    d_out = d_in - 2.0*cosI*normal
    return d_out


def phase(ifc, inc_pt, d_in, normal, ifc_cntxt):
    """ apply phase shift to incoming direction, d_in, about normal """
    try:
        d_out, dW = ifc.phase(inc_pt, d_in, normal, ifc_cntxt)
        return d_out, dW
    except ValueError:
        z_dir, wvl, n_in, n_out, interact_mode = ifc_cntxt
        raise TraceEvanescentRayError(ifc, inc_pt, d_in, normal, n_in, n_out)


def trace(seq_model, pt0, dir0, wvl, **kwargs):
    """ fundamental raytrace function

    Args:
        seq_model: the sequential model to be traced
        pt0: starting point in coords of first interface
        dir0: starting direction cosines in coords of first interface
        wvl: wavelength in nm

    Returns:
        (**ray**, **op_delta**, **wvl**)

        - **ray** is a list for each interface in **path_pkg** of these
          elements: [pt, after_dir, after_dst, normal]

            - pt: the intersection point of the ray
            - after_dir: the ray direction cosine following the interface
            - after_dst: after_dst: the geometric distance to the next
              interface
            - normal: the surface normal at the intersection point

        - **op_delta** - optical path wrt equally inclined chords to the
          optical axis
        - **wvl** - wavelength (in nm) that the ray was traced in
    """
    path = seq_model.path(wvl)
    kwargs['first_surf'] = kwargs.get('first_surf', 1)
    kwargs['last_surf'] = kwargs.get('last_surf',
                                     seq_model.get_num_surfaces()-2)
    return trace_raw(path, pt0, dir0, wvl, **kwargs)


def trace_raw(path, pt0, dir0, wvl, eps=1.0e-12, check_apertures=False, 
              intersect_obj=True, **kwargs):
    """ fundamental raytrace function

    Args:
        path: an iterator containing interfaces and gaps to be traced.
              for each iteration, the sequence or generator should return a
              list containing: **Intfc, Gap, Trfm, Index, Z_Dir**
        pt0: starting point in coords of object interface
        dir0: starting direction cosines in coords of object interface
        wvl: wavelength in nm
        eps: accuracy tolerance for surface intersection calculation
        check_apertures: if True, do point_inside() test on inc_pt
        intersect_obj: if True, intersect the ray with the object, otherwise 
                       trace input ray coords directly.
        pt_inside_fuzz: accuracy tolerance for aperture clipping check

    Returns:
        (**ray**, **op_delta**, **wvl**)

        - **ray** is a list for each interface in **path** of these
          elements: [pt, after_dir, after_dst, normal]

            - pt: the intersection point of the ray
            - after_dir: the ray direction cosine following the interface
            - after_dst: the geometric distance to the next interface
            - normal: the surface normal at the intersection point

        - **op_delta** - optical path wrt equally inclined chords to the
          optical axis
        - **wvl** - wavelength (in nm) that the ray was traced in
    """
    ray = []

    first_surf = kwargs.get('first_surf', 0)
    last_surf = kwargs.get('last_surf', None)
    pt_inside_fuzz = kwargs.get('pt_inside_fuzz', None)
    fuzz = {} if pt_inside_fuzz is None else {'fuzz': pt_inside_fuzz}

    def in_gap_range(gap_indx, include_last_surf=False):
        if first_surf == last_surf:
            return False
        if gap_indx < first_surf:
            return False
        if last_surf is None:
            return True
        else:
            return (gap_indx <= last_surf if include_last_surf 
                    else gap_indx < last_surf)

    def in_surface_range(s):
        if s < first_surf:
            return False
        if last_surf is None:
            return True
        elif s > last_surf:
            return False
        else:
            return True

    # trace object surface
    before = obj = next(path)
    if intersect_obj:
        srf_obj = obj[mc.Intfc]
        _, before_pt = srf_obj.intersect(pt0, dir0, z_dir=obj[mc.Zdir])
        before_normal = srf_obj.normal(before_pt)
    else:
        before_pt = pt0
        before_normal = np.array([0., 0., 1.])

    before_dir = dir0
    tfrm_from_before = before[mc.Tfrm]
    z_dir_before = before[mc.Zdir]

    op_delta = 0.0
    opl = 0.0
    surf = 0
    # loop of remaining surfaces in path
    while True:
        try:
            after = next(path)
            surf += 1

            # transform ray data from previous ifc coords to current one
            rt, t = tfrm_from_before
            b4_pt, b4_dir = rt.dot(before_pt - t), rt.dot(before_dir)

            pp_dst = -b4_pt.dot(b4_dir)
            pp_pt_before = b4_pt + pp_dst*b4_dir

            ifc = after[mc.Intfc]
            z_dir_after = after[mc.Zdir]

            # intersect ray with profile
            pp_dst_intrsct, inc_pt = ifc.intersect(pp_pt_before, b4_dir,
                                                   eps=eps, z_dir=z_dir_before)
            dst_b4 = pp_dst + pp_dst_intrsct
            
            # add *previous* intersection point, direction, etc., to ray
            ray.append([before_pt, before_dir, dst_b4, before_normal])

            if in_gap_range(surf-1):
                opl += before[mc.Indx] * dst_b4

            normal = ifc.normal(inc_pt)

            if check_apertures and in_surface_range(surf):
                if not ifc.point_inside(inc_pt[0], inc_pt[1], **fuzz):
                    raise TraceRayBlockedError(ifc, inc_pt)

            # if present, use the phase element to calculate after_dir
            if hasattr(ifc, 'phase_element'):
                ifc_cntxt = (z_dir_before, wvl, 
                             before[mc.Indx], after[mc.Indx],
                             ifc.interact_mode)
                after_dir, phs = phase(ifc, inc_pt, b4_dir, normal, ifc_cntxt)
                op_delta += phs
            else:  # refract or reflect ray at interface
                if ifc.interact_mode == 'reflect':
                    after_dir = reflect(b4_dir, normal)
                elif ifc.interact_mode == 'transmit':
                    after_dir = bend(b4_dir, normal, before[mc.Indx], after[mc.Indx])
                elif ifc.interact_mode == 'dummy':
                    after_dir = b4_dir
                else:  # no action, input becomes output
                    after_dir = b4_dir

            before_pt = inc_pt
            before_normal = normal
            before_dir = after_dir
            z_dir_before = z_dir_after
            before = after
            tfrm_from_before = before[mc.Tfrm]

        except TraceMissedSurfaceError as ray_miss:
            ray.append([before_pt, before_dir, pp_dst, before_normal])
            ray_miss.surf = surf
            ray_miss.ifc = ifc
            ray_miss.prev_tfrm = before[mc.Tfrm]
            ray_miss.ray_pkg = ray, opl, wvl
            raise ray_miss

        except TraceTIRError as ray_tir:
            ray.append([inc_pt, before_dir, 0.0, normal])
            ray_tir.surf = surf
            ray_tir.ifc = ifc
            ray_tir.int_pt = inc_pt
            ray_tir.ray_pkg = ray, opl, wvl
            raise ray_tir

        except TraceRayBlockedError as ray_blocked:
            ray.append([inc_pt, before_dir, 0.0, normal])
            ray_blocked.surf = surf
            ray_blocked.ray_pkg = ray, opl, wvl
            raise ray_blocked
            
        except TraceEvanescentRayError as ray_evn:
            ray.append([inc_pt, before_dir, 0.0, normal])
            ray_evn.surf = surf
            ray_evn.ray_pkg = ray, opl, wvl
            raise ray_evn

        except StopIteration:
            ray.append([inc_pt, after_dir, 0.0, normal])
            op_delta += opl
            break

    return ray, op_delta, wvl


def calc_optical_path(ray, path):
    """ computes equally inclined chords and path info for ray

    Args:
        ray: ray data for traced ray
        path: an iterator containing interfaces and gaps to be traced.
              for each iteration, the sequence or generator should return a
              list containing: **Intfc, Gap, Trfm, Index, Z_Dir**

    Returns:
        **ray_op**

        - **ray_op** - the optical path between the first and last optical
          surfaces
    """
    num_items = len(ray)
    ray_seq_iter = zip(ray, path)
    next(ray_seq_iter)
    ray_op = 0
    for i in range(1, num_items-2):
        after_ray_seg, surf = next(ray_seq_iter)

        n_after = surf[mc.Indx]
        dist = after_ray_seg[mc.dst]
        ray_op += n_after*dist

    return ray_op
