#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2024 Michael J. Hayford
""" Wide angle raytrace and ray aiming

This package was developed to ray trace fisheye, i.e. very wide field, lenses. These lenses have significant pupil spherical aberration. In order to trace highly oblique field angles, one must locate the actual entrance pupil location for a field angle. The function find_real_enp implements the search by parameterizing the z offset from the 1st interface along the object space optical axis. For extreme field angles, the rays that successfully reach the stop surface are far from the z offset for the paraxial entrance pupil. 

The function z_enp_coordinate is used to evaluate where the chief ray hits the stop surface. This function is evaluated at regular intervals, spaced between the z distance of the paraxial entrance pupil and the first surface vertex (i.e. z_enp = 0). When the range of z values is identified that pass rays through the complete system, find_z_enp is called to find the exact conjugate point to the center of the stop surface.

The entrance pupil for the wide angle package is taken as normal to the chief ray at that field angle. The set_vignetting function using bisection works well with this definition.

.. Created on Mon Oct 14 11:17 2024

.. codeauthor: Michael J. Hayford
"""
import warnings
import logging

import numpy as np
from scipy.optimize import newton, fsolve, bisect

import rayoptics.raytr.raytrace as rt
from rayoptics.raytr import trace
from rayoptics.raytr import RayResult, RayPkg
from rayoptics.raytr.traceerror import TraceError, TraceMissedSurfaceError
from rayoptics.util.misc_math import normalize, rot_v1_into_v2
import rayoptics.optical.model_constants as mc

logger = logging.getLogger(__name__)


def enp_z_coordinate(z_enp, *args):
    """ Trace a ray thru the center of the entrance pupil at z_enp. 
    
    Args:

        z_enp:      The z distance from the 1st interface of the 
                    entrance pupil for the field angle dir0
        seq_model:  The sequential model
        stop_idx:   index of the aperture stop interface
        dir0:       direction cosine vector in object space
        obj_dist:   object distance to first interface
        wvl:        wavelength of raytrace (nm)

    """
    seq_model, stop_idx, dir0, obj_dist, wvl = args
    obj2enp_dist = -(obj_dist + z_enp)
    pt1 = np.array([0., 0., obj2enp_dist])
    rot_mat = rot_v1_into_v2(np.array([0., 0., 1.]), dir0)
    pt0 = np.matmul(rot_mat, pt1) - pt1
    # handle case where entrance pupil is behind the object (issue #120)
    # if dir0[2] * seq_model.z_dir[0] < 0:
    #     dir0 = -dir0

    try:
        ray_pkg = RayPkg(*rt.trace(seq_model, pt0, dir0, wvl, 
                                   intersect_obj=False))

    except TraceError as ray_error:
        logger.debug(f'ray_error: "{type(ray_error).__name__}", '
                     f'{ray_error.surf=}')
        ray_pkg = RayPkg(*ray_error.ray_pkg)
        rr = RayResult(ray_pkg, ray_error)
        final_coord = np.array([0., 0., 0.])

    else:
        rr = RayResult(ray_pkg, None)
        final_coord = ray_pkg.ray[stop_idx][mc.p]

    return final_coord, rr


def find_z_enp(opt_model, stop_idx, z_enp_0, fld, wvl, **kwargs):
    """ iterates a ray to [0, 0] on interface stop_ifc, returning aim info
    
    Args:

        opt_model:  input OpticalModel
        stop_idx:   index of the aperture stop interface
        z_enp_0:    estimate of pupil location. this estimate must support 
                    a raytrace up to stop_ifc
        fld:        field point
        wvl:        wavelength of raytrace (nm)

    Returns z distance from 1st interface to the entrance pupil.

    If stop_ifc is None, i.e. a floating stop surface, returns paraxial entrance pupil.

    If the iteration fails, a TraceError will be raised
    """
    rr = None
    def eval_z_enp(z_enp, *args):
        nonlocal rr
        y_target = args[-1]
        final_coord, rr = enp_z_coordinate(z_enp, *args[:-1])
        # print(f"{final_coord}")
        return final_coord[1] - y_target
        
    seq_model = opt_model['seq_model']
    osp = opt_model['optical_spec']

    fod = opt_model['analysis_results']['parax_data'].fod
    z_enp = z_enp_0
    obj_dist = fod.obj_dist

    pt0, dir0 = osp.obj_coords(fld)
    der = 0.005
    y_target = 0.  # chief ray -> center of stop surface
    results = None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if stop_idx is not None:
            # do 1D iteration if field and target points are zero in x
            try:
                z_enp, results = newton(eval_z_enp, z_enp,
                                        args=(seq_model, stop_idx, dir0,
                                                obj_dist, wvl, y_target),
                                        x1=(1+der)*z_enp,
                                        disp=False, full_output=True)
            except RuntimeError as rte:
                # if we come here, start_y is a RuntimeResults object
                # print(rte)
                z_enp = results.root
            except TraceError as ray_err:
                logger.debug(f"trace error: {ray_err.surf}")
                z_enp = 0.0
            start_coords = np.array([0., 0., z_enp])

        else:  # floating stop surface - use entrance pupil for aiming
            start_coords = np.array([0., 0., fod.enp_dist])

    return start_coords, rr, results


def find_real_enp(opm, stop_idx, fld, wvl):
    """ Locate the z center of the real pupil for `fld`, wrt 1st ifc
    
    This function implements a 2 step process to finding the chief ray for `fld` and `wvl`.
    The first phase searches for the window of pupil locations by sampling the z coordinate from the paraxial pupil location towards the first interface vertex. Failed rays are discarded until a range of z coordinates is found where rays trace successfully. If only a single successful trace is in hand, a second, more finely subdivided search is conducted about the successful point.
    The outcome is a range, start_z -> end_z, that is divided in 3 and a ray iteration to find the center of the stop surface is done. Sometimes the start point doesn't produce a solution; use of the mid-point as a start is a reliable second try.
    """
    sm = opm['seq_model']
    osp = opm['osp']

    fod = opm['ar']['parax_data'].fod

    pt0, dir0 = osp.obj_coords(fld)
    #print(f"field angle={fld_ang}:   sine={dir0[1]:8.4f}")

    # If there is aim_info, try it and return if good.
    if fld.aim_info is not None:
        z_enp = fld.aim_info
        args = sm, stop_idx, dir0, fod.obj_dist, wvl
        final_coord, rr = enp_z_coordinate(z_enp, *args)
        tol = 1.48e-08
        if abs(final_coord[1])<tol:
            return z_enp, rr

    # filter on-axis chief ray. z_enp is the paraxial result.
    z_enp_0 = fod.enp_dist
    if dir0[2] == 1:  # axial chief ray
        args = sm, stop_idx, dir0, fod.obj_dist, wvl
        final_coord, rr = enp_z_coordinate(z_enp_0, *args)
        # print(f"axial chief {z_enp_0=:8.4f}  {rr.err is None}")
        return z_enp_0, rr

    start_z = None
    end_z = None
    del_z = -z_enp_0/16
    z_enp = z_enp_0
    keep_going = True
    # protect against infinite loops
    trial = 0
    # if the trace succeeds 5 times in a row, go on to the next phase
    successes = 0
    while keep_going and successes < 4 and trial < 64:
        args = sm, stop_idx, dir0, fod.obj_dist, wvl
        final_coord, rr = enp_z_coordinate(z_enp, *args)
        if rr.err is None:
            successes += 1
            if start_z is None:
                start_z = z_enp
            end_z = z_enp
 
        elif isinstance(rr.err, TraceMissedSurfaceError):
            # if the first surface was missed, then exit
            if start_z is not None:
                keep_going = False
        z_enp += del_z
        trial += 1

    # print(f"trials: {trial}")

    # If start and end are equal, then only one ray was successful.
    # Sample z_enp evenly 1 del_z to either side.
    if start_z == end_z:
        start_new = start_z - del_z
        end_new = end_z + del_z
        start_z = None
        end_z = None
        for z_enp in np.linspace(start_new, end_new, num=8):
            args = sm, stop_idx, dir0, fod.obj_dist, wvl
            final_coord, rr = enp_z_coordinate(z_enp, *args)
            if rr.err is None:
                if start_z is None:
                    start_z = z_enp
                end_z = z_enp
            # print(f"singular point {z_enp=:8.4f}  {rr.err is None}")

    # Now that candidate z_enps have been identified that trace without
    # ray failures, iterate to find the ray thru the stop center
    starting_pts = [start_z, (start_z + end_z)/2, end_z]
    for init_z in starting_pts:
        start_coords, rr, results = find_z_enp(opm, stop_idx, init_z, 
                                               fld, wvl)
        if rr.err is None:
            # print(f"iter start {init_z:8.4f},  z_enp {start_coords[2]:8.4f}")
            break
    z_enp = start_coords[2]

    final_coord = rr.pkg.ray[stop_idx][mc.p]
    y_ht = final_coord[1]
    # print(f"fld: {fld.y:3.1f}:   {start_z:8.4f}  {end_z:8.4f}  "
    #              f"{z_enp:8.4f}  {y_ht:10.2e}")
    logger.debug(f"{fld.y:3.1f}:   {start_z:8.4f}  {end_z:8.4f}  "
                 f"{z_enp:8.4f}  {y_ht:10.2e}")
    return z_enp, rr


def eval_real_image_ht(opt_model, fld, wvl):
    """Trace reverse ray from image point to get object space inputs. 
    
    This the implementation of obj_coords() for ('image', 'real height')
    """
    sm = opt_model['seq_model']
    osp = opt_model['optical_spec']
    fov = osp['fov']
    fod = opt_model['analysis_results']['parax_data'].fod

    rpath = sm.reverse_path(wl=wvl, start=len(sm.ifcs), stop=None, step=-1)
    rpath_list = list(rpath)
    obj2pup_dist = fod.exp_dist - fod.img_dist
    p_i = np.array([fld.x, fld.y, 0])
    not_wa = not fov.is_wide_angle
    if fov.is_relative:
        p_i *= fov.value
    p_exp = np.array([0, 0, obj2pup_dist])
    d_i = normalize(p_exp - p_i)
    xy_target = [0., 0.]
    eprad = fod.exp_radius
    ifcx = len(sm.ifcs) - sm.stop_surface - 1
    start_coords, rrev_cr = trace.iterate_ray_raw(rpath_list, ifcx, xy_target, 
                                                  p_i, d_i, obj2pup_dist, 
                                                  eprad, wvl, not_wa)
    p_k = rrev_cr.pkg.ray[-2][mc.p]
    d_k = rrev_cr.pkg.ray[-2][mc.d]
    d_o = -d_k
    d_k01 = np.sqrt(d_k[0]**2 + d_k[1]**2)
    if d_k01 == 0.:
        z_enp = fod.enp_dist
    else:
        z_enp = p_k[2] + p_k[1]*d_k[2]/d_k01

    obj2enp_dist = fod.obj_dist + z_enp
    enp_pt = np.array([0., 0., obj2enp_dist])
    p_o = enp_pt + obj2enp_dist * d_k
    return (p_o, d_o), z_enp
