#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2022 Michael J. Hayford
""" Vignetting and clear aperture setting operations

.. Created on Mon Apr 18 15:28:25 2022

.. codeauthor: Michael J. Hayford
"""
import logging

from typing import Optional, Sequence

from math import sqrt, copysign
import numpy as np
from numpy import sqrt
from scipy.optimize import newton

import rayoptics.optical.model_constants as mc

from rayoptics.raytr import trace, RayPkg, RaySeg
from rayoptics.raytr import traceerror as terr
from rayoptics.parax import etendue

logger = logging.getLogger(__name__)

# label for coordinate chooser
xy_str = 'xy'


def max_aperture_at_surf(rayset, i):   
    max_ap = -1.0e+10
    for f in rayset:
        for p in f:
            ray = p.ray
            if len(ray) > i:
                ap = sqrt(ray[i].p[0]**2 + ray[i].p[1]**2)
                if ap > max_ap:
                    max_ap = ap
            else:  # ray failed before this interface, don't update
                return None
    return max_ap


def set_clear_apertures(opt_model: 'OpticalModel', 
                        avoid_list: Optional[Sequence[int]]=None, 
                        include_list: Optional[Sequence[int]]=None):
    """ From existing fields and vignetting, calculate clear apertures. 

    Args:
        avoid_list: list of surfaces to skip when setting apertures.
        include_list: list of surfaces to include when setting apertures.

    If specified, only one of either `avoid_list` or `include_list` should be specified. If neither is specified, all surfaces are set. If both are specified, the `avoid_list` is used.

    If a surface is specified as the aperture stop, that surface's aperture is determined from the boundary rays of the first field.
    
    The avoid_list idea and implementation was contributed by Quentin Bécar
    """
    sm = opt_model['seq_model']
    num_surfs = sm.get_num_surfaces()
    if avoid_list is None:
        if include_list is None:
            include_list = range(num_surfs)
    else:
        include_list = [i for i in range(num_surfs) if i not in avoid_list]
    
    rayset = trace.trace_boundary_rays(opt_model, use_named_tuples=True)

    stop_surf = sm.stop_surface
    if stop_surf is not None and stop_surf in include_list:
            max_ap = max_aperture_at_surf([rayset[0]], stop_surf)
            if max_ap is not None:
                sm.ifcs[stop_surf].set_max_aperture(max_ap)

    for i in include_list:
        if i != stop_surf:
            max_ap = max_aperture_at_surf(rayset, i)
            if max_ap is not None:
                sm.ifcs[i].set_max_aperture(max_ap)


def set_ape(opt_model, avoid_list=None, include_list=None):
    """ From existing fields and vignetting, calculate clear apertures. 
    
    This function modifies the max_aperture maintained by the list of
    :class:`~.interface.Interface` in the 
    :class:`~.sequential.SequentialModel`. For each interface, the smallest 
    aperture that will pass all of the (vignetted) boundary rays, for each 
    field, is chosen.
    
    The change of the apertures is propagated to the 
    :class:`~.elements.ElementModel` via 
    :meth:`~.elements.ElementModel.sync_to_seq`.
    """
    set_clear_apertures(opt_model, avoid_list, include_list)

    # sync the element model with the new clear apertures
    opt_model['em'].sync_to_seq(opt_model['sm'])


def set_vig(opm, **kwargs):
    """ From existing fields and clear apertures, calculate vignetting. """
    osp = opm['osp']
    for fi in range(len(osp['fov'].fields)):
        fld, wvl, foc = osp.lookup_fld_wvl_focus(fi)
        logger.debug(f"set vig field {fi}:")
        calc_vignetting_for_field(opm, fld, wvl, **kwargs)


def set_stop_aperture(opm, **kwargs):
    """ Set the aperture on the stop surface to satisfy the pupil spec.

    The vignetting is recalculated after the stop aperture change.
    """
    sm = opm['seq_model']
    # clear the axial vignetting so the pupil_spec defines the axial marginal rays
    opm['osp']['fov']['axis'].clear_vignetting()
    # now set the aperture at the stop surface to the pupil spec defined size
    set_clear_apertures(opm, include_list=[sm.stop_surface])
    # set vignetting to account for the stop aperture change
    set_vig(opm)


def set_pupil(opm, use_parax=False):
    """ From existing stop size, calculate pupil spec and vignetting. 
    
    Use the upper Y marginal ray on-axis (field #0) and iterate until it
    goes through the edge of the stop surface. Use the object or image 
    segments of this ray to update the pupil specification value 
    e.g. EPD, NA or f/#.
    """
    sm = opm['sm']
    if sm.stop_surface is None:
        # Nope, the whole purpose here is to go from aperture stop to pupil
        print('floating stop surface')
        return
    idx_stop = sm.stop_surface
    osp = opm['osp']

    # iterate the on-axis marginal ray thru the edge of the stop.
    fld_0, cwl, foc = osp.lookup_fld_wvl_focus(0)
    stop_radius = sm.ifcs[idx_stop].surface_od()
    start_coords = iterate_pupil_ray(opm, sm.stop_surface, 1, 1.0, 
                                     stop_radius, fld_0, cwl)
    
    logger.debug(f"set_pupil edge of stop coords: {start_coords[0]:8.4f} "
                 f"{start_coords[1]:8.4f}")
    
    # trace the real axial marginal ray
    ray_result = trace.trace_safe(opm, start_coords, fld_0, cwl, 
                                  None, 'full', apply_vignetting=False)
    ray_pkg, ray_err = ray_result

    obj_img_key = osp['pupil'].key[0]
    pupil_spec = osp['pupil'].key[1]
    pupil_value_orig = osp['pupil'].value

    ax_ray, pr_ray, fod = opm['ar']['parax_data']
    if use_parax:
        scale_ratio = stop_radius/ax_ray[idx_stop][mc.ht]
        logger.debug(f"{scale_ratio=:8.5f} (parax)")
        if obj_img_key == 'object':
            if pupil_spec == 'epd':
                osp['pupil'].value = scale_ratio*(2*fod.enp_radius)
            else:
                slp0 = scale_ratio*ax_ray[0][mc.slp]
                if pupil_spec == 'NA':
                    n0 = sm.central_rndx[0]
                    rs0 = RaySeg(*ray_pkg[0][0])
                    osp['pupil'].value = n0*rs0.d[1]
                    # osp['pupil'].value = etendue.slp2na(slp0)
                elif pupil_spec == 'f/#':
                    osp['pupil'].value = 1/(2*slp0)
        elif obj_img_key == 'image':
            if pupil_spec == 'epd':
                osp['pupil'].value = scale_ratio*(2*fod.exp_radius)
            else:
                slpk = scale_ratio*ax_ray[-1][mc.slp]
                if pupil_spec == 'NA':
                    nk = sm.central_rndx[-1]
                    rsm2 = RaySeg(*ray_pkg[0][-2])
                    osp['pupil'].value = -nk*rsm2.d[1]
                    # osp['pupil'].value = etendue.slp2na(slpk)
                elif pupil_spec == 'f/#':
                    osp['pupil'].value = -1/(2*slpk)
    else:  # use real marginal ray
        scale_ratio = ray_pkg[0][1][0][1]/ax_ray[1][mc.ht]
        logger.debug(f"{scale_ratio=:8.5f}")
        if obj_img_key == 'object':
            if pupil_spec == 'epd':
                rs1 = RaySeg(*ray_pkg[0][1])
                ht = rs1.p[1]
                osp['pupil'].value *= scale_ratio
            else:
                rs0 = RaySeg(*ray_pkg[0][0])
                slp0 = rs0.d[1]/rs0.d[2]
                if pupil_spec == 'NA':
                    n0 = sm.central_rndx[0]
                    osp['pupil'].value = n0*rs0.d[1]
                elif pupil_spec == 'f/#':
                    osp['pupil'].value = 1/(2*slp0)
        elif obj_img_key == 'image':
            rsm2 = RaySeg(*ray_pkg[0][-2])
            if pupil_spec == 'epd':
                ht = rsm2.p[1]
                osp['pupil'].value = 2*ht
            else:
                slpk = scale_ratio*ax_ray[-1][mc.slp]
                if pupil_spec == 'NA':
                    nk = sm.central_rndx[-1]
                    osp['pupil'].value = -nk*rsm2.d[1]
                    # osp['pupil'].value = etendue.slp2na(slpk)
                elif pupil_spec == 'f/#':
                    osp['pupil'].value = -1/(2*slpk)
    
    # trace the real axial marginal ray with aperture clipping
    clipped_rr = trace.trace_safe(opm, start_coords, fld_0, cwl, 
                                  None, 'full', apply_vignetting=False, 
                                  check_apertures=True)
    clipped_ray_pkg, clipped_ray_err = clipped_rr
    if clipped_ray_err is not None:
        if isinstance(clipped_ray_err, terr.TraceRayBlockedError):
            print(f"Axial bundle limited by surface {clipped_ray_err.surf}, "
                  "not stop surface.")
            logger.warning("Axial bundle limited by surface "
                           f"{clipped_ray_err.surf}, not stop surface.")
    if pupil_value_orig != osp['pupil'].value:
        opm.update_model()
        set_vig(opm)


def calc_vignetting_for_field(opm, fld, wvl, **kwargs):
    """Calculate and set the vignetting parameters for `fld`. """
    vg_kwargs = {}
    if 'max_iter_count' in kwargs:
        vg_kwargs['max_iter_count'] = kwargs.get('max_iter_count')
    pupil_starts = opm['osp']['pupil'].pupil_rays[1:]
    vig_factors = [0.]*4
    for i in range(4):
        xy = i//2
        start = pupil_starts[i]
        vig, clip_indx, ray_pkg = calc_vignetted_ray(
            opm, xy, start, fld, wvl, **vg_kwargs)
        vig_factors[i] = vig

    # update the field's vignetting factors
    fld.vux = vig_factors[0]
    fld.vlx = vig_factors[1]
    fld.vuy = vig_factors[2]
    fld.vly = vig_factors[3]


def calc_vignetted_ray(opm, xy, start_dir, fld, wvl, max_iter_count=50):
    """ Find the limiting aperture and return the vignetting factor. 

    Args:
        opm: :class:`~.OpticalModel` instance
        xy: 0 or 1 depending on x or y axis as the pupil direction
        start_dir: the unit length starting pupil coordinates, e.g [1., 0.]. 
                   This establishes the radial direction of the ray iteration.
        fld: :class:`~.Field` point for wave aberration calculation
        wvl: wavelength of ray (nm)
        max_iter_count: fail-safe limit on aperture search

    Returns:
        (**vig**, **clip_indx**, **ray_pkg**)

        - **vig** - vignetting factor
        - **clip_indx** - the index of the limiting interface
        - **ray_pkg** - the vignetting-limited ray
 
    """
    logger.debug(f"fld={fld.yf:5.2f}, [{start_dir[0]:5.2f}, "
                 f"{start_dir[1]:5.2f}]")
    rel_p1 = np.array(start_dir)
    sm = opm['sm']
    still_iterating = True
    clip_indx = None
    iter_count = 0  # safe guard against runaway iteration
    while still_iterating and iter_count<max_iter_count:
        iter_count += 1
        try:
            ray_pkg = trace.trace_base(opm, rel_p1, fld, wvl, 
                                       apply_vignetting=False, 
                                       check_apertures=True,
                                       pt_inside_fuzz=1e-4)

        except terr.TraceError as ray_error:
            ray_pkg = RayPkg(*ray_error.ray_pkg)
            indx = ray_error.surf
            if indx == clip_indx:
                r_target = sm.ifcs[clip_indx].edge_pt_target(start_dir)
                p = ray_pkg[mc.ray][clip_indx][mc.p]
                r_ray = copysign(sqrt(p[0]**2 + p[1]**2), r_target[xy])
                r_error = r_ray - r_target[xy]
                logger.debug(f" A {xy_str[xy]} = {rel_p1[xy]:10.6f}:   "
                             f"blocked at {clip_indx}, del={r_error:8.1e}, "
                             "exiting")
                still_iterating = False
            else:
                r_target = sm.ifcs[indx].edge_pt_target(start_dir)
                logger.debug(f" B {xy_str[xy]} = {rel_p1[xy]:10.6f}:   "
                             f"blocked at {indx}. target={r_target[xy]:9.6f}")
                rel_p1 = iterate_pupil_ray(opm, indx, xy, rel_p1[xy], 
                                           r_target[xy], fld, wvl)
                still_iterating = True
                clip_indx = indx
        else:  # ray successfully traced.
            if clip_indx is not None:
                # fall through and exit
                r_target = sm.ifcs[clip_indx].edge_pt_target(start_dir)
                p = ray_pkg[mc.ray][clip_indx][mc.p]
                r_ray = copysign(sqrt(p[0]**2 + p[1]**2), r_target[xy])
                r_error = r_ray - r_target[xy]
                logger.debug(f" C {xy_str[xy]} = {rel_p1[xy]:10.6f}:   "
                             f"blocked at {clip_indx}, del={r_error:8.1e}, "
                             "exiting")
                still_iterating = False
            else: # this is the first time through
                # iterate to find the ray that goes through the edge
                # of the stop surface
                indx = stop_indx = sm.stop_surface
                if stop_indx is not None:
                    r_target = sm.ifcs[stop_indx].edge_pt_target(start_dir)
                    logger.debug(f" D {xy_str[xy]} = {rel_p1[xy]:10.6f}:   "
                                 f"passed first time, iterate to edge of stop, "
                                 f"ifcs[{stop_indx}]")
                    rel_p1 = iterate_pupil_ray(opm, indx, xy, rel_p1[xy], 
                                               r_target[xy], fld, wvl)
                    still_iterating = True
                    clip_indx = indx
                else: # floating stop, exit
                    still_iterating = False

    vig = 1.0 - (rel_p1[xy]/start_dir[xy])
    logger.info(f" ray: ({start_dir[0]:2.0f}, {start_dir[1]:2.0f}), "
                f"vig={vig:8.4f}, limited at ifcs[{clip_indx}]")
    return vig, clip_indx, ray_pkg


def calc_vignetted_ray_by_bisection(opm, xy, start_dir, fld, wvl, 
                                    max_iter_count=10):
    """ Find the limiting aperture and return the vignetting factor. 

    Args:
        opm: :class:`~.OpticalModel` instance
        xy: 0 or 1 depending on x or y axis as the pupil direction
        start_dir: the unit length starting pupil coordinates, e.g [1., 0.]. 
                   This establishes the radial direction of the ray iteration.
        fld: :class:`~.Field` point for wave aberration calculation
        wvl: wavelength of ray (nm)
        max_iter_count: fail-safe limit on aperture search

    Returns:
        (**vig**, **clip_indx**, **ray_pkg**)

        - **vig** - vignetting factor
        - **clip_indx** - the index of the limiting interface
        - **ray_pkg** - the vignetting-limited ray
 
    """
    logger.debug(f"fld={fld.yf:5.2f}, [{start_dir[0]:5.2f}, "
                 f"{start_dir[1]:5.2f}]")
    rel_p1 = np.array(start_dir)
    clip_indx = None
    iter_count = 0  # safe guard against runaway iteration
    step_size = 1.0
    while iter_count<max_iter_count:
        iter_count += 1
        try:
            step_size /= 2
            ray_pkg = trace.trace_base(opm, rel_p1, fld, wvl, 
                                       apply_vignetting=False, 
                                       check_apertures=True,
                                       pt_inside_fuzz=1e-4)

        except terr.TraceError as ray_error:
            ray_pkg = RayPkg(*ray_error.ray_pkg)
            clip_indx = ray_error.surf
            rel_p1 = -step_size*np.array(start_dir) + rel_p1
            logger.debug(f"{xy_str[xy]} = {rel_p1[xy]:10.6f}: "
                         f"blocked at {clip_indx}")
        else:  # ray successfully traced.
            rel_p1 = step_size*np.array(start_dir) + rel_p1
            logger.debug(f"{xy_str[xy]} = {rel_p1[xy]:10.6f}: passed")

    vig = 1.0 - (rel_p1[xy]/start_dir[xy])
    logger.debug(f"   {vig=:7.4f}, {clip_indx=}")
    return vig, clip_indx, ray_pkg


def iterate_pupil_ray(opt_model, indx, xy, start_r0, r_target, 
                      fld, wvl, **kwargs):
    """ iterates a ray to r_target on interface indx, returns aim points on
    the paraxial entrance pupil plane

    If indx is None, i.e. a floating stop surface, returns r_target.

    If the iteration fails, a :class:`~.traceerror.TraceError` will be raised

    Args:
        opm: :class:`~.OpticalModel` instance
        indx: index of interface whose edge is the iteration target
        xy: 0 or 1 depending on x or y axis as the pupil direction
        start_r0: iteration starting point
        r_target: clear aperture radius that is the iteration target.
        fld: :class:`~.Field` point for wave aberration calculation
        wvl: wavelength of ray (nm)

    Returns:
        start_coords: pupil coordinates for ray thru r_target on ifc indx.

    """

    def r_pupil_coordinate(xy_coord, *args):
        opt_model, indx, xy, fld, wvl, r_target = args

        rel_p1 = np.array([0., 0.])
        rel_p1[xy] = xy_coord
        try:
            ray_pkg = trace.trace_base(opt_model, rel_p1, fld, wvl, 
                                       apply_vignetting=False, 
                                       check_apertures=False)
        except terr.TraceError as ray_error:
            ray_pkg = ray_error.ray_pkg
            if isinstance(ray_error, terr.TraceMissedSurfaceError):
                # no surface intersection, so no ray data at indx
                if ray_error.surf <= indx:
                    ray_error.rel_p1 = rel_p1
                    raise ray_error
            else:
                if ray_error.surf < indx:
                    ray_error.rel_p1 = rel_p1
                    raise ray_error

        # compute the radial distance to the intersection point
        p = ray_pkg[mc.ray][indx][mc.p]
        r_ray = copysign(sqrt(p[0]**2 + p[1]**2), r_target)
        delta = r_ray - r_target
        logger.debug(f"  {xy_coord=:8.5f}   {r_ray=:8.5f}    "
                     f"delta={delta:9.2g}")
        return delta

    start_coords = np.array([0., 0.])
    if indx is not None:
        logging.captureWarnings(True)
        try:
            start_r, results = newton(r_pupil_coordinate, start_r0,
                                      args=(opt_model, indx, xy,
                                            fld, wvl, r_target), tol=1e-6,
                                      disp=False, full_output=True)
        except RuntimeError as rte:
            # if we come here, set start_r to a RuntimeResults object
            start_r = results.root
            print(f"vigcalc.iterate_pupil_ray {rte=}")
        except terr.TraceError as rt_err:
            logger.debug(f"  {type(rt_err).__name__}: surf={rt_err.surf}    "
                         f"rel_p1={rt_err.rel_p1[xy]=:8.5f}   ")
            start_r = 0.9*rt_err.rel_p1[xy]
        start_coords[xy] = start_r

    else:  # floating stop surface - use entrance pupil for aiming
        start_coords[xy] = r_target

    return start_coords
