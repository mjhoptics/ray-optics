#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2022 Michael J. Hayford
""" Vignetting and clear aperture setting operations

.. Created on Mon Apr 18 15:28:25 2022

.. codeauthor: Michael J. Hayford
"""
import logging

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


def set_ape(opm):
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
    rayset = trace.trace_boundary_rays(opm, use_named_tuples=True)

    for i, ifc in enumerate(opm['sm'].ifcs):
        max_ap = -1.0e+10
        update = True
        for f in rayset:
            for p in f:
                ray = p.ray
                if len(ray) > i:
                    ap = sqrt(ray[i].p[0]**2 + ray[i].p[1]**2)
                    if ap > max_ap:
                        max_ap = ap
                else:  # ray failed before this interface, don't update
                    update = False
        if update:
            ifc.set_max_aperture(max_ap)

    # sync the element model with the new clear apertures
    opm['em'].sync_to_seq(opm['sm'])


def set_vig(opm, **kwargs):
    """ From existing fields and clear apertures, calculate vignetting. """
    osp = opm['osp']
    for fi in range(len(osp['fov'].fields)):
        fld, wvl, foc = osp.lookup_fld_wvl_focus(fi)
        logger.debug(f"set vig field {fi}:")
        calc_vignetting_for_field(opm, fld, wvl, **kwargs)


def set_pupil(opm):
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
    osp = opm['osp']

    # iterate the on-axis marginal ray thru the edge of the stop.
    fld, wvl, foc = osp.lookup_fld_wvl_focus(0)
    stop_radius = sm.ifcs[sm.stop_surface].surface_od()
    start_coords = iterate_pupil_ray(opm, sm.stop_surface, 1, 1.0, 
                                     stop_radius, fld, wvl)

    logger.debug(f"set_pupil edge of stop coords: {start_coords[0]:8.4f} "
                 f"{start_coords[1]:8.4f}")
    
    # trace the real axial marginal ray
    ray_result = trace.trace_safe(opm, start_coords, fld, wvl, 
                                  None, None, apply_vignetting=False, 
                                  check_apertures=True)
    ray_pkg, ray_err = ray_result

    obj_img_key = osp['pupil'].key[0]
    pupil_spec = osp['pupil'].key[1]
    pupil_value_orig = osp['pupil'].value

    if obj_img_key == 'object':
        if pupil_spec == 'epd':
            rs1 = RaySeg(*ray_pkg[0][1])
            ht = rs1.p[1]
            osp['pupil'].value = 2*ht
        else:
            rs0 = RaySeg(*ray_pkg[0][0])
            slp0 = rs0.d[1]/rs0.d[2]
            if pupil_spec == 'NA':
                n0 = sm.rindx[0]
                osp['pupil'].value = n0*rs0.d[1]
                # osp['pupil'].value = etendue.slp2na(slp0)
            elif pupil_spec == 'f/#':
                osp['pupil'].value = 1/(2*slp0)
    elif obj_img_key == 'image':
        rsm2 = RaySeg(*ray_pkg[0][-2])
        if pupil_spec == 'epd':
            ht = rsm2.p[1]
            osp['pupil'].value = 2*ht
        else:
            slpk = rsm2.d[1]/rsm2.d[2]
            if pupil_spec == 'NA':
                nk = sm.rindx[-1]
                osp['pupil'].value = -nk*rsm2.d[1]
                # osp['pupil'].value = etendue.slp2na(slpk)
            elif pupil_spec == 'f/#':
                osp['pupil'].value = -1/(2*slpk)

    if pupil_value_orig != osp['pupil'].value:
        opm.update_model()
        set_vig(opm)


def calc_vignetting_for_field(opm, fld, wvl, **kwargs):
    """Calculate and set the vignetting parameters for `fld`. """
    use_bisection = kwargs.get('use_bisection', 
                               opm['osp']['fov'].is_wide_angle)
    pupil_starts = opm['osp']['pupil'].pupil_rays[1:]
    vig_factors = [0.]*4
    for i in range(4):
        xy = i//2
        start = pupil_starts[i]
        if use_bisection:
            vig, last_indx, ray_pkg = calc_vignetted_ray_by_bisection(
                opm, xy, start, fld, wvl)
        else:
            vig, last_indx, ray_pkg = calc_vignetted_ray(
                opm, xy, start, fld, wvl)
        logger.debug(f"ray: ({start[0]:2.0f}, {start[1]:2.0f}), "
                     f"vig={vig:8.4f}, limited at ifcs[{last_indx}]")
        vig_factors[i] = vig

    # update the field's vignetting factors
    fld.vux = vig_factors[0]
    fld.vlx = vig_factors[1]
    fld.vuy = vig_factors[2]
    fld.vly = vig_factors[3]


def calc_vignetted_ray(opm, xy, start_dir, fld, wvl, max_iter_count=10):
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
        (**vig**, **last_indx**, **ray_pkg**)

        - **vig** - vignetting factor
        - **last_indx** - the index of the limiting interface
        - **ray_pkg** - the vignetting-limited ray
 
    """
    rel_p1 = np.array(start_dir)
    sm = opm['sm']
    still_iterating = True
    last_indx = None
    iter_count = 0  # safe guard against runaway iteration
    while still_iterating and iter_count<max_iter_count:
        iter_count += 1
        try:
            ray_pkg = trace.trace_base(opm, rel_p1, fld, wvl, 
                                       apply_vignetting=False, 
                                       check_apertures=True,
                                       pt_inside_fuzz=1e-4)

        except terr.TraceError as ray_error:
            indx = ray_error.surf
            ray_pkg = ray_error.ray_pkg
            logger.debug(f"{xy_str[xy]} = {rel_p1[xy]:10.6f}: "
                         f"blocked at {indx}")
            if indx == last_indx:
                still_iterating = False
            else:
                r_target = sm.ifcs[indx].edge_pt_target(start_dir)
                rel_p1 = iterate_pupil_ray(opm, indx, xy, rel_p1[xy], 
                                           r_target[xy], fld, wvl)
                still_iterating = True
                last_indx = indx
        else:  # ray successfully traced.
            logger.debug(f"{xy_str[xy]} = {rel_p1[xy]:10.6f}: passed")
            if last_indx is not None:
                # fall through and exit
                still_iterating = False
            else: # this is the first time through
                # iterate to find the ray that goes through the edge
                # of the stop surface
                indx = stop_indx = sm.stop_surface
                if stop_indx is not None:
                    r_target = sm.ifcs[stop_indx].edge_pt_target(start_dir)
                    rel_p1 = iterate_pupil_ray(opm, indx, xy, rel_p1[xy], 
                                               r_target[xy], fld, wvl)
                    still_iterating = True
                    last_indx = indx
                else: # floating stop, exit
                    still_iterating = False

    vig = 1.0 - (rel_p1[xy]/start_dir[xy])
    return vig, last_indx, ray_pkg


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
        (**vig**, **last_indx**, **ray_pkg**)

        - **vig** - vignetting factor
        - **last_indx** - the index of the limiting interface
        - **ray_pkg** - the vignetting-limited ray
 
    """
    rel_p1 = np.array(start_dir)
    sm = opm['sm']
    still_iterating = True
    last_indx = None
    iter_count = 0  # safe guard against runaway iteration
    step_size = 1.0
    while still_iterating and iter_count<max_iter_count:
        iter_count += 1
        try:
            step_size /= 2
            ray_pkg = trace.trace_base(opm, rel_p1, fld, wvl, 
                                       apply_vignetting=False, 
                                       check_apertures=True,
                                       pt_inside_fuzz=1e-4)

        except terr.TraceError as ray_error:
            ray_pkg = RayPkg(*ray_error.ray_pkg)
            last_indx = ray_error.surf
            rel_p1 = -step_size*np.array(start_dir) + rel_p1
            logger.debug(f"{xy_str[xy]} = {rel_p1[xy]:10.6f}: "
                         f"blocked at {last_indx}")
        else:  # ray successfully traced.
            rel_p1 = step_size*np.array(start_dir) + rel_p1
            logger.debug(f"{xy_str[xy]} = {rel_p1[xy]:10.6f}: passed")

    vig = 1.0 - (rel_p1[xy]/start_dir[xy])
    return vig, last_indx, ray_pkg


def iterate_pupil_ray(opt_model, indx, xy, start_r0, r_target, fld, wvl, **kwargs):
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
                    raise ray_error
            else:
                if ray_error.surf < indx:
                    raise ray_error

        ray = ray_pkg[mc.ray]
        r_ray = ray[indx][mc.p][xy]
        return r_ray - r_target

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
        except terr.TraceError:
            start_r = 0.0
        start_coords[xy] = start_r
        # print(f"converged={results.converged} in {results.iterations} iterations")

    else:  # floating stop surface - use entrance pupil for aiming
        start_coords[xy] = r_target

    return start_coords
