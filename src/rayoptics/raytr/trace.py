#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Supports model ray tracing in terms of relative aperture and field.

.. Created on Mon Sep 17 23:10:59 2018

.. codeauthor: Michael J. Hayford
"""

import itertools
import logging
import warnings
import math
import numpy as np
from numpy.linalg import norm
from scipy.optimize import newton, fsolve
import pandas as pd

from . import raytrace as rt
from . import RayPkg, RaySeg, RayResult
from .waveabr import (wave_abr_full_calc, calculate_reference_sphere, 
                      transfer_to_exit_pupil)
from .wideangle import find_real_enp, enp_z_coordinate
from rayoptics.optical import model_constants as mc
from .traceerror import TraceError, TraceMissedSurfaceError, TraceTIRError
from rayoptics.util.misc_math import normalize

logger = logging.getLogger(__name__)

def ray_pkg(ray_pkg):
    """ return a |Series| containing a ray package (RayPkg) """
    return pd.Series(ray_pkg, index=['ray', 'op', 'wvl'])


def ray_df(ray):
    """ return a |DataFrame| containing ray data """
    r = pd.DataFrame(ray, columns=['inc_pt', 'after_dir',
                                   'after_dst', 'normal'])
    r.index.names = ['intrfc']
    return r


def list_ray(ray_obj, tfrms=None, start=0):
    """ pretty print a ray either in local or global coordinates 
    
    The input ray_obj can be either the return from trace_ray(), i.e.
    a (ray_pkg, ray_err) tuple or a ray_pkg, i.e. (`ray`, opl, wvl) 
    or a `ray` alone.
    """
    ray_err = None
    if isinstance(ray_obj, tuple):
        if len(ray_obj) == 2:
            ray_obj, ray_err = ray_obj
        ray = ray_obj[0]
    else:
        ray = ray_obj

    colHeader = "            X            Y            Z           L" \
                "            M            N               Len"
    print(colHeader)

    colFormats = "{:3d}: {:12.5f} {:12.5f} {:12.5g} {:12.6f} {:12.6f} " \
                 "{:12.6f} {:12.5g}"

    for i, r in enumerate(ray[start:], start=start):
        if tfrms is None:
            print(colFormats.format(i,
                                    r[mc.p][0], r[mc.p][1], r[mc.p][2],
                                    r[mc.d][0], r[mc.d][1], r[mc.d][2],
                                    r[mc.dst]))
        else:
            rot, trns = tfrms[i]
            p = rot.dot(r[mc.p]) + trns
            d = rot.dot(r[mc.d])
            print(colFormats.format(i, p[0], p[1], p[2], d[0], d[1], d[2],
                                    r[mc.dst]))
    if ray_err is not None:
        print(f"ray failure: {type(ray_err).__name__}")


def list_in_out_dir(path, ray):
    """ list the incident and exiting ray direction cosines.  """
    lcl_tfrms = [path_seg[mc.Tfrm] for path_seg in path]
    print("                  in_dir              |              out_dir")
    before_dir = ray[0][1]
    print(f"{0:2d}:                                   |"
          f"{before_dir[0]:10.6f} {before_dir[1]:10.6f} {before_dir[2]:10.6f}")
    for i, rst in enumerate(zip(ray[1:], lcl_tfrms), start=1):
        seg, tfrm = rst
        after_dir = seg[1]
        rt, t = tfrm
        b4_dir = rt.dot(before_dir)
        print(f"{i:2d}: {b4_dir[0]:10.6f} {b4_dir[1]:10.6f} {b4_dir[2]:10.6f}  |"
              f"{after_dir[0]:10.6f} {after_dir[1]:10.6f} {after_dir[2]:10.6f}")
        before_dir = after_dir

    rt, t = lcl_tfrms[-1]
    b4_dir = rt.dot(before_dir)
    print(f"{i+1:2d}: {b4_dir[0]:10.6f} {b4_dir[1]:10.6f} {b4_dir[2]:10.6f}  |")


def trace_ray(opt_model, pupil, fld, wvl, 
              output_filter=None, rayerr_filter='full',
               **kwargs) -> RayResult:
    """ Trace a single ray via pupil, field and wavelength specs.
    
    This function traces a single ray at a given wavelength, pupil and field specification. 

    Ray failures (miss surface, TIR) and aperture clipping are handled via RayError exceptions. If a failure occurs, a second item is returned (if  *rayerr_filter* is set to 'summary' or 'full') that contains information about the failure. Apertures are tested using the :meth:`~.seq.interface.Interface.point_inside` API when *check_apertures* is True. 

    The pupil coordinates by default are normalized to the vignetted pupil extent. Alternatively, the pupil coordinates can be taken as actual coordinates on the pupil plane (and similarly for ray direction) using the **pupil_type** keyword argument.

    The amount of output that is returned can range from the entire ray (default) to the image segment only or even the return from a user-supplied filtering function.

    Args:
        opt_model: :class:`~.OpticalModel` instance
        pupil: 2d vector of relative pupil coordinates
        fld: :class:`~.Field` point for wave aberration calculation
        wvl: wavelength of ray (nm)

        check_apertures: if True, do point_inside() test on inc_pt
        apply_vignetting: if True, apply the `fld` vignetting factors to **pupil**

        pupil_type: ::

            - 'rel pupil': relative pupil coordinates
            - 'aim pt': aim point on pupil plane
            - 'aim dir': aim direction in object space

        use_named_tuples: if True, returns data as RayPkg and RaySeg.

        output_filter: ::

            - if None, append entire ray
            - if 'last', append the last ray segment only
            - else treat as callable and append the return value

        rayerr_filter: ::

            - if None, on ray error append nothing
            - if 'summary', append the exception without ray data
            - if 'full', append the exception with ray data up to error
            - else append nothing

        eps: accuracy tolerance for surface intersection calculation

    Returns:
        tuple: ray_pkg, trace_error | None

    """
    unt = kwargs.pop('use_named_tuples', True)
    
    ray_result = trace_safe(opt_model, pupil, fld, wvl, 
                            output_filter, rayerr_filter, 
                            use_named_tuples=unt, **kwargs)
    return ray_result


def trace_safe(opt_model, pupil, fld, wvl,
               output_filter, rayerr_filter, **kwargs) -> RayResult:
    """Wrapper for trace_base that handles exceptions.
    
    Args:
        opt_model: :class:`~.OpticalModel` instance
        pupil: 2d vector of relative pupil coordinates
        fld: :class:`~.Field` point for wave aberration calculation
        wvl: wavelength of ray (nm)
        output_filter: ::

            - if None, append entire ray
            - if 'last', append the last ray segment only
            - else treat as callable and append the return value

        rayerr_filter: ::

            - if None, on ray error append nothing
            - if 'summary', append the exception without ray data
            - if 'full', append the exception with ray data up to error
            - else append nothing

    Returns:
        ray_result: see discussion of filters, above.

    """
    use_named_tuples = kwargs.get('use_named_tuples', False)

    ray_result = None, None

    try:
        ray_pkg = trace_base(opt_model, pupil, fld, wvl, **kwargs)

    except TraceError as rayerr:
        if rayerr_filter is None:
            pass
        elif rayerr_filter == 'full':
            ray, op_delta, wvl = rayerr.ray_pkg
            ray = [RaySeg(*rs) for rs in ray]
            rayerr.ray_pkg = RayPkg(ray, op_delta, wvl)
            ray_result = rayerr.ray_pkg, rayerr
        elif rayerr_filter == 'summary':
            rayerr.ray_pkg = None
            ray_result = rayerr.ray_pkg, rayerr
        else:
            pass
    else:
        if use_named_tuples:
            ray, op_delta, wvl = ray_pkg
            ray = [RaySeg(*rs) for rs in ray]
            ray_pkg = RayPkg(ray, op_delta, wvl)

        if output_filter is None:
            ray_result = ray_pkg, None
        elif output_filter == 'last':
            ray, op_delta, wvl = ray_pkg
            final_seg_pkg = RayPkg([ray[-1]], op_delta, wvl)
            ray_result = final_seg_pkg, None
        else:
            ray_result = output_filter(ray_pkg), None

    return RayResult(*ray_result)


def trace(seq_model, pt0, dir0, wvl, **kwargs) -> RayPkg:
    """ returns (ray, ray_opl, wvl)

    Args:
        seq_model: the :class:`~.SequentialModel` to be traced
        pt0: starting coordinate at object interface
        dir0: starting direction cosines following object interface
        wvl: ray trace wavelength in nm
        **kwargs: keyword arguments

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
    return RayPkg(*rt.trace(seq_model, pt0, dir0, wvl, **kwargs))


def trace_base(opt_model, pupil, fld, wvl, 
               apply_vignetting=True, pupil_type='rel pupil', 
               **kwargs) -> RayPkg:
    """Trace ray specified by relative aperture and field point.
    
    `pupil_type` controls how `pupil` data is interpreted when calculating the starting ray coordinates.

    Args:
        opt_model: instance of :class:`~.OpticalModel` to trace
        pupil: aperture coordinates of ray
        fld: instance of :class:`~.Field`
        wvl: ray trace wavelength in nm
        apply_vignetting: if True, apply the `fld` vignetting factors to **pupil**
        pupil_type: ::

            - 'rel pupil': relative pupil coordinates
            - 'aim pt': aim point on pupil plane
            - 'aim dir': aim direction in object space
    
        **kwargs: keyword arguments

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
    if pupil_type == 'rel pupil':
        pupil_coords = (fld.apply_vignetting(pupil) if apply_vignetting 
                        else pupil)
    else:
        pupil_coords = pupil

    pt0, dir0 = opt_model['optical_spec'].ray_start_from_osp(pupil_coords, fld,
                                                             pupil_type)

    # if wide_angle, don't try to intercept object and don't disallow 
    # propagation against z_dir; this will be the case for rays exceeding 
    # 90 degrees at the first surface.
    if opt_model['optical_spec']['fov'].is_wide_angle:
        kwargs['intersect_obj'] = False
    else:  # otherwise, if not wide angle, propagation against z_dir means 
        # a virtual object. To handle virtual object distances, always 
        # propagate from the object in a positive Z direction.
        if dir0[2] * opt_model['seq_model'].z_dir[0] < 0:
            dir0 = -dir0

    return rt.trace(opt_model['seq_model'], pt0, dir0, wvl, **kwargs)


def iterate_ray(opt_model, ifcx, xy_target, fld, wvl, **kwargs):
    """ iterates a ray to xy_target on interface ifcx, returns aim points on
    the paraxial entrance pupil plane

    If idcx is None, i.e. a floating stop surface, returns xy_target.

    If the iteration fails, a TraceError will be raised
    """
    rr = None
    def y_stop_coordinate(y1, *args):
        nonlocal rr
        seq_model, ifcx, pt0, obj2enp_dist, wvl, y_target, not_wa = args
        pt1 = np.array([0., y1, obj2enp_dist])
        dir0 = normalize(pt1 - pt0)
        # handle case where entrance pupil is behind the object (issue #120)
        if not_wa and dir0[2] * seq_model.z_dir[0] < 0:
            dir0 = -dir0

        try:
            ray_pkg = RayPkg(*rt.trace(seq_model, pt0, dir0, wvl))

        except TraceError as ray_error:
            logger.debug(f'ray_error: "{type(ray_error).__name__}", '
                        f'{ray_error.surf=}')
            ray_pkg = RayPkg(*ray_error.ray_pkg)
            rr = RayResult(ray_pkg, ray_error)
            final_coord = np.array([0., 0., 0.])
            if ray_error.surf < ifcx:
                raise ray_error
 
        else:
            rr = RayResult(ray_pkg, None)
            final_coord = ray_pkg.ray[ifcx][mc.p]

        y_ray = final_coord[1]
        return y_ray - y_target

    def surface_coordinate(coord, *args):
        nonlocal rr
        seq_model, ifcx, pt0, dist, wvl, target, not_wa = args
        pt1 = np.array([coord[0], coord[1], dist])
        dir0 = normalize(pt1 - pt0)
        if not_wa and dir0[2] * seq_model.z_dir[0] < 0:
            dir0 = -dir0
        
        try:
            ray_pkg = RayPkg(*rt.trace(seq_model, pt0, dir0, wvl))
            
        except TraceError as ray_error:
            ray_pkg = RayPkg(*ray_error.ray_pkg)
            rr = RayResult(ray_pkg, ray_error)
            final_coord = np.array([0., 0., 0.])
            if ray_error.surf < ifcx:
                raise ray_error
        else:
            rr = RayResult(ray_pkg, None)
            final_coord = ray_pkg.ray[ifcx][mc.p]
        xy_ray = np.array([final_coord[0], final_coord[1]])
        return xy_ray - target

    seq_model = opt_model['seq_model']
    osp = opt_model['optical_spec']

    fod = opt_model['analysis_results']['parax_data'].fod
    obj2enp_dist = fod.obj_dist + fod.enp_dist
    
    not_wa = not osp['fov'].is_wide_angle

    pt0, d0 = osp.obj_coords(fld)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if ifcx is not None:
            if pt0[0] == 0.0 and xy_target[0] == 0.0:
                # do 1D iteration if field and target points are zero in x
                y_target = xy_target[1]
                try:
                    start_y, results = newton(y_stop_coordinate, 0.,
                                              args=(seq_model, ifcx, pt0,
                                                    obj2enp_dist, wvl, y_target, not_wa),
                                              disp=False, full_output=True)
                except RuntimeError as rte:
                    # if we come here, start_y is a RuntimeResults object
                    # print(rte)
                    start_y = results.root
                except TraceError:
                    start_y = 0.0
                start_coords = np.array([0., start_y])
            else:
                # do 2D iteration. epsfcn is a parameter increment,
                #  make proportional to pupil radius
                try:
                    start_coords = fsolve(surface_coordinate, 
                                          np.array([0., 0.]),
                                          epsfcn=0.0001*fod.enp_radius,
                                          args=(seq_model, ifcx, pt0, 
                                                obj2enp_dist, wvl, 
                                                xy_target, not_wa))
                except TraceError:
                    start_coords = np.array([0., 0.])
        else:  # floating stop surface - use entrance pupil for aiming
            start_coords = np.array([0., 0.]) + xy_target

    return start_coords, rr


def trace_with_opd(opt_model, pupil, fld, wvl, foc, **kwargs):
    """ returns (ray, ray_opl, wvl, opd) """

    chief_ray_pkg = get_chief_ray_pkg(opt_model, fld, wvl, foc)
    image_pt_2d = kwargs.get('image_pt', None)
    image_delta = kwargs.get('image_delta', None)
    ref_sphere = calculate_reference_sphere(opt_model, fld, wvl, foc,
                                            chief_ray_pkg,
                                            image_pt_2d=image_pt_2d,
                                            image_delta=image_delta)

    ray_pkg, ray_err = trace_ray(opt_model, pupil, fld, wvl, **kwargs)

    fld.chief_ray = chief_ray_pkg
    fld.ref_sphere = ref_sphere

    fod = opt_model['analysis_results']['parax_data'].fod
    opd = wave_abr_full_calc(fod, fld, wvl, foc, ray_pkg,
                             chief_ray_pkg, ref_sphere)
    ray, ray_op, wvl = ray_pkg
    return ray, ray_op, wvl, opd


def trace_boundary_rays_at_field(opt_model, fld, wvl, 
                                 use_named_tuples=False, **kwargs):
    """ returns a list of RayPkgs for the boundary rays for field fld
    """
    kwargs['rayerr_filter'] = kwargs.get('rayerr_filter', 'full')
    rim_rays = []
    osp = opt_model.optical_spec
    for p in osp.pupil.pupil_rays:
        ray_pkg, ray_err = trace_ray(opt_model, p, fld, wvl, **kwargs)
        ray, op, wvl = ray_pkg

        if use_named_tuples:
            ray = [RaySeg(*rs) for rs in ray]
        rim_rays.append(RayPkg(ray, op, wvl))
    return rim_rays


def boundary_ray_dict(opt_model, rim_rays):
    pupil_rays = {}
    for ray, lbl in zip(rim_rays, opt_model.optical_spec.pupil.ray_labels):
        pupil_rays[lbl] = ray
    return pupil_rays


def trace_boundary_rays(opt_model, **kwargs):
    rayset = []
    wvl = opt_model.seq_model.central_wavelength()
    fov = opt_model.optical_spec.field_of_view
    for fi, fld in enumerate(fov.fields):
        rim_rays = trace_boundary_rays_at_field(opt_model, fld, wvl, **kwargs)
        fld.pupil_rays = boundary_ray_dict(opt_model, rim_rays)
        rayset.append(rim_rays)
    return rayset


def trace_ray_list_at_field(opt_model, ray_list, fld, wvl, foc, **kwargs):
    """ returns a list of ray |DataFrame| for the ray_list at field fld """
    rayset = []
    for p in ray_list:
        ray_pkg, ray_err = trace_ray(opt_model, p, fld, wvl, **kwargs)
        ray, op, wvl = ray_pkg
        rayset.append(ray)
    rdf_list = [ray_df(r) for r in rayset]
    return rdf_list


def trace_field(opt_model, fld, wvl, foc):
    """ returns a |DataFrame| with the boundary rays for field fld """
    osp = opt_model.optical_spec
    pupil_rays = osp.pupil.pupil_rays
    rdf_list = trace_ray_list_at_field(opt_model, pupil_rays, fld, wvl, foc)
    rset = pd.concat(rdf_list, keys=osp.pupil.ray_labels,
                     names=['pupil'])
    return rset


def trace_all_fields(opt_model):
    """ returns a |DataFrame| with the boundary rays for all fields """
    osp = opt_model.optical_spec
    fld, wvl, foc = osp.lookup_fld_wvl_focus(0)
    fset = []
    for f in osp.field_of_view.fields:
        rset = trace_field(opt_model, f, wvl, foc)
        fset.append(rset)

    fdf = pd.concat(fset, keys=osp.field_of_view.index_labels,
                    names=['field'])
    return fdf


def trace_chief_ray(opt_model, fld, wvl, foc):
    """Trace a chief ray at fld and wvl.

    Returns:
        tuple: **chief_ray**, **cr_exp_seg**
    
            - **chief_ray**: RayPkg of chief ray
            - **cr_exp_seg**: exp_pt, exp_dir, exp_dst, ifc, b4_pt, b4_dir
    """
    fod = opt_model['analysis_results']['parax_data'].fod

    ray_result = trace_safe(opt_model, [0., 0.], fld, wvl,
                            output_filter=None, rayerr_filter='full')

    cr = RayPkg(*ray_result.pkg)

    # cr_exp_pt: E upper bar prime: pupil center for pencils from Q
    # cr_exp_pt, cr_b4_dir, cr_exp_dist
    cr_exp_seg = transfer_to_exit_pupil(opt_model.seq_model.ifcs[-2],
                                        (cr.ray[-2][mc.p],
                                         cr.ray[-2][mc.d]), fod.exp_dist)
    return cr, cr_exp_seg


def trace_fan(opt_model, fan_rng, fld, wvl, foc, img_filter=None,
              **kwargs):
    output_filter = kwargs.pop('output_filter', None)
    rayerr_filter = kwargs.pop('rayerr_filter', None)
    start = np.array(fan_rng[0])
    stop = fan_rng[1]
    num = fan_rng[2]
    step = (stop - start)/(num - 1)
    fan = []
    for r in range(num):
        pupil = np.array(start)
        ray_result = trace_safe(opt_model, pupil, fld, wvl, 
                                output_filter, rayerr_filter, 
                                **kwargs)

        if ray_result.pkg is not None:
            if img_filter:
                result = img_filter(pupil, ray_result.pkg)
                fan.append([pupil, result])
            else:
                fan.append([pupil, ray_result.pkg])

        start += step
    return fan


def trace_grid(opt_model, grid_rng, fld, wvl, foc, img_filter=None,
               form='grid', append_if_none=True, **kwargs):
    output_filter = kwargs.pop('output_filter', None)
    rayerr_filter = kwargs.pop('rayerr_filter', None)
    start = np.array(grid_rng[0])
    stop = grid_rng[1]
    num = grid_rng[2]
    step = np.array((stop - start)/(num - 1))
    grid = []
    for i in range(num):
        if form == 'list':
            working_grid = grid
        elif form == 'grid':
            grid_row = []
            working_grid = grid_row

        for j in range(num):
            pupil = np.array(start)
            ray_result = trace_safe(opt_model, pupil, fld, wvl, 
                                    output_filter, rayerr_filter, 
                                    check_apertures=True, **kwargs)

            if ray_result.pkg is not None:
                if img_filter:
                    result = img_filter(pupil, ray_result.pkg)
                    working_grid.append(result)
                else:
                    working_grid.append([pupil[0], pupil[1], ray_result.pkg])
            else:  # ray outside pupil or failed
                if img_filter:
                    result = img_filter(pupil, None)
                    if result is not None or append_if_none:
                        working_grid.append(result)
                else:
                    if append_if_none:
                        working_grid.append([pupil[0], pupil[1], None])

            start[1] += step[1]
        if form == 'grid':
            grid.append(grid_row)
        start[0] += step[0]
        start[1] = grid_rng[0][1]
    return np.array(grid)


def setup_pupil_coords(opt_model, fld, wvl, foc, 
                       image_pt=None, image_delta=None):
    """Trace chief ray and setup reference sphere for `fld`.

    Returns:
        tuple: **ref_sphere**, **chief_ray_pkg**

            - **ref_sphere**: image_pt, ref_dir, ref_sphere_radius, lcl_tfrm_last
            - **chief_ray_pkg**: chief_ray, cr_exp_seg
    """
    chief_ray_pkg = get_chief_ray_pkg(opt_model, fld, wvl, foc)
    image_pt_2d = None if image_pt is None else image_pt[:2]
    ref_sphere = calculate_reference_sphere(opt_model, fld, wvl, foc,
                                            chief_ray_pkg,
                                            image_pt_2d=image_pt_2d,
                                            image_delta=image_delta)
    return ref_sphere, chief_ray_pkg


def aim_chief_ray(opt_model, fld, wvl=None):
    """ aim chief ray at center of stop surface and save results on **fld** """
    seq_model = opt_model['seq_model']
    osp = opt_model['optical_spec']
    if wvl is None:
        wvl = seq_model.central_wavelength()
    stop = seq_model.stop_surface
    if osp['fov'].is_wide_angle:
        aim_info, rr = find_real_enp(opt_model, stop, fld, wvl)

    else:
        aim_info, cr_ray = iterate_ray(opt_model, stop, np.array([0., 0.]),
                                       fld, wvl)
    return aim_info


def apply_paraxial_vignetting(opt_model):
    fov = opt_model.optical_spec.field_of_view
    pm = opt_model.parax_model
    max_field, jth = fov.max_field()
    for j, fld in enumerate(fov.fields):
        rel_fov = math.sqrt(fld.x**2 + fld.y**2)
        if not fov.is_relative and max_field != 0:
            rel_fov = rel_fov/max_field
        min_vly, min_vuy = pm.paraxial_vignetting(rel_fov)
        if min_vly[1] is not None:
            fld.vly = 1 - min_vly[0]
        if min_vuy[1] is not None:
            fld.vuy = 1 - min_vuy[0]
        # print("Field {:2d}: {:8.3f}, ly:{:8.3f} uy:{:8.3f}".format(
        #     j, rel_fov, fld.vly, fld.vuy))
        

def get_chief_ray_pkg(opt_model, fld, wvl, foc):
    """Get the chief ray package at **fld**, computing it if necessary.

    Args:
        opt_model: :class:`~.OpticalModel` instance
        fld: :class:`~.Field` point for wave aberration calculation
        wvl: wavelength of ray (nm)
        foc: defocus amount

    Returns:
        tuple: **chief_ray**, **cr_exp_seg**

            - **chief_ray**: chief_ray, chief_ray_op, wvl
            - **cr_exp_seg**: chief ray exit pupil segment (pt, dir, dist)

                - pt: chief ray intersection with exit pupil plane
                - dir: direction cosine of the chief ray in exit pupil space
                - dist: distance from interface to the exit pupil point

    """
    if fld.chief_ray is None:
        fld.aim_info = aim_chief_ray(opt_model, fld, wvl=wvl)
        chief_ray_pkg = trace_chief_ray(opt_model, fld, wvl, foc)
    elif fld.chief_ray[0][2] != wvl:
        chief_ray_pkg = trace_chief_ray(opt_model, fld, wvl, foc)
    else:
        chief_ray_pkg = fld.chief_ray
    return chief_ray_pkg


def refocus(opt_model):
    """ Compute a focus shift bringing the axial marginal ray to zero. """
    osp = opt_model['optical_spec']

    fld = osp['fov'].fields[0]      # assumed to be the axial field
    wvl = osp['wvls'].central_wvl

    ray_result = trace_safe(opt_model, [0., 1.], fld, wvl,
                            output_filter=None, rayerr_filter='full', 
                            use_named_tuples=True)

    df_ray, ray_op, wvl = ray_result.pkg

    defocus = -df_ray[-1].p[1]/(df_ray[-2].d[1]/df_ray[-2].d[2])

    return defocus


def trace_astigmatism_coddington_fan(opt_model, fld, wvl, foc):
    """ calculate astigmatism by Coddington trace at **fld** """
    cr_ray_pkg, ray_err = trace_ray(opt_model, [0., 0.], fld, wvl)
    s_dfoc, t_dfoc = trace_coddington_fan(opt_model, cr_ray_pkg, foc=foc)
    return s_dfoc, t_dfoc


def trace_coddington_fan(opt_model, ray_pkg, foc=None):
    """ astigmatism calculation via Coddington trace

.. note:: spherical surfaces only
    """
    seq_model = opt_model.seq_model
    wl = seq_model.index_for_wavelength(ray_pkg.wvl)

    path = itertools.zip_longest(ray_pkg.ray, seq_model.ifcs,
                                 [n[wl] for n in seq_model.rndx],
                                 seq_model.lcl_tfrms,
                                 seq_model.z_dir)

    before_rind = seq_model.rndx[0][wl]
    before_dir = None
    s_before, t_before = None, None
    for r, ifc, after_rind, tfrm, z_dir in path:
        after_rind = after_rind if after_rind is not None else before_rind
        pt, after_dir, after_dst, normal = r
        if before_dir is not None:
            normal_len = norm(normal)
            cosI_prime = np.dot(after_dir, normal)/normal_len
            sinI_prime = math.sqrt(1.0 - cosI_prime**2)
            sinI = after_rind*sinI_prime/before_rind
            cosI = math.sqrt(1.0 - sinI**2)

            obl_power = ifc.optical_power
            if obl_power != 0.0:
                obl_power *= ((after_rind*cosI_prime - before_rind*cosI) /
                              (after_rind - before_rind))
#                print("pwr, obl_pwr, after_dst:",
#                      ifc.optical_power/(after_rind - before_rind),
#                      obl_power, after_dst)
#            else:
#                print("pwr, obl_pwr, after_dst:",
#                      ifc.optical_power, obl_power, after_dst)

            n_by_s_prime = before_rind/s_before + obl_power
            s_prime = after_rind/n_by_s_prime
#            print("s, s':", s_before, s_prime)
            s_before = s_prime - after_dst

            n_cosIp2_by_t_prime = before_rind*cosI**2/t_before + obl_power
            t_prime = after_rind*cosI_prime**2/n_cosIp2_by_t_prime
#            print("t, t':", t_before, t_prime)
            t_before = t_prime - after_dst
        else:
            s_before = -after_dst
            t_before = -after_dst

        before_rind = after_rind
        before_dir = after_dir

    s_dfoc = s_prime*after_dir[2] + pt[2]
    t_dfoc = t_prime*after_dir[2] + pt[2]
    if foc is not None:
        focus_shift = foc
        s_dfoc -= focus_shift
        t_dfoc -= focus_shift

#    print("delta s, t:", s_dfoc, t_dfoc)
    return s_dfoc, t_dfoc


def intersect_2_lines(P1, V1, P2, V2):
    """ intersect 2 non-parallel lines, returning distance from P1

    s = ((P2 - P1) x V1).(V1 x V2)/\|(V1 x V2)\|**2

    `Weisstein, Eric W. "Line-Line Intersection." From MathWorld--A Wolfram Web
    Resource. <http://mathworld.wolfram.com/Line-LineIntersection.html>`_
    """
    Vx = np.cross(V1, V2)
    s = np.dot(np.cross(P2 - P1, V1), Vx)/np.dot(Vx, Vx)
    return s


def trace_astigmatism_curve(opt_model, num_points=21, **kwargs):
    """ Trace a fan of fields and collect astigmatism data for each field.

    Args:
        opt_model: the optical model
        num_points: the number of FOV sampling points
        kwargs: keyword args for ray trace

    Returns:
        tuple: field point, sagittal and tangential focus shifts
    """
    from rayoptics.raytr.opticalspec import Field
    s_data = []
    t_data = []
    field_data = []

    osp = opt_model.optical_spec
    _, wvl, foc = osp.lookup_fld_wvl_focus(0)
    fld = Field()
    max_field = osp['fov'].max_field()[0]
    for f in np.linspace(0., max_field, num=num_points):
        fld.y = f
        ref_sphere, cr_pkg = setup_pupil_coords(opt_model, fld, wvl, foc)
        fld.chief_ray = cr_pkg
        fld.ref_sphere = ref_sphere
        
        s_foc, t_foc = trace_astigmatism(opt_model, fld, wvl, foc, **kwargs)
        s_data.append(s_foc)
        t_data.append(t_foc)
        field_data.append(f)
    return field_data, s_data, t_data


def trace_astigmatism(opt_model, fld, wvl, foc, dx=0.001, dy=0.001):
    """ calculate astigmatism by tracing close rays about the chief ray at **fld**

    This function implicitly assumes that the **fld** point is on a plane of
    symmetry, i.e. the system is rotationally symmetric, bilaterally symmetric,
    or quad symmetric. No check is done to ensure this.

    Args:
        opt_model: the optical model
        fld: a Field object
        wvl: wavelength in nm
        foc: defocus amount
        dx: delta in pupil coordinates for x/sagittal direction
        dy: delta in pupil coordinates for y/tangential direction

    Returns:
        tuple: sagittal and tangential focus shifts at **fld**
    """
    rlist = []
    rlist.append(trace_ray(opt_model, [0., 0.], fld, wvl)[0])
    rlist.append(trace_ray(opt_model, [dx, 0.], fld, wvl)[0])
    rlist.append(trace_ray(opt_model, [0., dy], fld, wvl)[0])
    rlist.append(trace_ray(opt_model, [-dx, 0.], fld, wvl)[0])
    rlist.append(trace_ray(opt_model, [0., -dy], fld, wvl)[0])

    s = intersect_2_lines(rlist[1].ray[-1][mc.p], rlist[1].ray[-1][mc.d],
                          rlist[3].ray[-1][mc.p], rlist[3].ray[-1][mc.d])
    s_foc = s * rlist[1].ray[-1][mc.d][2]

    t = intersect_2_lines(rlist[2].ray[-1][mc.p], rlist[2].ray[-1][mc.d],
                          rlist[4].ray[-1][mc.p], rlist[4].ray[-1][mc.d])
    t_foc = t * rlist[2].ray[-1][mc.d][2]

    if foc is not None:
        focus_shift = foc
        s_foc -= focus_shift
        t_foc -= focus_shift
    return s_foc, t_foc


def iterate_ray_raw(pthlist, ifcx, xy_target, pt0, d0, 
                    obj2pup_dist, eprad, wvl, not_wa, **kwargs):
    """ iterates a ray to xy_target on interface ifcx, returns aim points on
    the paraxial entrance pupil plane

    If idcx is None, i.e. a floating stop surface, returns xy_target.

    If the iteration fails, a TraceError will be raised
    """
    rr = None
    def y_stop_coordinate(y1, *args):
        nonlocal rr
        pthlist, ifcx, pt0, obj2pup_dist, wvl, y_target, not_wa = args
        pt1 = np.array([0., y1, obj2pup_dist])
        dir0 = normalize(pt1 - pt0)
        # handle case where entrance pupil is behind the object (issue #120)
        if not_wa and dir0[2] * pthlist[0][mc.Zdir] < 0:
            dir0 = -dir0

        try:
            ray_pkg = RayPkg(*rt.trace_raw(iter(pthlist), pt0, dir0, wvl))

        except TraceError as ray_error:
            logger.debug(f'ray_error: "{type(ray_error).__name__}", '
                        f'{ray_error.surf=}')
            ray_pkg = RayPkg(*ray_error.ray_pkg)
            rr = RayResult(ray_pkg, ray_error)
            final_coord = np.array([0., 0., 0.])
            if ray_error.surf < ifcx:
                raise ray_error
 
        else:
            rr = RayResult(ray_pkg, None)
            final_coord = ray_pkg.ray[ifcx][mc.p]

        y_ray = final_coord[1]
        return y_ray - y_target

    def surface_coordinate(coord, *args):
        nonlocal rr
        pthlist, ifcx, pt0, dist, wvl, target, not_wa = args
        pt1 = np.array([coord[0], coord[1], dist])
        dir0 = normalize(pt1 - pt0)
        if not_wa and dir0[2] * pthlist[0][mc.Zdir] < 0:
            dir0 = -dir0
        
        try:
            ray_pkg = RayPkg(*rt.trace_raw(iter(pthlist), pt0, dir0, wvl))
            
        except TraceError as ray_error:
            ray_pkg = RayPkg(*ray_error.ray_pkg)
            rr = RayResult(ray_pkg, ray_error)
            final_coord = np.array([0., 0., 0.])
            if ray_error.surf < ifcx:
                raise ray_error
        else:
            rr = RayResult(ray_pkg, None)
            final_coord = ray_pkg.ray[ifcx][mc.p]
        xy_ray = np.array([final_coord[0], final_coord[1]])
        return xy_ray - target

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if ifcx is not None:
            if pt0[0] == 0.0 and xy_target[0] == 0.0:
                # do 1D iteration if field and target points are zero in x
                y_target = xy_target[1]
                try:
                    start_y, results = newton(y_stop_coordinate, 0.,
                                              args=(pthlist, ifcx, 
                                                    pt0, obj2pup_dist, wvl, 
                                                    y_target, not_wa),
                                              disp=False, full_output=True)
                except RuntimeError as rte:
                    # if we come here, start_y is a RuntimeResults object
                    # print(rte)
                    start_y = results.root
                except TraceError:
                    start_y = 0.0
                start_coords = np.array([0., start_y])
            else:
                # do 2D iteration. epsfcn is a parameter increment,
                #  make proportional to pupil radius
                try:
                    start_coords = fsolve(surface_coordinate, 
                                          np.array([0., 0.]),
                                          epsfcn=0.0001*eprad,
                                          args=(pthlist, ifcx, pt0, 
                                                obj2pup_dist, wvl, 
                                                xy_target, not_wa))
                except TraceError:
                    start_coords = np.array([0., 0.])
        else:  # floating stop surface - use entrance pupil for aiming
            start_coords = np.array([0., 0.]) + xy_target

    return start_coords, rr
