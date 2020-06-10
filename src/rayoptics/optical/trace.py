#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Supports model ray tracing in terms of relative aperture and field.

.. Created on Mon Sep 17 23:10:59 2018

.. codeauthor: Michael J. Hayford
"""

import itertools
import logging
import math
import numpy as np
from numpy.linalg import norm
from scipy.optimize import newton, fsolve
from collections import namedtuple
import pandas as pd
import attr

from . import raytrace as rt
from rayoptics.optical.analyses import (wave_abr_full_calc,
                                        get_chief_ray_pkg,
                                        setup_exit_pupil_coords)
from . import model_constants as mc
from .traceerror import TraceError, TraceMissedSurfaceError, TraceTIRError
from rayoptics.util.misc_math import normalize


RayPkg = namedtuple('RayPkg', ['ray', 'op', 'wvl'])
RaySeg = namedtuple('RaySeg', ['p', 'd', 'dst', 'nrml'])


def ray_pkg(ray_pkg):
    """ return a |Series| containing a ray package (RayPkg) """
    return pd.Series(ray_pkg, index=['ray', 'op', 'wvl'])


def ray_df(ray):
    """ return a |DataFrame| containing ray data """
    r = pd.DataFrame(ray, columns=['inc_pt', 'after_dir',
                                   'after_dst', 'normal'])
    r.index.names = ['intrfc']
    return r


@attr.s
class RSeg():
    inc_pt = attr.ib()
    after_dir = attr.ib()
    after_dst = attr.ib()
    normal = attr.ib()
    phase = attr.ib(default=0.0)


def list_ray(ray_obj, tfrms=None):
    """ pretty print a ray either in local or global coordinates """
    if isinstance(ray_obj, tuple):
        ray = ray_obj[0]
    else:
        ray = ray_obj

    colHeader = "            X            Y            Z           L" \
                "            M            N               Len"
    print(colHeader)

    colFormats = "{:3d}: {:12.5f} {:12.5f} {:12.5g} {:12.6f} {:12.6f} " \
                 "{:12.6f} {:12.5g}"

    for i, r in enumerate(ray):
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


def trace(seq_model, pt0, dir0, wvl, **kwargs):
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
    return rt.trace(seq_model, pt0, dir0, wvl, **kwargs)


def trace_base(opt_model, pupil, fld, wvl, **kwargs):
    """Trace ray specified by relative aperture and field point.

    Args:
        opt_model: instance of :class:`~.OpticalModel` to trace
        pupil: relative pupil coordinates of ray
        fld: instance of :class:`~.Field`
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
    vig_pupil = fld.apply_vignetting(pupil)
    osp = opt_model.optical_spec
    fod = osp.parax_data.fod
    eprad = fod.enp_radius
    aim_pt = np.array([0., 0.])
    if hasattr(fld, 'aim_pt') and fld.aim_pt is not None:
        aim_pt = fld.aim_pt
    pt1 = np.array([eprad*vig_pupil[0]+aim_pt[0], eprad*vig_pupil[1]+aim_pt[1],
                    fod.obj_dist+fod.enp_dist])
    pt0 = osp.obj_coords(fld)
    dir0 = pt1 - pt0
    length = norm(dir0)
    dir0 = dir0/length
    return rt.trace(opt_model.seq_model, pt0, dir0, wvl, **kwargs)


def iterate_ray(opt_model, ifcx, xy_target, fld, wvl, **kwargs):
    """ iterates a ray to xy_target on interface ifcx, returns aim points on
    the paraxial entrance pupil plane

    If idcx is None, i.e. a floating stop surface, returns xy_target.

    If the iteration fails, a TraceError will be raised
    """
    def y_stop_coordinate(y1, *args):
        seq_model, ifcx, pt0, dist, wvl, y_target = args
        pt1 = np.array([0., y1, dist])
        dir0 = pt1 - pt0
        length = norm(dir0)
        dir0 = dir0/length
        try:
            ray, _, _ = rt.trace(seq_model, pt0, dir0, wvl)
        except TraceMissedSurfaceError as ray_miss:
            ray = ray_miss.ray
            if ray_miss.surf <= ifcx:
                raise ray_miss
        except TraceTIRError as ray_tir:
            ray = ray_tir.ray
            if ray_tir.surf < ifcx:
                raise ray_tir
        y_ray = ray[ifcx][mc.p][1]
#        print(y1, y_ray)
        return y_ray - y_target

    def surface_coordinate(coord, *args):
        seq_model, ifcx, pt0, dist, wvl, target = args
        pt1 = np.array([coord[0], coord[1], dist])
        dir0 = pt1 - pt0
        length = norm(dir0)
        dir0 = dir0/length
        ray, _, _ = rt.trace(seq_model, pt0, dir0, wvl)
        xy_ray = np.array([ray[ifcx][mc.p][0], ray[ifcx][mc.p][1]])
#        print(coord[0], coord[1], xy_ray[0], xy_ray[1])
        return xy_ray - target

    seq_model = opt_model.seq_model
    osp = opt_model.optical_spec

    fod = osp.parax_data.fod
    dist = fod.obj_dist + fod.enp_dist

    pt0 = osp.obj_coords(fld)
    if ifcx is not None:
        if pt0[0] == 0.0 and xy_target[0] == 0.0:
            # do 1D iteration if field and target points are zero in x
            y_target = xy_target[1]
            logging.captureWarnings(True)
            try:
                start_y, results = newton(y_stop_coordinate, 0.,
                                          args=(seq_model, ifcx, pt0,
                                                dist, wvl, y_target),
                                          disp=False, full_output=True)
            except RuntimeError as rte:
                # if we come here, start_y is a RuntimeResults object
                print(rte)
                start_y = results.root
            except TraceError:
                start_y = 0.0
            start_coords = np.array([0., start_y])
        else:
            # do 2D iteration. epsfcn is a parameter increment,
            #  make proportional to pupil radius
            try:
                start_coords = fsolve(surface_coordinate, np.array([0., 0.]),
                                      epsfcn=0.0001*fod.enp_radius,
                                      args=(seq_model, ifcx, pt0, dist,
                                            wvl, xy_target))
            except TraceError:
                start_coords = np.array([0., 0.])
    else:  # floating stop surface - use entrance pupil for aiming
        start_coords = np.array([0., 0.]) + xy_target

    return start_coords


def trace_with_opd(opt_model, pupil, fld, wvl, foc, **kwargs):
    """ returns (ray, ray_opl, wvl, opd) """
    chief_ray_pkg = get_chief_ray_pkg(opt_model, fld, wvl, foc)
    image_pt_2d = None if 'image_pt' not in kwargs else kwargs['image_pt'][:2]
    ref_sphere = setup_exit_pupil_coords(opt_model, fld, wvl, foc,
                                         chief_ray_pkg,
                                         image_pt_2d=image_pt_2d)

    ray, op, wvl = trace_base(opt_model, pupil, fld, wvl, **kwargs)
    # opl = rt.calc_optical_path(ray, opt_model.seq_model.path())
    ray_pkg = ray, op, wvl

    fld.chief_ray = chief_ray_pkg
    fld.ref_sphere = ref_sphere

    fod = opt_model.optical_spec.parax_data.fod
    opd = wave_abr_full_calc(fod, fld, wvl, foc, ray_pkg,
                             chief_ray_pkg, ref_sphere)
    ray, ray_op, wvl = ray_pkg
    return ray, ray_op, wvl, opd


def trace_boundary_rays_at_field(opt_model, fld, wvl, use_named_tuples=False):
    """ returns a list of RayPkgs for the boundary rays for field fld
    """
    rim_rays = []
    osp = opt_model.optical_spec
    for p in osp.pupil.pupil_rays:
        try:
            ray, op, wvl = trace_base(opt_model, p, fld, wvl)
        except TraceError as ray_error:
            ray = ray_error.ray
            op = 0.

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


def trace_ray_list_at_field(opt_model, ray_list, fld, wvl, foc):
    """ returns a list of ray |DataFrame| for the ray_list at field fld """
    rayset = []
    for p in ray_list:
        ray, op, wvl = trace_base(opt_model, p, fld, wvl)
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
    """Trace a chief ray for fld and wvl, returning the ray_pkg and exit pupil segment."""
    osp = opt_model.optical_spec
    fod = osp.parax_data.fod

    ray, op, wvl = trace_base(opt_model, [0., 0.], fld, wvl)
    # op = rt.calc_optical_path(ray, opt_model.seq_model.path())
    cr = RayPkg(ray, op, wvl)

    # cr_exp_pt: E upper bar prime: pupil center for pencils from Q
    # cr_exp_pt, cr_b4_dir, cr_exp_dist
    cr_exp_seg = rt.transfer_to_exit_pupil(opt_model.seq_model.ifcs[-2],
                                           (cr.ray[-2][mc.p],
                                            cr.ray[-2][mc.d]), fod.exp_dist)
    return cr, cr_exp_seg


def trace_fan(opt_model, fan_rng, fld, wvl, foc, img_filter=None,
              **kwargs):
    start = np.array(fan_rng[0])
    stop = fan_rng[1]
    num = fan_rng[2]
    step = (stop - start)/(num - 1)
    fan = []
    for r in range(num):
        pupil = np.array(start)
        ray, op, wvl = trace_base(opt_model, pupil, fld, wvl, **kwargs)
        # opl = rt.calc_optical_path(ray, opt_model.seq_model.path())
        ray_pkg = ray, op, wvl

        if img_filter:
            result = img_filter(pupil, ray_pkg)
            fan.append([pupil, result])
        else:
            fan.append([pupil, ray_pkg])

        start += step
    return fan


def trace_grid(opt_model, grid_rng, fld, wvl, foc, img_filter=None,
               form='grid', append_if_none=True, **kwargs):
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
            if (pupil[0]**2 + pupil[1]**2) < 1.0:
                ray_pkg = trace_base(opt_model, pupil, fld, wvl, **kwargs)
                if img_filter:
                    result = img_filter(pupil, ray_pkg)
                    working_grid.append(result)
                else:
                    working_grid.append([pupil[0], pupil[1], ray_pkg])
            else:  # ray outside pupil
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


def aim_chief_ray(opt_model, fld, wvl=None):
    """ aim chief ray at center of stop surface and save results on **fld** """
    seq_model = opt_model.seq_model
    if wvl is None:
        wvl = seq_model.central_wavelength()
    stop = seq_model.stop_surface
    aim_pt = iterate_ray(opt_model, stop, np.array([0., 0.]), fld, wvl)
    return aim_pt


def setup_pupil_coords(opt_model, fld, wvl, foc, image_pt=None):
    chief_ray_pkg = get_chief_ray_pkg(opt_model, fld, wvl, foc)
    image_pt_2d = None if image_pt is None else image_pt[:2]
    ref_sphere = setup_exit_pupil_coords(opt_model, fld, wvl, foc,
                                         chief_ray_pkg,
                                         image_pt_2d=image_pt_2d)
    return ref_sphere, chief_ray_pkg


def setup_canonical_coords(opt_model, fld, wvl, image_pt=None):
    osp = opt_model.optical_spec
    seq_model = opt_model.seq_model
    fod = osp.parax_data.fod

    if fld.chief_ray is None:
        ray, op, wvl = trace_base(opt_model, [0., 0.], fld, wvl)
        fld.chief_ray = RayPkg(ray, op, wvl)
    cr = fld.chief_ray

    if image_pt is None:
        image_pt = cr.ray[-1][mc.p]

    # cr_exp_pt: E upper bar prime: pupil center for pencils from Q
    # cr_exp_pt, cr_b4_dir, cr_dst
    cr_exp_seg = rt.transfer_to_exit_pupil(seq_model.ifcs[-2],
                                           (cr.ray[-2][mc.p],
                                            cr.ray[-2][mc.d]),
                                           fod.exp_dist)
    cr_exp_pt = cr_exp_seg[mc.p]
    cr_exp_dist = cr_exp_seg[mc.dst]

    img_dist = seq_model.gaps[-1].thi
    img_pt = np.array(image_pt)
    img_pt[2] += img_dist

    # R' radius of reference sphere for O'
    ref_sphere_vec = img_pt - cr_exp_pt
    ref_sphere_radius = np.linalg.norm(ref_sphere_vec)
    ref_dir = normalize(ref_sphere_vec)

    ref_sphere = (image_pt, cr_exp_pt, cr_exp_dist,
                  ref_dir, ref_sphere_radius)

    z_dir = seq_model.z_dir[-1]
    wl = seq_model.index_for_wavelength(wvl)
    n_obj = seq_model.rndx[0][wl]
    n_img = seq_model.rndx[-1][wl]
    ref_sphere_pkg = (ref_sphere, osp.parax_data, n_obj, n_img, z_dir)
    fld.ref_sphere = ref_sphere_pkg
    return ref_sphere_pkg, cr


def trace_astigmatism_coddington_fan(opt_model, fld, wvl, foc):
    """ calculate astigmatism by Coddington trace at **fld** """
    cr = RayPkg(*trace_base(opt_model, [0., 0.], fld, wvl))
    s_dfoc, t_dfoc = trace_coddington_fan(opt_model, cr, foc=foc)
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
    for r, ifc, after_rind, tfrm, z_dir in path:
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
    rlist.append(RayPkg(*trace_base(opt_model, [0., 0.], fld, wvl)))
    rlist.append(RayPkg(*trace_base(opt_model, [dx, 0.], fld, wvl)))
    rlist.append(RayPkg(*trace_base(opt_model, [0., dy], fld, wvl)))
    rlist.append(RayPkg(*trace_base(opt_model, [-dx, 0.], fld, wvl)))
    rlist.append(RayPkg(*trace_base(opt_model, [0., -dy], fld, wvl)))

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
