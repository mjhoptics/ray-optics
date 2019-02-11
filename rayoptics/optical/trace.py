#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Supports model ray tracing in terms of relative aperture and field.

.. Created on Mon Sep 17 23:10:59 2018

.. codeauthor: Michael J. Hayford
"""

import itertools
import math
import numpy as np
from numpy.linalg import norm
from collections import namedtuple
import pandas as pd
import attr

from . import raytrace as rt
from . import model_constants as mc
from rayoptics.util.misc_math import normalize


RayPkg = namedtuple('RayPkg', ['ray', 'op', 'wvl'])


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
class RaySeg():
    inc_pt = attr.ib()
    after_dir = attr.ib()
    after_dst = attr.ib()
    phase = attr.ib(default=0.0)


def list_ray(ray):
    print("            X            Y            Z           L"
          "            M            N               Len")
    r = None
    for i, r in enumerate(ray):
        print("{:3d}: {:12.5f} {:12.5f} {:12.5f} {:12.6f} {:12.6f} "
              "{:12.6f} {:12.5g}".format(i,
                                         r[mc.p][0], r[mc.p][1], r[mc.p][2],
                                         r[mc.d][0], r[mc.d][1], r[mc.d][2],
                                         r[mc.dst]))


def trace(sequence, pt0, dir0, wvl, **kwargs):
    """ returns (ray, ray_opl, wvl)

    Args:
        sequence: a Sequence or generator that returns a list containing:
            Intfc, Gap, Index, Trfm, Z_Dir
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
    return rt.trace(sequence, pt0, dir0, wvl, **kwargs)


def trace_base(opt_model, pupil, fld, wvl, **kwargs):
    """ returns (ray, ray_opl, wvl) """
    vig_pupil = fld.apply_vignetting(pupil)
    osp = opt_model.optical_spec
    fod = osp.parax_data.fod
    eprad = fod.enp_radius
    pt1 = np.array([eprad*vig_pupil[0], eprad*vig_pupil[1],
                    fod.obj_dist+fod.enp_dist])
    pt0 = osp.obj_coords(fld)
    dir0 = pt1 - pt0
    length = norm(dir0)
    dir0 = dir0/length
    return rt.trace(opt_model.seq_model, pt0, dir0, wvl, **kwargs)


def trace_with_opd(opt_model, pupil, fld, wvl, foc, **kwargs):
    """ returns (ray, ray_opl, wvl, opd) """
    ray_pkg = trace_base(opt_model, pupil, fld, wvl, **kwargs)

    rs_pkg, cr_pkg = setup_pupil_coords(opt_model, fld, wvl, foc)
    fld.chief_ray = cr_pkg
    fld.ref_sphere = rs_pkg

    opd_pkg = rt.wave_abr(fld, wvl, ray_pkg)
    ray, ray_op, wvl = ray_pkg
    return ray, ray_op, wvl, opd_pkg[0]


def trace_boundary_rays_at_field(opt_model, fld, wvl):
    """ returns a list of (ray, opl, wvl) for the boundary rays
        for field fld
    """
    rim_rays = []
    osp = opt_model.optical_spec
    for p in osp.pupil.pupil_rays:
        ray, op, wvl = trace_base(opt_model, p, fld, wvl)
        rim_rays.append((ray, op, wvl))
    return rim_rays


def trace_ray_list_at_field(opt_model, ray_list, fld, wvl, foc):
    rayset = pd.DataFrame(data=np.nan)
    for p in ray_list:
        ray, op, wvl = trace_base(opt_model, p, fld, wvl)
        rayset[(fld, wvl, foc, p)] = ray
    return rayset


def trace_field(opt_model, fld, wvl):
    """ returns a |DataFrame| with the boundary rays for field fld """
    osp = opt_model.optical_spec
    rayset = trace_boundary_rays_at_field(opt_model, fld, wvl)
    rdf_list = [ray_df(r[0]) for r in rayset]
    rset = pd.concat(rdf_list, keys=osp.pupil.ray_labels,
                     names=['pupil'])
    return rset


def trace_all_fields(opt_model):
    """ returns a |DataFrame| with the boundary rays for all fields """
    osp = opt_model.optical_spec
    fld, wvl, foc = osp.lookup_fld_wvl_focus(0)
    fset = []
    for f in osp.field_of_view.fields:
        rset = trace_field(opt_model, f, wvl)
        fset.append(rset)

    fdf = pd.concat(fset, keys=osp.field_of_view.index_labels,
                    names=['field'])
    return fdf


def trace_chief_ray(opt_model, fld, wvl, foc):
    osp = opt_model.optical_spec
    fod = osp.parax_data.fod

    ray, op, wvl = trace_base(opt_model, [0., 0.], fld, wvl)
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
        ray_pkg = trace_base(opt_model, pupil, fld, wvl, **kwargs)

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


def setup_pupil_coords(opt_model, fld, wvl, foc,
                       chief_ray_pkg=None, image_pt=None):
    if chief_ray_pkg is None:
        chief_ray_pkg = trace_chief_ray(opt_model, fld, wvl, foc)
    elif chief_ray_pkg[2] != wvl:
        chief_ray_pkg = trace_chief_ray(opt_model, fld, wvl, foc)

    cr, cr_exp_seg = chief_ray_pkg

    if image_pt is None:
        image_pt = cr.ray[-1][mc.p]

    # cr_exp_pt: E upper bar prime: pupil center for pencils from Q
    # cr_exp_pt, cr_b4_dir, cr_dst
    cr_exp_pt = cr_exp_seg[mc.p]
    cr_exp_dist = cr_exp_seg[mc.dst]

    seq_model = opt_model.seq_model
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
    n_obj = seq_model.rndx[wvl].iloc[0]
    n_img = seq_model.rndx[wvl].iloc[-1]
    ref_sphere_pkg = (ref_sphere, opt_model.optical_spec.parax_data,
                      n_obj, n_img, z_dir)

    return ref_sphere_pkg, chief_ray_pkg


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
    n_obj = seq_model.rndx[wvl].iloc[0]
    n_img = seq_model.rndx[wvl].iloc[-1]
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
    wvl = ray_pkg.wvl

    path = itertools.zip_longest(ray_pkg.ray, seq_model.ifcs,
                                 seq_model.rndx[wvl], seq_model.lcl_tfrms,
                                 seq_model.z_dir)

    before_rind = seq_model.rndx[wvl][0]
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
        s focus shift, t focus shift: sagittal and tangential focus shifts at **fld**
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
