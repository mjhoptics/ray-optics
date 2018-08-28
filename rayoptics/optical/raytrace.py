#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Functions to support ray tracing a sequential optical model

Created on Thu Jan 25 11:01:04 2018

@author: Michael J. Hayford
"""

import itertools
import numpy as np
from collections import namedtuple
from numpy.linalg import norm
from math import sqrt, copysign
from . import transform as trns
from util.misc_math import normalize
#from optical.model_constants import Surf, Gap
import attr

Intfc, Gap, Index, Trfm, Z_Dir = range(5)
pt, dcs = range(2)

RayPkg = namedtuple('RayPkg', ['ray', 'op', 'wvl'])


@attr.s
class RaySeg():
    inc_pt = attr.ib()
    after_dir = attr.ib()
    dst_before = attr.ib()
    phase = attr.ib(default=0.0)


def list_ray(ray):
    print("          X            Y            Z           L"
          "            M            N               Len")
    for i, r in enumerate(ray):
        print("{}: {:12.5f} {:12.5f} {:12.5f} {:12.6f} {:12.6f} "
              "{:12.6f} {:12.5g}".format(i,
                                         r[0][0], r[0][1], r[0][2],
                                         r[1][0], r[1][1], r[1][2], r[2]))


def bend(d_in, normal, n_in, n_out):
    """ refract incoming direction, d_in, about normal """
    normal_len = norm(normal)
    cosI = np.dot(d_in, normal)/normal_len
    sinI_sqr = 1.0 - cosI*cosI
    n_cosIp = copysign(sqrt(n_out*n_out - n_in*n_in*sinI_sqr), cosI)
    alpha = n_cosIp - n_in*cosI
    d_out = (n_in*d_in + alpha*normal)/n_out
    return d_out


def reflect(d_in, normal):
    """ reflect incoming direction, d_in, about normal """
    normal_len = norm(normal)
    cosI = np.dot(d_in, normal)/normal_len
    d_out = d_in - 2.0*cosI*normal
    return d_out


def phase(intrfc, inc_pt, d_in, normal, wvl, n_in, n_out):
    """ apply phase shift to incoming direction, d_in, about normal """
    d_out, dW = intrfc.phase(inc_pt, d_in, normal, wvl)
    return d_out, dW


def trace(seq_model, pt0, dir0, wvl, eps=1.0e-12):
    """ fundamental raytrace function

    inputs:
        seq_model: the sequential model to be traced
        pt0: starting point in coords of first interface
        dir0: starting direction cosines in coords of first interface
        wvl: wavelength in nm
        eps: accuracy tolerance for surface intersection calculation

    returns ray, op_delta
    where ray is:
        [pt, after_dir, dst_before]
        where
        pt: the intersection point of the ray in interface coordinates
        after_dir: the ray direction cosine following the interface in
                   interface coordinates
        dst_before: the geometric distance from the previous interface
    and
        op_delta: optical path wrt equally inclined chords to the optical axis
    """
    path = itertools.zip_longest(seq_model.ifcs, seq_model.gaps,
                                 seq_model.rndx[wvl], seq_model.lcl_tfrms,
                                 seq_model.z_dir)
    path_pkg = (path, seq_model.get_num_surfaces())
    return trace_raw(path_pkg, pt0, dir0, wvl, eps)


def trace_raw(path_pkg, pt0, dir0, wvl, eps=1.0e-12):
    """ fundamental raytrace function

    inputs:
        path_pkg: an iterator containing interfaces and gaps to be traced
        pt0: starting point in coords of first interface
        dir0: starting direction cosines in coords of first interface
        wvl: wavelength in nm
        eps: accuracy tolerance for surface intersection calculation

    returns ray, op_delta, wvl
    where ray is:
        [pt, after_dir, dst_b4]
        where
        pt: the intersection point of the ray in interface coordinates
        after_dir: the ray direction cosine following the interface in
                   interface coordinates
        dst_b4: the geometric distance from the previous interface
    and
        op_delta: optical path wrt equally inclined chords to the optical axis
    """
    ray = []
    eic = []

    path, path_length = path_pkg

    # trace object surface
    obj = next(path)
    srf_obj = obj[Intfc]
    dst_b4, pt_obj = srf_obj.intersect(pt0, dir0)
    ray.append([pt_obj, dir0, dst_b4])

    before = obj
    before_pt = pt_obj
    before_dir = dir0
    tfrm_from_before = before[Trfm]
    z_dir_before = before[Z_Dir]
    n_before = before[Index] if z_dir_before > 0.0 else -before[Index]

    op_delta = 0.0
    surf = 0
    # loop of remaining surfaces in path
    while True:
        try:
            after = next(path)

            rt, t = tfrm_from_before
            b4_pt, b4_dir = rt.dot(before_pt - t), rt.dot(before_dir)

            pp_dst = -b4_pt.dot(b4_dir)
            pp_pt_before = b4_pt + pp_dst*b4_dir

            ifc = after[Intfc]
            z_dir_after = after[Z_Dir]
            n_after = after[Index] if z_dir_after > 0.0 else -after[Index]

            # intersect ray with profile
            pp_dst_intrsct, inc_pt = ifc.intersect(pp_pt_before, b4_dir, eps)
            normal = ifc.normal(inc_pt)

            eic_dst_before = ((inc_pt.dot(b4_dir) + z_dir_before*inc_pt[2]) /
                              (1.0 + z_dir_before*b4_dir[2]))

            # refract or reflect ray at interface
            if ifc.refract_mode == 'REFL':
                after_dir = reflect(b4_dir, normal)
            elif ifc.refract_mode == 'PHASE':
                after_dir, phs = phase(ifc, inc_pt, b4_dir, normal, wvl,
                                       n_before, n_after)
                op_delta += phs
            else:
                after_dir = bend(b4_dir, normal, n_before, n_after)

            eic_dst_after = ((inc_pt.dot(after_dir) + z_dir_after*inc_pt[2]) /
                             (1.0 + z_dir_after*after_dir[2]))

            surf += 1

            dW = n_after*eic_dst_after - n_before*eic_dst_before
            eic.append([n_before, eic_dst_before,
                        n_after, eic_dst_after, dW])

            dst_b4 = pp_dst + pp_dst_intrsct
            ray.append([inc_pt, after_dir, dst_b4])
#            print("after:", surf, inc_pt, after_dir)
#            print("e{}= {:12.5g} e{}'= {:12.5g} dW={:10.8g} n={:8.5g}"
#                  " n'={:8.5g}".format(surf, eic_dst_before,
#                                       surf, eic_dst_after,
#                                       dW, before[Index], after[Index]))
            before_pt = inc_pt
            before_dir = after_dir
            n_before = n_after
            z_dir_before = z_dir_after
            before = after
            tfrm_from_before = before[Trfm]

        except StopIteration:
            P, P1k, Ps = calc_path_length(eic, offset=1)
            op_delta += P
            break

    return ray, op_delta, wvl


def eic_distance(r, r0):
    e = (np.dot(r[dcs] + r0[dcs], r[pt] - r0[pt]) /
         (1. + np.dot(r[dcs], r0[dcs])))
    return e


def wave_abr(seq_model, fld, wvl, ray_pkg):
    return wave_abr_real_coord(seq_model, fld, wvl, ray_pkg)
#    return wave_abr_HHH(fld.ref_sphere_pkg, fld.chief_ray_pkg, ray_pkg)


def wave_abr_real_coord(seq_model, fld, wvl, ray_pkg):
    ref_sphere, parax_data, n_obj, n_img, z_dir = fld.ref_sphere
    image_pt, cr_exp_pt, cr_exp_dist, ref_dir, ref_sphere_radius = ref_sphere
    chief_ray, chief_ray_op, wvl = fld.chief_ray[0]
    ray, ray_op, wvl = ray_pkg
    fod = parax_data[2]
    k = -2  # last interface in sequence

    # eq 3.12
    e1 = eic_distance((ray[1][pt], ray[0][dcs]),
                      (chief_ray[1][pt], chief_ray[0][dcs]))
    # eq 3.13
    ekp = eic_distance((ray[k][pt], ray[k][dcs]),
                       (chief_ray[k][pt], chief_ray[k][dcs]))

    dst = ekp - cr_exp_dist

    eic_exp_pt = ray[k][pt] - dst*ray[k][dcs]
#    eic_exp_pt[2] -= cr_exp_dist
    p_coord = eic_exp_pt - cr_exp_pt
    F = ref_dir.dot(ray[k][dcs]) - ray[k][dcs].dot(p_coord)/ref_sphere_radius
    J = p_coord.dot(p_coord)/ref_sphere_radius - 2.0*ref_dir.dot(p_coord)
    ep = J/(F + sqrt(F**2 + J/ref_sphere_radius))

    opd = -n_obj*e1 - ray_op + n_img*ekp + chief_ray_op - n_img*ep
    return opd, e1, ekp, ep


def wave_abr_HHH(seq_model, fld, wvl, ray_pkg):
    ref_sphere, parax_data, n_obj, n_img, z_dir = fld.ref_sphere
    image_pt, cr_exp_pt, ref_dir, ref_sphere_radius = ref_sphere
    chief_ray, chief_ray_op, wvl = fld.chief_ray
    ray, ray_op, wvl = ray_pkg
    ax_ray, pr_ray, fod = parax_data
    k = -2  # last interface in sequential model
    ax_k = ax_ray[k]
    pr_k = pr_ray[k]
    ht, slp = range(2)
    H = n_img*(pr_k[ht]*ax_k[slp] - ax_k[ht]*pr_k[slp])

    # eq 3.12
    e1 = eic_distance((ray[1][pt], ray[0][dcs]),
                      (chief_ray[1][pt], chief_ray[0][dcs]))
    # eq 3.13
    ekp = eic_distance((ray[k][pt], ray[k][dcs]),
                       (chief_ray[k][pt], chief_ray[k][dcs]))

    # eq 4.33
    eic_pt = ray[k][pt] - ekp*ray[k][dcs]
    print("eic_pt", eic_pt)

    Nk_cr = chief_ray[k][dcs][2]
    Zk_cr = chief_ray[k][pt][2]

    def reduced_pupil_coord(X, L, e):
        coef1 = -n_img/(H * Nk_cr)
        xp = coef1*(pr_k[slp]*(Nk_cr*(X - L*e) - L*Zk_cr) +
                    pr_k[ht]*L)
        return xp

    # eq 5.4
    xp_ray = reduced_pupil_coord(ray[k][pt][0], ray[k][dcs][0], ekp)
    yp_ray = reduced_pupil_coord(ray[k][pt][1], ray[k][dcs][1], ekp)
    # eq 5.5
    xp_cr = reduced_pupil_coord(chief_ray[k][pt][0], chief_ray[k][dcs][0], 0.)
    yp_cr = reduced_pupil_coord(chief_ray[k][pt][1], chief_ray[k][dcs][1], 0.)
    # eq 5.6
    zp_ray = -(((ray[k][dcs][0] + chief_ray[k][dcs][0])*(xp_ray - xp_cr) +
                (ray[k][dcs][1] + chief_ray[k][dcs][1])*(yp_ray - yp_cr)) /
                (ray[k][dcs][2] + chief_ray[k][dcs][2]))

    rpc_ray = np.array([xp_ray, yp_ray, zp_ray])
    rpc_cr = np.array([xp_cr, yp_cr, 0.])
    print("rpc ray", xp_ray, yp_ray, zp_ray)
    print("rpc cr", xp_cr, yp_cr, 0.)
    # eq 4.11
    G0_ref = n_img*ax_k[slp]*image_pt[0]
    H0_ref = n_img*ax_k[slp]*image_pt[1]

    def reduced_image_coord(X, L, Z, N):
        G0 = (n_img/N)*(ax_k[slp]*(N*X - L*Z) + ax_k[ht]*L)
        return G0

    # eq 5.13
    G0_ray = reduced_image_coord(ray[k][pt][0], ray[k][dcs][0],
                                 ray[k][pt][2], ray[k][dcs][2])
    H0_ray = reduced_image_coord(ray[k][pt][1], ray[k][dcs][1],
                                 ray[k][pt][2], ray[k][dcs][2])
    # eq 5.14
    G0_cr = reduced_image_coord(chief_ray[k][pt][0], chief_ray[k][dcs][0],
                                Zk_cr, Nk_cr)
    H0_cr = reduced_image_coord(chief_ray[k][pt][1], chief_ray[k][dcs][1],
                                Zk_cr, Nk_cr)
    print("G0, H0_ref; G0, H0_cr:", G0_ref, H0_ref, G0_cr, H0_cr)
    # eq 4.17
    a = pr_k[slp]*G0_ref/H + ax_k[slp]*xp_cr
    b = pr_k[slp]*H0_ref/H + ax_k[slp]*yp_cr
    g = z_dir/sqrt(1. - a**2 + b**2)
    # eq 4.18
    ref_dir = np.array([-a*g, -b*g, g])
    print("ref_dir, cr_dir", ref_dir, chief_ray[k][dcs])
    # eq 4.25
    F = (np.dot(ref_dir, ray[k][dcs]) + ref_dir[2]*ax_k[slp] *
         np.dot(chief_ray[k][dcs], (rpc_ray - rpc_cr)))

    # eq 4.28
    Ja = (ref_dir[2]*ax_k[slp] *
          np.dot((eic_pt - chief_ray[k][0]), (rpc_ray - rpc_cr)))
    Jb = -(2.0*(ref_dir[2]/Nk_cr)*(ax_k[ht] - ax_k[slp]*Zk_cr) *
           np.dot(chief_ray[k][dcs], (rpc_ray - rpc_cr)))
    Jc = (2.0*(ref_dir[2]/n_img)*((G0_cr - G0_ref)*(xp_ray - xp_cr) +
                                  (H0_cr - H0_ref)*(yp_ray - yp_cr)))
    J = Ja + Jb + Jc
    print("F, J, Ja, Jb, Jc", F, J, Ja, Jb, Jc)
#    J = ((ref_dir[2]*ax_k[slp] *
#         np.dot((eic_pt - chief_ray[k][0]), (rpc_ray - rpc_cr))) -
#         2.0*(ref_dir[2]/Nk_cr)*(ax_k[ht] - ax_k[slp]*Zk_cr) *
#         np.dot(chief_ray[k][dcs], (rpc_ray - rpc_cr)) +
#         2.0*(ref_dir[2]/n_img)*((G0_cr - G0_ref)*(xp_ray - xp_cr) +
#                                 (H0_cr - H0_ref)*(yp_ray - yp_cr)))

#    # eq 4.29 Q' = image_pt
#    F = (np.dot(chief_ray[k][dcs], ray[k][dcs]) + Nk_cr*ax_k[slp] *
#         np.dot(chief_ray[k][dcs], (rpc_ray - rpc_cr)))
#    J = ((Nk_cr*ax_k[slp] *
#         np.dot((eic_pt - chief_ray[k][0]), (rpc_ray - rpc_cr))) -
#         2.0*(ax_k[ht] - ax_k[slp]*Zk_cr) *
#         np.dot(chief_ray[k][dcs], (rpc_ray - rpc_cr)))

    # eq 4.21
    ep = J/(F + sqrt(F**2 + J*(n_img*ax_k[slp]*pr_k[slp]*ref_dir[2] / H)))

    print("F, J, ep (canon)", F, J, ep)

    # eq 3.14/3.22, using 3.15
    opd = -n_obj*e1 - ray_op + n_img*ekp + chief_ray_op - n_img*ep
    return opd, e1, ekp, ep


def transfer_to_exit_pupil(interface, ray_seg, exp_dst_parax):
    if interface.decenter:
        # get transformation info after surf
        r, t = interface.decenter.tform_after_surf()
        rt = r.transpose()
        b4_pt, b4_dir = rt.dot(ray_seg[0] - t), rt.dot(ray_seg[1])
    else:
        b4_pt, b4_dir = ray_seg[0], ray_seg[1]

    h = b4_pt[pt]**2 + b4_pt[dcs]**2
    u = b4_dir[pt]**2 + b4_dir[dcs]**2
    if u == 0.0:
        dst = exp_dst_parax
    else:
        dst = -sqrt(h/u)

    exp_pt = b4_pt + dst*b4_dir

    return exp_pt, b4_dir, dst


def eic_path_accumulation(ray, rndx, lcl_tfrms, z_dir):
    """ fundamental raytrace function

    inputs:
        path: an iterator containing interfaces and gaps to be traced
        pt0: starting point in coords of first interface
        dir0: starting direction cosines in coords of first interface
        wl: wavelength in nm
        eps: accuracy tolerance for surface intersection calculation

    returns ray, op_delta
    where ray is:
        [pt, after_dir, dst_b4]
        where
        pt: the intersection point of the ray in interface coordinates
        after_dir: the ray direction cosine following the interface in
                   interface coordinates
        dst_b4: the geometric distance from the previous interface
    and
        op_delta: optical path wrt equally inclined chords to the optical axis
    """
    eic = []

    z_dir_before = z_dir[0]
    n_before = z_dir_before*rndx[0]

    before_dir = ray[0][1]

    for i, r in enumerate(ray):
        rotT, _ = lcl_tfrms[i]
        b4_dir = rotT.dot(before_dir)

        z_dir_after = z_dir[i]
        n_after = z_dir_after*rndx[i]

        inc_pt = ray[i][0]
        eic_dst_before = ((inc_pt.dot(b4_dir) + z_dir_before*inc_pt[2]) /
                          (1.0 + z_dir_before*b4_dir[2]))

        after_dir = ray[i][1]
        eic_dst_after = ((inc_pt.dot(after_dir) + z_dir_after*inc_pt[2]) /
                         (1.0 + z_dir_after*after_dir[2]))

        dW = n_after*eic_dst_after - n_before*eic_dst_before

        eic.append([n_before, eic_dst_before, n_after, eic_dst_after, dW])
        print("e{}= {:12.5g} e{}'= {:12.5g} dW={:10.8g} n={:8.5g}"
              " n'={:8.5g}".format(i, eic_dst_before,
                                   i, eic_dst_after,
                                   dW, n_before, n_after))

        n_before = n_after
        before_dir = after_dir
        z_dir_before = z_dir_after

    P, P1k, Ps = calc_path_length(eic)
    return eic, P


def calc_path_length(eic, offset=0):
    """ given eic array, compute path length between outer surfaces
        offset is beginning index of eic array wrt the object interface
        """
    P1k = -eic[1-offset][2]*eic[1-offset][3] + eic[-2][0]*eic[-2][1]
    Ps = 0.
    for i in range(2-offset, len(eic)-2):
        Ps -= eic[i][4]
#        Ps -= eic[i][2]*eic[i][3] - eic[i][0]*eic[i][1]
    P = P1k + Ps
    return P, P1k, Ps
