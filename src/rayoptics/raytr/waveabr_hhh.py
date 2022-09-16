#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Wavefront aberration calculations - **Experimental**

    Functions for implementing HH Hopkins canonical coordinates and path 
    difference calculations using equally inclined chords (eic). There is a
    problem somewhere in that a simple fold mirror added to the system isn't
    handled correctly by the eic calculation. See the `Newtonian with diagonal`
    lens model.

    The main references for the calculations are in the H. H. Hopkins paper
    `Calculation of the Aberrations and Image Assessment for a General Optical
    System <https://doi.org/10.1080/713820605>`_

.. Created on Tue Apr 12 22:04:59 2022

.. codeauthor: Michael J. Hayford
"""
from math import sqrt
import numpy as np

from rayoptics.optical import model_constants as mc

from . import RayPkg
from .trace import trace_base
from .waveabr import eic_distance, transfer_to_exit_pupil

from rayoptics.util.misc_math import normalize
from rayoptics.elem.transform import (transform_before_surface,
                                      transform_after_surface)


def eic_distance_from_axis(r, z_dir):
    """ calculate equally inclined chord distance between a ray and the axis

    Args:
        r: (p, d), where p is a point on the ray r and d is the direction
           cosine of r
        z_dir: direction of propagation of ray segment, +1 or -1

    Returns:
        float: distance along r from equally inclined chord point to p
    """
    # eq 3.20/3.21
    dir0 = np.array([0., 0., z_dir])
    e = np.dot(r[mc.p], r[mc.d] + dir0) / (1.0 + z_dir*r[mc.d][2])
    return e


def calc_delta_op_via_eic(ray, path):
    """ computes equally inclined chords and path info for ray

    Args:
        ray: ray data for traced ray
        path: an iterator containing interfaces and gaps to be traced.
              for each iteration, the sequence or generator should return a
              list containing: **Intfc, Gap, Trfm, Index, Z_Dir**

    Returns:
        (**eic**, **op_delta**)

        - **eic** - list of [n_before, eic_dst_before, n_after, eic_dst_after,
          dW]
        - **op_delta** - optical path wrt equally inclined chords to the
          optical axis
    """
    eic = []

    ray_seq_iter = zip(ray, path)
    before = next(ray_seq_iter)
    before_ray_seg, obj_surf = before

    z_dir_before = obj_surf[mc.Zdir]
    n_before = obj_surf[mc.Indx]

    before_dir = before_ray_seg[mc.d]

    for i, item in enumerate(ray_seq_iter):
        after_ray_seg, surf = item

        inc_pt = after_ray_seg[mc.p]
        after_dir = after_ray_seg[mc.d]
        z_dir_after = (surf[mc.Zdir] if surf[mc.Zdir] is not None
                       else z_dir_before)

        # calculate the eic values in the interface coords
        #  1) before_dir is in previous coord sys. transform it into the 
        #     current coord system.
        b4_pt, b4_dir = transform_before_surface(surf[mc.Intfc],
                                                 (inc_pt, before_dir))
        e = eic_distance_from_axis((inc_pt, b4_dir), z_dir_before)

        # 2) the after coords are already in the current coord system.
        ep = eic_distance_from_axis((inc_pt, after_dir), z_dir_after)
        
        print(" ")
        print(f"before pt{i}: {inc_pt} -> {b4_pt}")
        print(f"before dir{i}: {before_dir} -> {b4_dir}")

        # Per `Hopkins, 1981 <https://dx.doi.org/10.1080/713820605>`_, the
        #  propagation direction is given by the direction cosines of the ray
        #  and therefore doesn't require the use of a negated refractive index
        #  following a reflection. Thus we use the (positive) refractive indices
        #  from the seq_model.rndx array.
        n_after = surf[mc.Indx] if surf[mc.Indx] is not None else n_before
        dW = n_after*ep - n_before*e

        eic.append([n_before, e, n_after, ep, dW])
        print(f"e{i}= {e:12.5g} e{i}'= {ep:12.5g} dW={dW:<14.8g} "
              f"n={n_before:8.5g} n'={n_after:8.5g}")

        n_before = n_after
        before_dir = after_dir
        z_dir_before = z_dir_after

    P, P1k, Ps = calc_path_length(eic, offset=1)
    return eic, P, P1k, Ps


def eic_path_accumulation(ray, rndx, lcl_tfrms, z_dir):
    """ computes equally inclined chords and path info for ray

    Args:
        ray: ray data for traced ray
        rndx: refractive index array
        lcl_tfrms: local surface interface transformation data
        z_dir: z direction array

    Returns:
        (**eic**, **op_delta**)

        - **eic** - list of [n_before, eic_dst_before, n_after, eic_dst_after,
          dW]
        - **op_delta** - optical path wrt equally inclined chords to the
          optical axis
    """
    eic = []

    z_dir_before = z_dir[0]
    n_before = rndx[0]

    before_dir = ray[0][1]

    img_idx = len(z_dir)
    
    for i, r in enumerate(ray):
        rot, _ = lcl_tfrms[i]
        if r is None:
            b4_dir = before_dir
        else:
            rotT = rot.transpose()
            b4_dir = rotT.dot(before_dir)

        z_dir_after = z_dir[i] if i < img_idx else z_dir[-1]

        inc_pt = ray[i][0]
        eic_dst_before = eic_distance_from_axis((inc_pt, b4_dir), z_dir_before)

        after_dir = ray[i][1]
        eic_dst_after = eic_distance_from_axis((inc_pt, after_dir),
                                               z_dir_after)

        # Per `Hopkins, 1981 <https://dx.doi.org/10.1080/713820605>`_, the
        #  propagation direction is given by the direction cosines of the ray
        #  and therefore doesn't require the use of a negated refractive index
        #  following a reflection. Thus we use the (positive) refractive indices
        #  from the seq_model.rndx array.
        n_after = z_dir_after*rndx[i] if i < img_idx else z_dir_before*rndx[-1]
        dW = n_after*eic_dst_after - n_before*eic_dst_before

        eic.append([n_before, eic_dst_before, n_after, eic_dst_after, dW])
        print(f"e{i}= {eic_dst_before:12.5g} e{i}'= {eic_dst_after:12.5g} "
              f"dW={dW:10.8g} n={n_before:8.5g} n'={n_after:8.5g}")

        n_before = n_after
        before_dir = after_dir
        z_dir_before = z_dir_after

    P, P1k, Ps = calc_path_length(eic)
    return eic, P


def calc_path_length(eic, offset=0):
    """ given eic array, compute path length between outer surfaces

    Args:
        eic: equally inclined chord array
        offset (int): beginning index of eic array wrt the object interface

    Returns:
        float: path length
    """
    num_ifcs = len(eic)
    # path calc needs at least 2 interfaces, plus object and image
    if num_ifcs + offset > 3:
        # eq 3.18/3.19
        P1k = -eic[1-offset][2]*eic[1-offset][3] + eic[-2][0]*eic[-2][1]
        Ps = 0.
        for i in range(2-offset, num_ifcs-2):
            Ps -= eic[i][4]
            # Ps -= eic[i][2]*eic[i][3] - eic[i][0]*eic[i][1]
        P = P1k + Ps
        return P, P1k, Ps
    else:
        return 0., 0., 0.


# *****************
# not validated yet
# *****************
def wave_abr_full_calc_HHH(opm, fod, fld, wvl, foc, path,
                           ray_pkg, chief_ray_pkg, ref_sphere):
    """Given a ray, a chief ray and an image pt, evaluate the OPD.

    The main references for the calculations are in the H. H. Hopkins paper
    `Calculation of the Aberrations and Image Assessment for a General Optical
    System <https://doi.org/10.1080/713820605>`_

    Args:
        fod: :class:`~.FirstOrderData` for object and image space refractive
             indices
        fld: :class:`~.Field` point for wave aberration calculation
        wvl: wavelength of ray (nm)
        foc: defocus amount
        ray_pkg: input tuple of ray, ray_op, wvl
        chief_ray_pkg: input tuple of chief_ray, cr_exp_seg
        ref_sphere: input tuple of image_pt, ref_dir, ref_sphere_radius

    Returns:
        opd: OPD of ray wrt chief ray at **fld**
    """
    image_pt, ref_dir, ref_sphere_radius = ref_sphere
    cr, cr_exp_seg = chief_ray_pkg
    cr_eic, cr_P, cr_P1k, cr_Ps = calc_delta_op_via_eic(cr, opm['sm'].path())
    chief_ray, chief_ray_op, wvl = cr
    cr_exp_pt, cr_exp_dir, cr_exp_dist, ifc, cr_b4_pt, cr_b4_dir = cr_exp_seg

    ray, ray_op, wvl = ray_pkg
    eic, P, P1k, Ps = calc_delta_op_via_eic(ray, opm['sm'].path())

    k = -2  # last interface in sequence

    # eq 3.12
    e1 = eic_distance((ray[1][mc.p], ray[0][mc.d]),
                      (chief_ray[1][mc.p], chief_ray[0][mc.d]))
    # eq 3.13
    ekp = eic_distance((ray[k][mc.p], ray[k][mc.d]),
                       (chief_ray[k][mc.p], chief_ray[k][mc.d]))

    b4_pt, b4_dir = transform_after_surface(ifc, (ray[k][mc.p], ray[k][mc.d]))
    dst = ekp - cr_exp_dist
    eic_exp_pt = b4_pt - dst*b4_dir
    p_coord = eic_exp_pt - cr_exp_pt

    F = ref_dir.dot(b4_dir) - b4_dir.dot(p_coord)/ref_sphere_radius
    J = p_coord.dot(p_coord)/ref_sphere_radius - 2.0*ref_dir.dot(p_coord)

    sign_soln = -1 if ref_dir[2]*cr.ray[-1][mc.d][2] < 0 else 1
    denom = F + sign_soln*sqrt(F**2 + J/ref_sphere_radius)
    ep = 0 if denom == 0 else J/denom

    n_obj = abs(fod.n_obj)
    n_img = abs(fod.n_img)
    opd = -n_obj*e1 - P + n_img*ekp + cr_P - n_img*ep

    return opd


# HH Hopkins Canonical Coordinates - experimental
def setup_canonical_coords(opt_model, fld, wvl, image_pt=None):
    seq_model = opt_model.seq_model
    parax_data = opt_model['analysis_results']['parax_data']
    fod = parax_data.fod

    if fld.chief_ray is None:
        ray, op, wvl = trace_base(opt_model, [0., 0.], fld, wvl)
        fld.chief_ray = RayPkg(ray, op, wvl)
    cr = fld.chief_ray

    if image_pt is None:
        image_pt = cr.ray[-1][mc.p]

    # cr_exp_pt: E upper bar prime: pupil center for pencils from Q
    # cr_exp_pt, cr_b4_dir, cr_dst
    cr_exp_seg = transfer_to_exit_pupil(seq_model.ifcs[-2],
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
    ref_sphere_pkg = (ref_sphere, parax_data, n_obj, n_img, z_dir)
    fld.ref_sphere = ref_sphere_pkg
    return ref_sphere_pkg, cr


def wave_abr_HHH(fld, wvl, foc, ray_pkg):
    """ computes optical path difference (OPD) for ray_pkg at fld and wvl

.. deprecated:: 0.4.9
    """
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
    e1 = eic_distance((ray[1][mc.p], ray[0][mc.d]),
                      (chief_ray[1][mc.p], chief_ray[0][mc.d]))
    # eq 3.13
    ekp = eic_distance((ray[k][mc.p], ray[k][mc.d]),
                       (chief_ray[k][mc.p], chief_ray[k][mc.d]))

    # eq 4.33
    eic_pt = ray[k][mc.p] - ekp*ray[k][mc.d]
    print("eic_pt", eic_pt)

    Nk_cr = chief_ray[k][mc.d][2]
    Zk_cr = chief_ray[k][mc.p][2]

    def reduced_pupil_coord(X, L, e):
        coef1 = -n_img/(H * Nk_cr)
        xp = coef1*(pr_k[slp]*(Nk_cr*(X - L*e) - L*Zk_cr) +
                    pr_k[ht]*L)
        return xp

    # eq 5.4
    xp_ray = reduced_pupil_coord(ray[k][mc.p][0], ray[k][mc.d][0], ekp)
    yp_ray = reduced_pupil_coord(ray[k][mc.p][1], ray[k][mc.d][1], ekp)
    # eq 5.5
    xp_cr = reduced_pupil_coord(chief_ray[k][mc.p][0], chief_ray[k][mc.d][0], 0.)
    yp_cr = reduced_pupil_coord(chief_ray[k][mc.p][1], chief_ray[k][mc.d][1], 0.)
    # eq 5.6
    zp_ray = -(((ray[k][mc.d][0] + chief_ray[k][mc.d][0])*(xp_ray - xp_cr) +
                (ray[k][mc.d][1] + chief_ray[k][mc.d][1])*(yp_ray - yp_cr)) /
                (ray[k][mc.d][2] + chief_ray[k][mc.d][2]))

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
    G0_ray = reduced_image_coord(ray[k][mc.p][0], ray[k][mc.d][0],
                                 ray[k][mc.p][2], ray[k][mc.d][2])
    H0_ray = reduced_image_coord(ray[k][mc.p][1], ray[k][mc.d][1],
                                 ray[k][mc.p][2], ray[k][mc.d][2])
    # eq 5.14
    G0_cr = reduced_image_coord(chief_ray[k][mc.p][0], chief_ray[k][mc.d][0],
                                Zk_cr, Nk_cr)
    H0_cr = reduced_image_coord(chief_ray[k][mc.p][1], chief_ray[k][mc.d][1],
                                Zk_cr, Nk_cr)
    print("G0, H0_ref; G0, H0_cr:", G0_ref, H0_ref, G0_cr, H0_cr)
    # eq 4.17
    a = pr_k[slp]*G0_ref/H + ax_k[slp]*xp_cr
    b = pr_k[slp]*H0_ref/H + ax_k[slp]*yp_cr
    g = z_dir/sqrt(1. - a**2 + b**2)
    # eq 4.18
    ref_dir = np.array([-a*g, -b*g, g])
    print("ref_dir, cr_dir", ref_dir, chief_ray[k][mc.d])
    # eq 4.25
    F = (np.dot(ref_dir, ray[k][mc.d]) + ref_dir[2]*ax_k[slp] *
         np.dot(chief_ray[k][mc.d], (rpc_ray - rpc_cr)))

    # eq 4.28
    Ja = (ref_dir[2]*ax_k[slp] *
          np.dot((eic_pt - chief_ray[k][0]), (rpc_ray - rpc_cr)))
    Jb = -(2.0*(ref_dir[2]/Nk_cr)*(ax_k[ht] - ax_k[slp]*Zk_cr) *
           np.dot(chief_ray[k][mc.d], (rpc_ray - rpc_cr)))
    Jc = (2.0*(ref_dir[2]/n_img)*((G0_cr - G0_ref)*(xp_ray - xp_cr) +
                                  (H0_cr - H0_ref)*(yp_ray - yp_cr)))
    J = Ja + Jb + Jc
    print("F, J, Ja, Jb, Jc", F, J, Ja, Jb, Jc)
#    J = ((ref_dir[2]*ax_k[slp] *
#         np.dot((eic_pt - chief_ray[k][0]), (rpc_ray - rpc_cr))) -
#         2.0*(ref_dir[2]/Nk_cr)*(ax_k[ht] - ax_k[slp]*Zk_cr) *
#         np.dot(chief_ray[k][mc.d], (rpc_ray - rpc_cr)) +
#         2.0*(ref_dir[2]/n_img)*((G0_cr - G0_ref)*(xp_ray - xp_cr) +
#                                 (H0_cr - H0_ref)*(yp_ray - yp_cr)))

#    # eq 4.29 Q' = image_pt
#    F = (np.dot(chief_ray[k][mc.d], ray[k][mc.d]) + Nk_cr*ax_k[slp] *
#         np.dot(chief_ray[k][mc.d], (rpc_ray - rpc_cr)))
#    J = ((Nk_cr*ax_k[slp] *
#         np.dot((eic_pt - chief_ray[k][0]), (rpc_ray - rpc_cr))) -
#         2.0*(ax_k[ht] - ax_k[slp]*Zk_cr) *
#         np.dot(chief_ray[k][mc.d], (rpc_ray - rpc_cr)))

    # eq 4.21
    ep = J/(F + sqrt(F**2 + J*(n_img*ax_k[slp]*pr_k[slp]*ref_dir[2] / H)))

    print("F, J, ep (canon)", F, J, ep)

    # eq 3.14/3.22, using 3.15
    opd = -n_obj*e1 - ray_op + n_img*ekp + chief_ray_op - n_img*ep
    return opd, e1, ekp, ep
