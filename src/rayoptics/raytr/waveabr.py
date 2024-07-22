#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2022 Michael J. Hayford
"""Wavefront aberration calculations

   Functions for setting up and calculating wavefront aberrations for 
   (fld, wvl, foc), including focus and image shift.

.. Created on Thu Mar 31 22:28:44 2022

.. codeauthor: Michael J. Hayford
"""
from math import sqrt
import numpy as np

from rayoptics.optical import model_constants as mc
from rayoptics.elem.transform import transform_after_surface

from rayoptics.util.misc_math import normalize, is_kinda_big


def calculate_reference_sphere(opt_model, fld, wvl, foc, 
                               chief_ray_pkg, 
                               image_pt_2d=None, 
                               image_delta=None):
    """Compute the reference sphere for a defocussed image point at **fld**.

        The local transform from the final interface to the image interface is 
        included to facilitate infinite refernce sphere calculations.

    Args:
        opt_model: :class:`~.OpticalModel` instance
        fld: :class:`~.Field` point for wave aberration calculation
        wvl: wavelength of ray (nm)
        foc: defocus amount
        chief_ray_pkg: input tuple of chief_ray, cr_exp_seg
        image_pt_2d: x, y image point in (defocussed) image plane, if None, use
                     the chief ray coordinate.
        image_delta: x, y displacements from image_pt_2d in (defocussed) 
                     image plane, if not None.

    Returns:
        ref_sphere: tuple of image_pt, ref_dir, ref_sphere_radius, lcl_tfrm_last
    """
    cr, cr_exp_seg = chief_ray_pkg
    # cr_exp_pt: E upper bar prime: pupil center for pencils from Q
    # cr_exp_pt, cr_b4_dir, cr_dst
    # cr_exp_pt = cr_exp_seg[mc.p]

    if image_pt_2d is None:
        # get distance along cr corresponding to a z shift of the defocus
        dist = foc / cr.ray[-1][mc.d][2]
        image_pt = cr.ray[-1][mc.p] + dist*cr.ray[-1][mc.d]
    else:
        image_pt = np.array([image_pt_2d[0], image_pt_2d[1], foc])

    if image_delta is not None:
        image_pt[:2] += image_delta

    # get the image point wrt the final surface
    seq_model = opt_model['seq_model']
    lcl_tfrm_last = seq_model.lcl_tfrms[-2]
    image_thi = seq_model.gaps[-1].thi
    img_pt = np.array(image_pt)
    img_pt[2] += image_thi

    # R' radius of reference sphere for O'
    ref_sphere_vec = img_pt - cr_exp_seg[mc.p]
    ref_sphere_radius = np.linalg.norm(ref_sphere_vec)
    ref_dir = normalize(ref_sphere_vec)

    ref_sphere = (image_pt, ref_dir, ref_sphere_radius, lcl_tfrm_last)

    return ref_sphere


def transfer_to_exit_pupil(interface, ray_seg, exp_dst_parax):
    """Given the exiting interface and chief ray data, return exit pupil ray coords.

    Args:
        interface: the exiting :class:'~.Interface' for the path sequence
        ray_seg: ray segment exiting from **interface**
        exp_dst_parax: z distance to the paraxial exit pupil

    Returns:
        (**exp_pt**, **exp_dir**, **exp_dst**)

        - **exp_pt** - ray intersection with exit pupil plane
        - **exp_dir** - direction cosine of the ray in exit pupil space
        - **exp_dst** - distance from interface to exit pupil pt
    """
    b4_pt, b4_dir = transform_after_surface(interface, ray_seg)

    # h = b4_pt[0]**2 + b4_pt[1]**2
    # u = b4_dir[0]**2 + b4_dir[1]**2
    # handle field points in the YZ plane
    h = b4_pt[1]
    u = b4_dir[1]
    if abs(u) < 1e-14:
        exp_dst = exp_dst_parax
    else:
        # exp_dst = -np.sign(b4_dir[2])*sqrt(h/u)
        exp_dst = -h/u

    exp_pt = b4_pt + exp_dst*b4_dir
    exp_dir = b4_dir

    return exp_pt, exp_dir, exp_dst, interface, b4_pt, b4_dir


# --- Wavefront aberration
def eic_distance(r, r0):
    """ calculate equally inclined chord distance between 2 rays

    Args:
        r: (p, d), where p is a point on the ray r and d is the direction
           cosine of r
        r0: (p0, d0), where p0 is a point on the ray r0 and d0 is the direction
            cosine of r0

    Returns:
        float: distance along r from equally inclined chord point to p
    """
    # eq 3.9
    e = (np.dot(r[mc.d] + r0[mc.d], r[mc.p] - r0[mc.p]) /
         (1. + np.dot(r[mc.d], r0[mc.d])))
    return e


def ray_dist_to_perp_from_pt(r, pt):
#def ray_dist_to_perp_from_pt(p, d, pt):
    """ compute distance along ray to perpendicular to `pt`. 

    Args:
        p, d: a ray, defined by point p and unit direction d
        pt: point to cast perpendicular to

    Returns:
        t: distance from p to perpendicular to pt
    """
    #return np.dot(r[mc.d], (pt - r[mc.p]))
    p, d = r
    return np.dot(d, (pt - p))


def ray_dist_to_perp_from_origin(r):
    """ compute distance along ray to perpendicular to the origin.

    Args:
        p, d: a ray, defined by point p and unit direction d

    Returns:
        t: distance from p to perpendicular to the origin
    """
    p, d = r
    return np.dot(d, -p)


def dist_to_shortest_join(r1, r2):
    """ compute distance to pts at the closest join between 2 rays. 
    
    Args:
        r1: ray 1, defined by point p1 and unit direction d1
        r2: ray 2, defined by point p2 and unit direction d2

    Returns: (p1_min, t1), (p2_min, t2)
        p1_min: point on ray 1 at the closest join
        t1: distance from p1 to p1_min
        p2_min: point on ray 2 at the closest join
        t2: distance from p2 to p2_min
    """
    p1, d1 = r1
    p2, d2 = r2
    del_p = p2 - p1
    n = np.cross(d1, d2)
    nn = np.dot(n, n)
    t1 = np.dot( np.cross(d2, n), del_p) / nn
    t2 = np.dot( np.cross(d1, n), del_p) / nn
    p1_min = p1 + t1*d1
    p2_min = p2 + t2*d2
    return (p1_min, t1), (p2_min, t2)


def wave_abr_full_calc(fod, fld, wvl, foc, ray_pkg, chief_ray_pkg, ref_sphere):
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
    image_pt, ref_dir, ref_sphere_radius, lcl_tfrm_last = ref_sphere

    if is_kinda_big(ref_sphere_radius):
        opd = wave_abr_full_calc_inf_ref(fod, fld, wvl, foc, ray_pkg,
                                         chief_ray_pkg, ref_sphere)
    else:
        opd = wave_abr_full_calc_finite_pup(fod, fld, wvl, foc, ray_pkg,
                                chief_ray_pkg, ref_sphere)

    return opd


def wave_abr_pre_calc(fod, fld, wvl, foc, ray_pkg, chief_ray_pkg, ref_sphere):
    """Pre-calculate the part of the OPD calc independent of focus."""
    image_pt, ref_dir, ref_sphere_radius, lcl_tfrm_last = ref_sphere

    if is_kinda_big(ref_sphere_radius):
        pre_opd_pkg = wave_abr_pre_calc_inf_ref(fod, fld, wvl, foc, ray_pkg,
                                                chief_ray_pkg, ref_sphere)
    else:
        pre_opd_pkg = wave_abr_pre_calc_finite_pup(fod, fld, wvl, foc, ray_pkg,
                                                   chief_ray_pkg, ref_sphere)

    return pre_opd_pkg


def wave_abr_calc(fod, fld, wvl, foc, ray_pkg, chief_ray_pkg,
                  pre_opd_pkg, ref_sphere):
    """Given pre-calculated info and a ref. sphere, return the ray's OPD."""
    image_pt, ref_dir, ref_sphere_radius, lcl_tfrm_last = ref_sphere

    if is_kinda_big(ref_sphere_radius):
        opd = wave_abr_calc_inf_ref(fod, fld, wvl, foc, ray_pkg,
                                    chief_ray_pkg, pre_opd_pkg, ref_sphere)
    else:
        opd = wave_abr_calc_finite_pup(fod, fld, wvl, foc, ray_pkg,
                                       chief_ray_pkg, pre_opd_pkg, ref_sphere)

    return opd


def wave_abr_full_calc_finite_pup(fod, fld, wvl, foc, ray_pkg, chief_ray_pkg, ref_sphere):
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
    image_pt, ref_dir, ref_sphere_radius, lcl_tfrm_last = ref_sphere
    cr, cr_exp_seg = chief_ray_pkg
    cr_ray, cr_op, wvl = cr
    cr_exp_pt, cr_exp_dir, cr_exp_dist, ifc, cr_b4_pt, cr_b4_dir = cr_exp_seg

    ray, ray_op, wvl = ray_pkg

    k = -2  # last interface in sequence

    # eq 3.12
    e1 = eic_distance((ray[1][mc.p], ray[0][mc.d]),
                      (cr_ray[1][mc.p], cr_ray[0][mc.d]))
    # eq 3.13
    ekp = eic_distance((ray[k][mc.p], ray[k][mc.d]),
                       (cr_ray[k][mc.p], cr_ray[k][mc.d]))

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
    opd = -n_obj*e1 - ray_op + n_img*ekp + cr_op - n_img*ep

    return opd


def wave_abr_pre_calc_finite_pup(fod, fld, wvl, foc, ray_pkg, chief_ray_pkg, ref_sphere):
    """Pre-calculate the part of the OPD calc independent of focus."""
    cr, cr_exp_seg = chief_ray_pkg
    cr_ray, cr_op, wvl = cr
    cr_exp_pt, cr_exp_dir, cr_exp_dist, ifc, cr_b4_pt, cr_b4_dir = cr_exp_seg

    ray, ray_op, wvl = ray_pkg

    k = -2  # last interface in sequence

    # eq 3.12
    e1 = eic_distance((ray[1][mc.p], ray[0][mc.d]),
                      (cr_ray[1][mc.p], cr_ray[0][mc.d]))
    # eq 3.13
    ekp = eic_distance((ray[k][mc.p], ray[k][mc.d]),
                       (cr_ray[k][mc.p], cr_ray[k][mc.d]))

    pre_opd = -abs(fod.n_obj)*e1 - ray_op + abs(fod.n_img)*ekp + cr_op

    b4_pt, b4_dir = transform_after_surface(ifc, (ray[k][mc.p], ray[k][mc.d]))
    dst = ekp - cr_exp_dist
    eic_exp_pt = b4_pt - dst*b4_dir
    p_coord = eic_exp_pt - cr_exp_pt

    return pre_opd, p_coord, b4_pt, b4_dir


def wave_abr_calc_finite_pup(fod, fld, wvl, foc, ray_pkg, chief_ray_pkg,
                  pre_opd_pkg, ref_sphere):
    """Given pre-calculated info and a ref. sphere, return the ray's OPD."""
    cr, cr_exp_seg = chief_ray_pkg
    image_pt, ref_dir, ref_sphere_radius, lcl_tfrm_last = ref_sphere
    pre_opd, p_coord, b4_pt, b4_dir = pre_opd_pkg
    ray, ray_op, wvl = ray_pkg

    F = ref_dir.dot(b4_dir) - b4_dir.dot(p_coord)/ref_sphere_radius
    J = p_coord.dot(p_coord)/ref_sphere_radius - 2.0*ref_dir.dot(p_coord)

    sign_soln = -1 if ref_dir[2]*cr.ray[-1][mc.d][2] < 0 else 1
    denom = F + sign_soln*sqrt(F**2 + J/ref_sphere_radius)
    ep = 0 if denom == 0 else J/denom

    opd = pre_opd - abs(fod.n_img)*ep
    return opd


def wave_abr_full_calc_inf_ref(fod, fld, wvl, foc, ray_pkg, 
                               chief_ray_pkg, ref_sphere):
    """Given a ray, a chief ray and an image pt, evaluate the OPD.

    This set of functions calculates the wavefront aberration using an
    infinite reference sphere. 
    The main references for the calculations are in the paper
    `Dependence of the wave-front aberration on the radius of the reference sphere <https://doi.org/10.1364/JOSAA.19.001187>`_ by Antonı́n Mikš.

    Args:
        fod: :class:`~.FirstOrderData` for object and image space refractive
             indices
        fld: :class:`~.Field` point for wave aberration calculation
        wvl: wavelength of ray (nm)
        foc: defocus amount
        ray_pkg: input tuple of ray, ray_op, wvl
        chief_ray_pkg: input tuple of chief_ray, cr_exp_seg
        ref_sphere: input tuple of image_pt, ref_dir, ref_sphere_radius, lcl_tfrm_last

    Returns:
        opd: OPD of ray wrt chief ray at **fld**
    """
    image_pt, ref_dir, ref_sphere_radius, lcl_tfrm_last = ref_sphere
    cr, cr_exp_seg = chief_ray_pkg
    cr_ray, cr_op, wvl = cr

    ray, ray_op, wvl = ray_pkg

    k = -2  # last interface in sequence
    n_obj = abs(fod.n_obj)
    n_img = abs(fod.n_img)

    # eq 3.12
    e1 = eic_distance((ray[1][mc.p], ray[0][mc.d]),
                      (cr_ray[1][mc.p], cr_ray[0][mc.d]))
    # eq 3.13
    ekp = eic_distance((ray[k][mc.p], ray[k][mc.d]),
                       (cr_ray[k][mc.p], cr_ray[k][mc.d]))

    if lcl_tfrm_last is not None:
        rt, t = lcl_tfrm_last   # sm.lcl_tfrms[k]

        p_b4, d_b4 = rt.dot(ray[k][mc.p] - t), rt.dot(ray[k][mc.d])
        p_cr_b4, d_cr_b4 = (rt.dot(cr_ray[k][mc.p] - t), 
                            rt.dot(cr_ray[k][mc.d]))
    else:
        p_b4, d_b4 = ray[k][mc.p], ray[k][mc.d]
        p_cr_b4, d_cr_b4 = cr_ray[k][mc.p], cr_ray[k][mc.d]

    op_b4 = ray_dist_to_perp_from_origin((p_b4, d_b4))
    op_cr_b4 = ray_dist_to_perp_from_origin((p_cr_b4, d_cr_b4))

    P1, P2 = dist_to_shortest_join((cr_ray[-1][mc.p], cr_ray[-1][mc.d]), 
                                   (ray[-1][mc.p], ray[-1][mc.d]))
    rF0 = (P1[0] + P2[0])/2

    V_B = ray_op + op_b4
    V_BE = cr_op + op_cr_b4

    W0 = V_B - V_BE + n_img * np.dot((d_b4 - d_cr_b4), rF0)

    ta = ray[-1][mc.p] - image_pt
    numer = np.dot(d_cr_b4 - d_b4*np.dot(d_b4, d_cr_b4), ta)
    denom = 1 + np.dot(d_b4, d_cr_b4)
    W_inf = W0 + n_img * numer / denom

    opd = -n_obj*e1 - W_inf

    return opd


def wave_abr_pre_calc_inf_ref(fod, fld, wvl, foc, 
                              ray_pkg, chief_ray_pkg, ref_sphere):
    """Pre-calculate the part of the OPD calc independent of focus."""
    cr, cr_exp_seg = chief_ray_pkg
    cr_ray, cr_op, wvl = cr
    cr_exp_pt, cr_exp_dir, cr_exp_dist, ifc, cr_b4_pt, cr_b4_dir = cr_exp_seg
    image_pt, ref_dir, ref_sphere_radius, lcl_tfrm_last = ref_sphere

    ray, ray_op, wvl = ray_pkg

    k = -2  # last interface in sequence
    n_obj = abs(fod.n_obj)
    n_img = abs(fod.n_img)

    # eq 3.12
    e1 = eic_distance((ray[1][mc.p], ray[0][mc.d]),
                      (cr_ray[1][mc.p], cr_ray[0][mc.d]))

    if lcl_tfrm_last is not None:
        rt, t = lcl_tfrm_last   # sm.lcl_tfrms[k]

        p_b4, d_b4 = rt.dot(ray[k][mc.p] - t), rt.dot(ray[k][mc.d])
        p_cr_b4, d_cr_b4 = (rt.dot(cr_ray[k][mc.p] - t), 
                            rt.dot(cr_ray[k][mc.d]))
    else:
        p_b4, d_b4 = ray[k][mc.p], ray[k][mc.d]
        p_cr_b4, d_cr_b4 = cr_ray[k][mc.p], cr_ray[k][mc.d]

    op_b4 = ray_dist_to_perp_from_origin((p_b4, d_b4))
    op_cr_b4 = ray_dist_to_perp_from_origin((p_cr_b4, d_cr_b4))

    P1, P2 = dist_to_shortest_join((cr_ray[-1][mc.p], cr_ray[-1][mc.d]), 
                                   (ray[-1][mc.p], ray[-1][mc.d]))
    rF0 = (P1[0] + P2[0])/2

    V_B = ray_op + op_b4
    V_BE = cr_op + op_cr_b4

    W0 = V_B - V_BE + n_img * np.dot((d_b4 - d_cr_b4), rF0)

    pre_opd = -n_obj*e1 - W0

    return pre_opd, W0, p_b4, d_b4, p_cr_b4, d_cr_b4


def wave_abr_calc_inf_ref(fod, fld, wvl, foc, ray_pkg, chief_ray_pkg,
                          pre_opd_pkg, ref_sphere):
    """Given pre-calculated info and a ref. sphere, return the ray's OPD."""
    cr, cr_exp_seg = chief_ray_pkg
    image_pt, ref_dir, ref_sphere_radius, lcl_tfrm_last = ref_sphere
    pre_opd, W0, p_b4, d_b4, p_cr_b4, d_cr_b4 = pre_opd_pkg
    ray, ray_op, wvl = ray_pkg

    n_img = abs(fod.n_img)

    ta = ray[-1][mc.p] - image_pt
    numer = np.dot(d_cr_b4 - d_b4*np.dot(d_b4, d_cr_b4), ta)
    denom = 1 + np.dot(d_b4, d_cr_b4)
    W_inf = n_img * numer / denom

    opd = pre_opd - W_inf
    return opd
