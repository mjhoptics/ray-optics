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

from rayoptics.util.misc_math import normalize


def calculate_reference_sphere(opt_model, fld, wvl, foc,
                            chief_ray_pkg, image_pt_2d=None):
    """Compute the reference sphere for a defocussed image point at **fld**.

    Args:
        opt_model: :class:`~.OpticalModel` instance
        fld: :class:`~.Field` point for wave aberration calculation
        wvl: wavelength of ray (nm)
        foc: defocus amount
        chief_ray_pkg: input tuple of chief_ray, cr_exp_seg
        image_pt_2d: x, y image point in (defocussed) image plane, if None, use
                     the chief ray coordinate.

    Returns:
        ref_sphere: tuple of image_pt, ref_dir, ref_sphere_radius
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

    # get the image point wrt the final surface
    image_thi = opt_model['seq_model'].gaps[-1].thi
    img_pt = np.array(image_pt)
    img_pt[2] += image_thi

    # R' radius of reference sphere for O'
    ref_sphere_vec = img_pt - cr_exp_seg[mc.p]
    ref_sphere_radius = np.linalg.norm(ref_sphere_vec)
    ref_dir = normalize(ref_sphere_vec)

    ref_sphere = (image_pt, ref_dir, ref_sphere_radius)

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
    image_pt, ref_dir, ref_sphere_radius = ref_sphere
    cr, cr_exp_seg = chief_ray_pkg
    chief_ray, chief_ray_op, wvl = cr
    cr_exp_pt, cr_exp_dir, cr_exp_dist, ifc, cr_b4_pt, cr_b4_dir = cr_exp_seg

    ray, ray_op, wvl = ray_pkg

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
    opd = -n_obj*e1 - ray_op + n_img*ekp + chief_ray_op - n_img*ep

    return opd


def wave_abr_pre_calc(fod, fld, wvl, foc, ray_pkg, chief_ray_pkg):
    """Pre-calculate the part of the OPD calc independent of focus."""
    cr, cr_exp_seg = chief_ray_pkg
    chief_ray, chief_ray_op, wvl = cr
    cr_exp_pt, cr_exp_dir, cr_exp_dist, ifc, cr_b4_pt, cr_b4_dir = cr_exp_seg

    ray, ray_op, wvl = ray_pkg

    k = -2  # last interface in sequence

    # eq 3.12
    e1 = eic_distance((ray[1][mc.p], ray[0][mc.d]),
                      (chief_ray[1][mc.p], chief_ray[0][mc.d]))
    # eq 3.13
    ekp = eic_distance((ray[k][mc.p], ray[k][mc.d]),
                       (chief_ray[k][mc.p], chief_ray[k][mc.d]))

    pre_opd = -abs(fod.n_obj)*e1 - ray_op + abs(fod.n_img)*ekp + chief_ray_op

    b4_pt, b4_dir = transform_after_surface(ifc, (ray[k][mc.p], ray[k][mc.d]))
    dst = ekp - cr_exp_dist
    eic_exp_pt = b4_pt - dst*b4_dir
    p_coord = eic_exp_pt - cr_exp_pt

    return pre_opd, p_coord, b4_pt, b4_dir


def wave_abr_calc(fod, fld, wvl, foc, ray_pkg, chief_ray_pkg,
                  pre_opd_pkg, ref_sphere):
    """Given pre-calculated info and a ref. sphere, return the ray's OPD."""
    cr, cr_exp_seg = chief_ray_pkg
    image_pt, ref_dir, ref_sphere_radius = ref_sphere
    pre_opd, p_coord, b4_pt, b4_dir = pre_opd_pkg
    ray, ray_op, wvl = ray_pkg

    F = ref_dir.dot(b4_dir) - b4_dir.dot(p_coord)/ref_sphere_radius
    J = p_coord.dot(p_coord)/ref_sphere_radius - 2.0*ref_dir.dot(p_coord)

    sign_soln = -1 if ref_dir[2]*cr.ray[-1][mc.d][2] < 0 else 1
    denom = F + sign_soln*sqrt(F**2 + J/ref_sphere_radius)
    ep = 0 if denom == 0 else J/denom

    opd = pre_opd - abs(fod.n_img)*ep
    return opd
