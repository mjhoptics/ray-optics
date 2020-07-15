#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Functions to support ray tracing a sequential optical model

.. Created on Thu Jan 25 11:01:04 2018

.. codeauthor: Michael J. Hayford
"""

import numpy as np
from numpy.linalg import norm
from math import sqrt, copysign

import rayoptics.optical.model_constants as mc
from rayoptics.elem.transform import (transform_before_surface,
                                      transform_after_surface)
from rayoptics.optical.model_constants import Intfc, Gap, Indx, Tfrm, Zdir
from .traceerror import (TraceMissedSurfaceError, TraceTIRError,
                         TraceEvanescentRayError)


def bend(d_in, normal, n_in, n_out):
    """ refract incoming direction, d_in, about normal """
    try:
        normal_len = norm(normal)
        cosI = np.dot(d_in, normal)/normal_len
        sinI_sqr = 1.0 - cosI*cosI
        n_cosIp = copysign(sqrt(n_out*n_out - n_in*n_in*sinI_sqr), cosI)
        alpha = n_cosIp - n_in*cosI
        d_out = (n_in*d_in + alpha*normal)/n_out
        return d_out
    except ValueError:
        raise TraceTIRError(d_in, normal, n_in, n_out)


def reflect(d_in, normal):
    """ reflect incoming direction, d_in, about normal """
    normal_len = norm(normal)
    cosI = np.dot(d_in, normal)/normal_len
    d_out = d_in - 2.0*cosI*normal
    return d_out


def phase(ifc, inc_pt, d_in, normal, wvl, n_in, n_out):
    """ apply phase shift to incoming direction, d_in, about normal """
    try:
        d_out, dW = ifc.phase(inc_pt, d_in, normal, wvl)
        return d_out, dW
    except ValueError:
        raise TraceEvanescentRayError(ifc, inc_pt, d_in, normal, n_in, n_out)


def trace(seq_model, pt0, dir0, wvl, **kwargs):
    """ fundamental raytrace function

    Args:
        seq_model: the sequential model to be traced
        pt0: starting point in coords of first interface
        dir0: starting direction cosines in coords of first interface
        wvl: wavelength in nm
        eps: accuracy tolerance for surface intersection calculation

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
    path = seq_model.path(wvl)
    kwargs['first_surf'] = kwargs.get('first_surf', 1)
    kwargs['last_surf'] = kwargs.get('last_surf',
                                     seq_model.get_num_surfaces()-2)
    return trace_raw(path, pt0, dir0, wvl, **kwargs)


def trace_raw(path, pt0, dir0, wvl, eps=1.0e-12, **kwargs):
    """ fundamental raytrace function

    Args:
        path: an iterator containing interfaces and gaps to be traced.
              for each iteration, the sequence or generator should return a
              list containing: **Intfc, Gap, Trfm, Index, Z_Dir**
        pt0: starting point in coords of first interface
        dir0: starting direction cosines in coords of first interface
        wvl: wavelength in nm
        eps: accuracy tolerance for surface intersection calculation

    Returns:
        (**ray**, **op_delta**, **wvl**)

        - **ray** is a list for each interface in **path** of these
          elements: [pt, after_dir, after_dst, normal]

            - pt: the intersection point of the ray
            - after_dir: the ray direction cosine following the interface
            - after_dst: the geometric distance to the next interface
            - normal: the surface normal at the intersection point

        - **op_delta** - optical path wrt equally inclined chords to the
          optical axis
        - **wvl** - wavelength (in nm) that the ray was traced in
    """
    ray = []
    eic = []

    print_details = kwargs.get('print_details', False)

    first_surf = kwargs.get('first_surf', 0)
    last_surf = kwargs.get('last_surf', None)

    def in_surface_range(s, include_last_surf=False):
        if first_surf == last_surf:
            return False
        if s < first_surf:
            return False
        if last_surf is None:
            return True
        else:
            return s <= last_surf if include_last_surf else s < last_surf

    # trace object surface
    obj = next(path)
    srf_obj = obj[Intfc]
    dst_b4, pt_obj = srf_obj.intersect(pt0, dir0)

    before = obj
    before_pt = pt_obj
    before_dir = dir0
    before_normal = srf_obj.normal(before_pt)
    tfrm_from_before = before[Tfrm]
    z_dir_before = before[Zdir]

    op_delta = 0.0
    opl = 0.0
    opl_eic = 0.0
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
            z_dir_after = after[Zdir]

            # intersect ray with profile
            pp_dst_intrsct, inc_pt = ifc.intersect(pp_pt_before, b4_dir,
                                                   eps=eps, z_dir=z_dir_before)
            dst_b4 = pp_dst + pp_dst_intrsct
            ray.append([before_pt, before_dir, dst_b4, before_normal])

            if in_surface_range(surf):
                opl += before[Indx] * dst_b4

            normal = ifc.normal(inc_pt)

            eic_dst_before = eic_distance_from_axis((inc_pt, b4_dir),
                                                    z_dir_before)

            # if the interface has a phase element, process that first
            if hasattr(ifc, 'phase_element'):
                doe_dir, phs = phase(ifc, inc_pt, b4_dir, normal, wvl,
                                     before[Indx], after[Indx])
                # the output of the phase element becomes the input for the
                #  refraction/reflection calculation
                b4_dir = doe_dir
                op_delta += phs

            # refract or reflect ray at interface
            if ifc.interact_mode == 'reflect':
                after_dir = reflect(b4_dir, normal)
            elif ifc.interact_mode == 'transmit':
                after_dir = bend(b4_dir, normal, before[Indx], after[Indx])
            else:  # no action, input becomes output
                after_dir = b4_dir

            eic_dst_after = eic_distance_from_axis((inc_pt, after_dir),
                                                   z_dir_after)

            surf += 1

            # Per `Hopkins, 1981 <https://dx.doi.org/10.1080/713820605>`_, the
            #  propagation direction is given by the direction cosines of the
            #  ray and therefore doesn't require the use of a negated
            #  refractive index following a reflection. Thus we use the
            #  (positive) refractive indices from the seq_model.rndx array.
            dW = after[Indx]*eic_dst_after - before[Indx]*eic_dst_before
            eic.append([before[Indx], eic_dst_before,
                        after[Indx], eic_dst_after, dW])
            if in_surface_range(surf, include_last_surf=True):
                opl_eic += dW

            if print_details:
                print("after:", surf, inc_pt, after_dir)
                print("e{}= {:12.5g} e{}'= {:12.5g} dW={:10.8g} n={:8.5g}"
                      " n'={:8.5g} zdb4={:2.0f} zdaft={:2.0f}"
                      .format(surf, eic_dst_before, surf, eic_dst_after, dW,
                              before[Indx], after[Indx],
                              z_dir_before, z_dir_after))

            before_pt = inc_pt
            before_normal = normal
            before_dir = after_dir
            z_dir_before = z_dir_after
            before = after
            tfrm_from_before = before[Tfrm]

        except TraceMissedSurfaceError as ray_miss:
            ray.append([before_pt, before_dir, pp_dst, before_normal])
            ray_miss.surf = surf+1
            ray_miss.ifc = ifc
            ray_miss.prev_tfrm = before[Tfrm]
            ray_miss.ray = ray
            raise ray_miss

        except TraceTIRError as ray_tir:
            ray.append([inc_pt, before_dir, 0.0, normal])
            ray_tir.surf = surf+1
            ray_tir.ifc = ifc
            ray_tir.int_pt = inc_pt
            ray_tir.ray = ray
            raise ray_tir

        except TraceEvanescentRayError as ray_evn:
            ray.append([inc_pt, before_dir, 0.0, normal])
            ray_evn.surf = surf+1
            ray_evn.ifc = ifc
            ray_evn.int_pt = inc_pt
            ray_evn.ray = ray
            raise ray_evn

        except StopIteration:
            ray.append([inc_pt, after_dir, 0.0, normal])
            op_delta += opl
            break

    return ray, op_delta, wvl


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
    e = ((np.dot(r[mc.p], r[mc.d]) + z_dir*r[mc.p][2]) /
         (1.0 + z_dir*r[mc.d][2]))
    return e


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


def calc_optical_path(ray, path):
    """ computes equally inclined chords and path info for ray

    Args:
        ray: ray data for traced ray
        path: an iterator containing interfaces and gaps to be traced.
              for each iteration, the sequence or generator should return a
              list containing: **Intfc, Gap, Trfm, Index, Z_Dir**

    Returns:
        **ray_op**

        - **ray_op** - the optical path between the first and last optical
          surfaces
    """
    num_items = len(ray)
    ray_seq_iter = zip(ray, path)
    next(ray_seq_iter)
    ray_op = 0
    for i in range(1, num_items-2):
        after_ray_seg, surf = next(ray_seq_iter)

        n_after = surf[Indx]
        dist = after_ray_seg[mc.dst]
        ray_op += n_after*dist

    return ray_op


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

    z_dir_before = obj_surf[Zdir]
    n_before = obj_surf[Indx]

    before_dir = before_ray_seg[mc.d]

    for i, item in enumerate(ray_seq_iter):
        after_ray_seg, surf = item

        inc_pt = after_ray_seg[mc.p]
        after_dir = after_ray_seg[mc.d]

        b4_pt, b4_dir = transform_before_surface(surf[Intfc],
                                                 (inc_pt, before_dir))
        e = eic_distance_from_axis((b4_pt, before_dir), z_dir_before)

        z_dir_after = surf[Zdir]
        aft_pt, aft_dir = transform_after_surface(surf[Intfc],
                                                  (inc_pt, after_dir))
        ep = eic_distance_from_axis((aft_pt, aft_dir), z_dir_after)

        # eic_dst_before = ((inc_pt.dot(b4_dir) + z_dir_before*inc_pt[2]) /
        #                   (1.0 + z_dir_before*b4_dir[2]))

        # Per `Hopkins, 1981 <https://dx.doi.org/10.1080/713820605>`_, the
        #  propagation direction is given by the direction cosines of the ray
        #  and therefore doesn't require the use of a negated refractive index
        #  following a reflection. Thus we use the (positive) refractive indices
        #  from the seq_model.rndx array.
        n_after = surf[Indx]
        dW = n_after*ep - n_before*e

        eic.append([n_before, e, n_after, ep, dW])
        print("e{}= {:12.5g} e{}'= {:12.5g} dW={:10.8g} n={:8.5g}"
              " n'={:8.5g}".format(i, e, i, ep, dW, n_before, n_after))

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

    for i, r in enumerate(ray):
        rotT, _ = lcl_tfrms[i]
        b4_dir = rotT.dot(before_dir)

        z_dir_after = z_dir[i]

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
        n_after = rndx[i]
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


def wave_abr(fld, wvl, foc, ray_pkg):
    """ computes optical path difference (OPD) for ray_pkg at fld and wvl

    The main references for the calculations are in the H. H. Hopkins paper
    `Calculation of the Aberrations and Image Assessment for a General Optical
    System <https://doi.org/10.1080/713820605>`_

    Args:
        fld: :class:`~.Field` point for wave aberration calculation
        wvl: wavelength of ray
        foc: :class:`~.FocusRange` instance, to specify defocus
        ray_pkg: input tuple of ray, ray_op, wvl

    Returns:
        (**opd**, **e1**, **ekp**, **ep**)

        - **opd** - OPD of ray wrt chief ray at **fld**
        - **e1** - eic in object space, prior to first interface
        - **ekp** - eic in image space, following final interface
        - **ep** - eic to reference sphere intersection

.. deprecated:: 0.4.9
    """
    return wave_abr_real_coord(fld, wvl, foc, ray_pkg)
#    return wave_abr_HHH(fld.ref_sphere_pkg, fld.chief_ray_pkg, ray_pkg)


def wave_abr_real_coord(fld, wvl, foc, ray_pkg):
    """ computes optical path difference (OPD) for ray_pkg at fld and wvl

.. deprecated:: 0.4.9
    """
    ref_sphere, parax_data, n_obj, n_img, z_dir = fld.ref_sphere
    image_pt, cr_exp_pt, cr_exp_dist, ref_dir, ref_sphere_radius = ref_sphere
    chief_ray, chief_ray_op, wvl = fld.chief_ray[0]
    ray, ray_op, wvl = ray_pkg

    k = -2  # last interface in sequence

    # eq 3.12
    e1 = eic_distance((ray[1][mc.p], ray[0][mc.d]),
                      (chief_ray[1][mc.p], chief_ray[0][mc.d]))
    # eq 3.13
    ekp = eic_distance((ray[k][mc.p], ray[k][mc.d]),
                       (chief_ray[k][mc.p], chief_ray[k][mc.d]))

    dst = ekp - cr_exp_dist

    eic_exp_pt = ray[k][mc.p] - dst*ray[k][mc.d]
    p_coord = eic_exp_pt - cr_exp_pt
    F = ref_dir.dot(ray[k][mc.d]) - ray[k][mc.d].dot(p_coord)/ref_sphere_radius
    J = p_coord.dot(p_coord)/ref_sphere_radius - 2.0*ref_dir.dot(p_coord)
    ep = J/(F + sqrt(F**2 + J/ref_sphere_radius))

    opd = -n_obj*e1 - ray_op + n_img*ekp + chief_ray_op - n_img*ep
    return opd, e1, ekp, ep


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


# def transfer_to_exit_pupil(interface, ray_seg, exp_dst_parax):
#     """Given the exiting interface and chief ray data, return exit pupil ray coords.

#     Args:
#         interface: the exiting :class:'~.Interface' for the path sequence
#         ray_seg: ray segment exiting from **interface**
#         exp_dst_parax: z distance to the paraxial exit pupil

#     Returns:
#         (**exp_pt**, **exp_dir**, **exp_dst**)

#         - **exp_pt** - ray intersection with exit pupil plane
#         - **exp_dir** - direction cosine of the ray in exit pupil space
#         - **exp_dst** - distance from interface to exit pupil pt
#     """
#     if interface.decenter:
#         # get transformation info after surf
#         r, t = interface.decenter.tform_after_surf()
#         if r is None:
#             b4_pt, b4_dir = (ray_seg[0] - t), ray_seg[1]
#         else:
#             rt = r.transpose()
#             b4_pt, b4_dir = rt.dot(ray_seg[0] - t), rt.dot(ray_seg[1])
#     else:
#         b4_pt, b4_dir = ray_seg[0], ray_seg[1]

#     h = b4_pt[0]**2 + b4_pt[1]**2
#     u = b4_dir[0]**2 + b4_dir[1]**2
#     if u == 0.0:
#         exp_dst = exp_dst_parax
#     else:
#         exp_dst = -sqrt(h/u)

#     exp_pt = b4_pt + exp_dst*b4_dir
#     exp_dir = b4_dir

#     return exp_pt, exp_dir, exp_dst
