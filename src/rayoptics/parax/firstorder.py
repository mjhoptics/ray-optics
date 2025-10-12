#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Functions to support paraxial ray tracing a sequential optical model

.. Created on Tue Feb 13 10:48:19 2018

.. codeauthor: Michael J. Hayford
"""
import math
import numpy as np
from collections import namedtuple
import rayoptics.optical.model_constants as mc
from rayoptics.coord_geometry_types import V2d
from rayoptics.typing import Path
from rayoptics.parax.idealimager import ideal_imager_setup
from rayoptics.util import misc_math

ParaxData = namedtuple('ParaxData', ['ax_ray', 'pr_ray', 'fod'])
ParaxData.ax_ray.__doc__ = "axial marginal ray data, y, u, i"
ParaxData.pr_ray.__doc__ = "chief ray data, y, u, i"
ParaxData.fod.__doc__ = "instance of :class:`~.FirstOrderData`"


""" tuple grouping together paraxial rays and first order properties

    Attributes:
        ax_ray: axial marginal ray data, y, u, i
        pr_ray: chief ray data, y, u, i
        fod: instance of :class:`~.FirstOrderData`
"""


class FirstOrderData:
    """ Container class for first order optical properties

    All quantities are based on paraxial ray tracing. The last interface is
    the image-1 interface.

    Attributes:
        opt_inv: optical invariant
        power: optical power of system
        efl: effective focal length
        fl_obj: object space focal length, f
        fl_img: image space focal length, f'
        pp1: distance from the 1st interface to the front principle plane
        ppk: distance from the last interface to the rear principle plane
        pp_sep: distance from the front principle plane to the rear principle 
                plane
        ffl: front focal length, distance from the 1st interface to the front 
             focal point
        bfl: back focal length, distance from the last interface to the back
             focal point
        fno: focal ratio at working conjugates, f/#
        m: transverse magnification
        red: reduction ratio, -1/m
        n_obj: refractive index at central wavelength in object space
        n_img: refractive index at central wavelength in image space
        obj_dist: object distance
        img_dist: paraxial image distance
        obj_ang: maximum object angle (degrees)
        img_ht: image height
        enp_dist: entrance pupil distance from 1st interface
        enp_radius: entrance pupil radius
        exp_dist: exit pupil distance from last interface
        exp_radius: exit pupil radius
        obj_na: numerical aperture in object space
        img_na: numerical aperture in image space
    """

    def __init__(self):
        self.opt_inv: float
        self.power: float
        self.efl: float
        self.fl_obj: float
        self.fl_img: float
        self.pp1: float
        self.ppk: float
        self.pp_sep: float
        self.ffl: float
        self.bfl: float
        self.fno: float
        self.m: float
        self.red: float
        self.n_obj: float
        self.n_img: float
        self.obj_dist: float
        self.img_dist: float
        self.obj_ang: float
        self.img_ht: float
        self.enp_dist: float
        self.enp_radius: float
        self.exp_dist: float
        self.exp_radius: float
        self.obj_na: float
        self.img_na: float

    def listobj_str(self):
        o_str = f"efl        {self.efl:12.4g}\n"
        o_str += f"f          {self.fl_obj:12.4g}\n"
        o_str += f"f'         {self.fl_img:12.4g}\n"
        o_str += f"ffl        {self.ffl:12.4g}\n"
        o_str += f"pp1        {self.pp1:12.4g}\n"
        o_str += f"bfl        {self.bfl:12.4g}\n"
        o_str += f"ppk        {self.ppk:12.4g}\n"
        o_str += f"pp sep     {self.pp_sep:12.4g}\n"
        o_str += f"f/#        {self.fno:12.4g}\n"
        o_str += f"m          {self.m:12.4g}\n"
        o_str += f"red        {self.red:12.4g}\n"
        o_str += f"obj_dist   {self.obj_dist:12.4g}\n"
        o_str += f"obj_ang    {self.obj_ang:12.4g}\n"
        o_str += f"enp_dist   {self.enp_dist:12.4g}\n"
        o_str += f"enp_radius {self.enp_radius:12.4g}\n"
        o_str += f"na obj     {self.obj_na:12.4g}\n"
        o_str += f"n obj      {self.n_obj:12.4g}\n"
        o_str += f"img_dist   {self.img_dist:12.4g}\n"
        o_str += f"img_ht     {self.img_ht:12.4g}\n"
        o_str += f"exp_dist   {self.exp_dist:12.4g}\n"
        o_str += f"exp_radius {self.exp_radius:12.4g}\n"
        o_str += f"na img     {self.img_na:12.4g}\n"
        o_str += f"n img      {self.n_img:12.4g}\n"
        o_str += f"optical invariant {self.opt_inv:12.4g}"
        return o_str
    
    def list_first_order_data(self):
        """ list the first order properties """
        print(self.listobj_str())


# paraxial_trace() - This routine performs a paraxial raytrace from object
#                    (surface 0) to image.
def paraxial_trace(path, start, start_yu, start_yu_bar):
    """ perform a paraxial raytrace of 2 linearly independent rays """
    p_ray = []
    p_ray_bar = []

    b4_ifc, b4_gap, _, b4_rndx, z_dir_before = next(path)
    n_before = b4_rndx if z_dir_before > 0 else -b4_rndx

    b4_yui = start_yu
    b4_yui_bar = start_yu_bar
    if start == 1:
        # compute object coords from 1st surface data
        if np.isinf(t0 := b4_gap.thi):
            obj_ht = 0.
            obj_htb = -np.inf
        else:
            obj_ht = start_yu[mc.ht] - t0*start_yu[mc.slp]
            obj_htb = start_yu_bar[mc.ht] - t0*start_yu_bar[mc.slp]
        b4_yui = [obj_ht, start_yu[mc.slp]]
        b4_yui_bar = [obj_htb, start_yu_bar[mc.slp]]

    cv = b4_ifc.profile_cv
    # calculate angle of incidence (aoi)
    aoi = b4_yui[mc.slp] + b4_yui[mc.ht] * cv
    aoi_bar = b4_yui_bar[mc.slp] + b4_yui_bar[mc.ht] * cv
    b4_yui.append(aoi)
    b4_yui_bar.append(aoi_bar)

    p_ray.append(b4_yui)
    p_ray_bar.append(b4_yui_bar)

    # loop over remaining surfaces in path
    while True:
        try:
            ifc, gap, _, rndx, z_dir_after = next(path)
            if rndx is None:
                rndx = abs(n_before)
            if z_dir_after is None:
                z_dir_after = z_dir_before

            # Transfer
            t = b4_gap.thi
            cur_ht = b4_yui[mc.ht] + t * b4_yui[mc.slp]
            cur_htb = b4_yui_bar[mc.ht] + t * b4_yui_bar[mc.slp]

            # Refraction/Reflection
            if (ifc.interact_mode == 'dummy' or 
                ifc.interact_mode == 'phantom'):
                cur_slp = b4_yui[mc.slp]
                cur_slpb = b4_yui_bar[mc.slp]
            else:
                n_after = rndx if z_dir_after > 0 else -rndx

                k = n_before/n_after
    
                # calculate slope after refraction/reflection
                pwr = ifc.optical_power
                cur_slp = k * b4_yui[mc.slp] - cur_ht * pwr/n_after
                cur_slpb = k * b4_yui_bar[mc.slp] - cur_htb * pwr/n_after

                n_before = n_after
                z_dir_before = z_dir_after
    
            # calculate angle of incidence (aoi)
            cv = ifc.profile_cv
            aoi = cur_slp + cur_ht * cv
            aoi_bar = cur_slpb + cur_htb * cv

            yu = [cur_ht, cur_slp, aoi]
            yu_bar = [cur_htb, cur_slpb, aoi_bar]

            p_ray.append(yu)
            p_ray_bar.append(yu_bar)

            b4_yui = yu
            b4_yui_bar = yu_bar

            b4_gap = gap

        except StopIteration:
            break

    return p_ray, p_ray_bar


def paraxial_trace_ynu(path: Path, n_before: float, 
                       start_ynu: V2d, start_ynu_bar: V2d) -> tuple[list[V2d], list[V2d]]:
    """ Trace paraxial ynu & ynu_bar rays """
    ynu_ray = []
    ynu_bar_ray = []
    b4_ynu = start_ynu
    b4_ynu_bar = start_ynu_bar
    cur_ht = 1.
    cur_htb = 0.

    # loop over remaining surfaces in path
    for ifc, gap, _, rndx, z_dir_after in path:
        if rndx is None:
            rndx = abs(n_before)

        # Refraction/Reflection
        if (ifc.interact_mode == 'dummy' or 
            ifc.interact_mode == 'phantom'):
            cur_slp = b4_ynu[mc.slp]
            cur_slpb = b4_ynu_bar[mc.slp]
        else:
            if z_dir_after is not None:
                n_after = rndx if z_dir_after > 0 else -rndx

            # calculate slope after refraction/reflection
            pwr = ifc.optical_power
            cur_slp = b4_ynu[mc.slp] - cur_ht * pwr
            cur_slpb = b4_ynu_bar[mc.slp] - cur_htb * pwr

            n_before = n_after

        ynu = [cur_ht, cur_slp]
        ynu_bar = [cur_htb, cur_slpb]

        # Transfer
        if gap is not None:
            t = gap.thi
            cur_ht += t * cur_slp / n_after
            cur_htb += t * cur_slpb / n_after

        ynu_ray.append(ynu)
        ynu_bar_ray.append(ynu_bar)

        b4_ynu = ynu
        b4_ynu_bar = ynu_bar

    return ynu_ray, ynu_bar_ray


def get_parax_matrix(p_ray, q_ray, kth, n_k):
    """ Calculate transfer matrix and inverse from 1st to kth surface. """
    ak1 = p_ray[kth][mc.ht]
    bk1 = q_ray[kth][mc.ht]
    ck1 = n_k*p_ray[kth][mc.slp]
    dk1 = n_k*q_ray[kth][mc.slp]
    Mk1 = np.array([[ak1, bk1], [ck1, dk1]])
    M1k = np.array([[dk1, -bk1], [-ck1, ak1]])
    return Mk1, M1k


def compute_first_order(opt_model, stop, wvl, src_model=None):
    """ Returns paraxial axial and chief rays, plus first order data. """
    sm = opt_model['seq_model']
    osp = opt_model['optical_spec']
    start = 1
    oal = sm.overall_length()
    n_0 = sm.z_dir[start-1]*sm.central_rndx(start-1)
    n_k = sm.z_dir[-1]*sm.central_rndx(-1)
    p_ray, q_ray, pp_info = compute_principle_points(sm.path(wl=wvl), 
                                                     oal, n_0, n_k)
    img = -2 if sm.get_num_surfaces() > 2 else -1
    ak1 = p_ray[img][mc.ht]
    bk1 = q_ray[img][mc.ht]
    ck1 = n_k*p_ray[img][mc.slp]
    dk1 = n_k*q_ray[img][mc.slp]
    Mk1 = np.array([[ak1, bk1], [ck1, dk1]])
    M1k = np.array([[dk1, -bk1], [-ck1, ak1]])

    # print(p_ray[-2][mc.ht], q_ray[-2][mc.ht], n_k*p_ray[-2][mc.slp], n_k*q_ray[-2][mc.slp])
    # print(ak1, bk1, ck1, dk1)

    # The code below computes the object yu and yu_bar values
    orig_stop = stop
    if stop is None:
        # check for previously computed paraxial data and
        # use that to float the stop
        if (pm := opt_model['parax_model']) == src_model:
            if pm.pr[0][mc.slp] == 0:
                enp_dist = misc_math.infinity_guard(np.inf)
            else:
                enp_dist = -pm.pr[1][mc.ht]/pm.pr[0][mc.slp]
        elif (opt_model['analysis_results'] is not None
            and 'parax_data' in opt_model['analysis_results'] 
            and opt_model['analysis_results']['parax_data'] is not None):
                pr = opt_model['analysis_results']['parax_data'].pr_ray
                enp_dist = -pr[1][mc.ht]/(n_0*pr[0][mc.slp])
        else:  # nothing pre-computed, assume 1st surface
            stop = 1
            ybar1 = 0.0
            ubar1 = 1.0
            enp_dist = 0.0

    if stop is not None:
        n_s = sm.z_dir[stop]*sm.central_rndx(stop)
        as1 = p_ray[stop][mc.ht]
        bs1 = q_ray[stop][mc.ht]
        cs1 = n_s*p_ray[stop][mc.slp]
        ds1 = n_s*q_ray[stop][mc.slp]

        # find entrance pupil location w.r.t. first surface
        ybar1 = -bs1
        ubar1 = as1
        n_0 = sm.gaps[0].medium.rindex(wvl)
        enp_dist = -ybar1/(n_0*ubar1)
    else:
        as1 = p_ray[1][mc.ht]
        bs1 = q_ray[1][mc.ht]
        ybar1 = -bs1
        ubar1 = as1


    thi0 = sm.gaps[0].thi

    # calculate reduction ratio for given object distance
    red = dk1 + thi0*ck1
    obj2enp_dist = thi0 + enp_dist

    yu = [0., 1.]
    pupil = osp['pupil']
    aperture_spec = osp['pupil'].derive_parax_params()
    pupil_oi_key, pupil_key, pupil_value = aperture_spec
    if pupil_oi_key == 'object':
        if pupil_key == 'height':
            slp0 = pupil_value/obj2enp_dist
        elif pupil_key == 'slope':
            slp0 = pupil_value
        elif pupil_key == 'epd':
            slp0 = 0.5*pupil.value/obj2enp_dist
        elif pupil_key == 'f/#':
            slp0 = -1./(2.0*pupil.value)
        elif pupil_key == 'NA':
            slp0 = pupil.value/n_0
    elif pupil_oi_key == 'image':
        if pupil_key == 'height':
            slpk = pupil_value/obj2enp_dist
        elif pupil_key == 'slope':
            slpk = pupil_value
        elif pupil_key == 'f/#':
            slpk = -1./(2.0*pupil.value)
        elif pupil_key == 'NA':
            slpk = pupil.value/n_k
        slp0 = slpk/red
    yu = [0., slp0]

    yu_bar = [1., 0.]
    field_spec = osp['fov'].derive_parax_params()
    fov_oi_key, field_key, field_value = field_spec
    if fov_oi_key == 'object':
        if field_key == 'slope':
            slpbar0 = field_value
            ybar0 = -slpbar0*obj2enp_dist
        elif field_key == 'height':
            ybar0 = field_value
            slpbar0 = -ybar0/obj2enp_dist
    elif fov_oi_key == 'image':
        Mi1, M1i = get_parax_matrix(p_ray, q_ray, -1, n_k)
        q_ray_1 = np.array([ybar1, ubar1])
        q_ray_k = np.matmul(Mk1, q_ray_1)
        exp_dist = -q_ray_k[mc.ht]/(n_k*q_ray_k[mc.slp])
        img2exp_dist = exp_dist - sm.gaps[-1].thi

        if field_key == 'height':
            ht_i = field_value
            slp_i = -ht_i/img2exp_dist
            pr_ray_i = np.array([ht_i, slp_i])
            pr_ray_1 = np.matmul(M1i, pr_ray_i)
            slpbar0 = pr_ray_1[mc.slp]
            ybar0 = -slpbar0*obj2enp_dist
        if field_key == 'slope':
            slp_k = field_value
            ht_k = -slp_k*exp_dist
            pr_ray_k = np.array([ht_k, slp_k])
            pr_ray_1 = np.matmul(M1k, pr_ray_k)
            slpbar0 = pr_ray_1[mc.slp]
            ybar0 = -slpbar0*obj2enp_dist
    yu_bar = [ybar0, slpbar0]

    stop = orig_stop
    idx = 0

    # We have the starting coordinates, now trace the rays
    ax_ray, pr_ray = paraxial_trace(sm.path(wl=wvl), idx, yu, yu_bar)

    # Calculate the optical invariant
    opt_inv = n_0*(ax_ray[1][mc.ht]*pr_ray[0][mc.slp] - pr_ray[1][mc.ht]*ax_ray[0][mc.slp])

    # Fill in the contents of the FirstOrderData struct
    fod = FirstOrderData()
    fod.opt_inv = opt_inv
    fod.obj_dist = obj_dist = sm.gaps[0].thi
    if ck1 == 0.0:
        fod.img_dist = img_dist = 1e10
        fod.power = power = 0.0
        fod.fl_obj = fl_obj = 0.0
        fod.fl_img = fl_img = 0.0
        fod.efl = 0.0
        fod.pp1 = 0.0
        fod.ppk = 0.0
    else:
        if ax_ray[img][mc.slp] != 0:
            fod.img_dist = img_dist = -ax_ray[img][mc.ht]/ax_ray[img][mc.slp]
        else:
            fod.img_dist = img_dist = math.copysign(1e10, sm.gaps[-1].thi)
        fod.power = power = -ck1
        fod.fl_obj = fl_obj = n_0/power
        fod.fl_img = fl_img = n_k/power
        fod.efl = fl_img
        fod.pp1 = (1.0 - dk1)*(fl_obj)
        fod.ppk = (ak1 - 1.0)*(fl_img)

    fod.ffl = fod.pp1 + (-fl_obj)
    fod.bfl = fod.ppk + fl_img
    fod.pp_sep = oal - fod.pp1 + fod.ppk

    if ax_ray[img][mc.slp] != 0:
        fod.fno = -1.0/(2.0*n_k*ax_ray[-1][mc.slp])
        fod.img_ht = -fod.opt_inv/(n_k*ax_ray[-1][mc.slp])
    else:
        fod.fno = 1e10
        fod.img_ht = 1e10

    fod.m = ak1 + ck1*img_dist/n_k
    fod.red = dk1 + ck1*obj_dist
    fod.n_obj = n_0
    fod.n_img = n_k
    fod.obj_ang = math.degrees(math.atan(pr_ray[0][mc.slp]))
    if pr_ray[0][mc.slp] != 0:
        nu_pr0 = n_0*pr_ray[0][mc.slp]
        fod.enp_dist = -pr_ray[1][mc.ht]/nu_pr0
        fod.enp_radius = abs(fod.opt_inv/nu_pr0)
    else:
        fod.enp_dist = -1e10
        fod.enp_radius = 1e10

    if pr_ray[-1][mc.slp] != 0:
        fod.exp_dist = -(pr_ray[-1][mc.ht]/pr_ray[-1][mc.slp] - fod.img_dist)
        fod.exp_radius = abs(fod.opt_inv/(n_k*pr_ray[-1][mc.slp]))
    else:
        fod.exp_dist = -1e10
        fod.exp_radius = 1e10

    # compute object and image space numerical apertures
    fod.obj_na = n_0*sm.z_dir[0]*ax_ray[0][mc.slp]
    fod.img_na = n_k*sm.z_dir[-1]*ax_ray[-1][mc.slp]

    return ParaxData(ax_ray, pr_ray, fod)


def compute_principle_points(path: Path, oal: float, 
                             n_0: float=1.0, n_k: float=1.0,
                             os_idx: int=1, is_idx: int|None=None):
    """ Returns paraxial p and q rays, plus partial first order data.

    Args:
        path: an iterator containing interfaces and gaps to be traced.
              for each iteration, the sequence or generator should return a
              list containing: **Intfc, Gap, Trfm, Index, Z_Dir**
        oal: overall geometric length of the gaps in `path`
        n_0: refractive index preceding the first interface
        n_k: refractive index following last interface

    Returns:
        (p_ray, q_ray, (efl, fl_obj, fl_img, pp1, ppk, pp_sep, ffl, bfl))

        - p_ray: [ht, slp, aoi], [1, 0, -]
        - q_ray: [ht, slp, aoi], [0, 1, -]
        - power: optical power of system
        - efl: effective focal length
        - fl_obj: object space focal length, f
        - fl_img: image space focal length, f'
        - pp1: distance from the 1st interface to the front principle plane
        - ppk: distance from the last interface to the rear principle plane
        - pp_sep: distance from the front principle plane to the rear 
                  principle plane
        - ffl: front focal length, distance from the 1st interface to the 
               front focal point
        - bfl: back focal length, distance from the last interface to the back
               focal point
    """
    uq0 = 1/n_0
    p_ray, q_ray = paraxial_trace(path, os_idx, [1., 0.], [0., uq0])

    if is_idx is None:
        img = -2 if len(p_ray) > 2 else -1
    else:
        img = is_idx
    ak1 = p_ray[img][mc.ht]
    bk1 = q_ray[img][mc.ht]
    ck1 = n_k*p_ray[img][mc.slp]
    dk1 = n_k*q_ray[img][mc.slp]

    # print(p_ray[-2][mc.ht], q_ray[-2][mc.ht], n_k*p_ray[-2][mc.slp], n_k*q_ray[-2][mc.slp])
    # print(ak1, bk1, ck1, dk1)

    if ck1 == 0.0:
        power = 0.0
        fl_obj = 0.0
        fl_img = 0.0
        efl = 0.0
        pp1 = 0.0
        ppk = 0.0
    else:
        power = -ck1
        fl_obj = n_0/power
        fl_img = n_k/power
        efl = fl_img
        pp1 = (1.0 - dk1)*(fl_obj)
        ppk = (ak1 - 1.0)*(fl_img)

    ffl = pp1 + (-fl_obj)
    bfl = ppk + fl_img
    pp_sep = oal - pp1 + ppk

    return p_ray, q_ray, (power, efl, fl_obj, fl_img, 
                          pp1, ppk, pp_sep, ffl, bfl)


def compute_principle_points_for_fragment(path: Path, oal: float):
    """ Returns paraxial p and q rays, plus partial first order data.

    
    Args:
        path: an iterator containing interfaces and gaps to be traced.
              for each iteration, the sequence or generator should return a
              list containing: **Intfc, Gap, Trfm, Index, Z_Dir**
        oal: overall geometric length of the gaps in `path`

    Returns:
        (p_ray, q_ray, (efl, fl_obj, fl_img, pp1, ppk, pp_sep, ffl, bfl))

        - p_ray: [ht, n*slp], [1, 0, -]
        - q_ray: [ht, n*slp], [0, 1, -]
        - power: optical power of system
        - efl: effective focal length
        - fl_obj: object space focal length, f
        - fl_img: image space focal length, f'
        - pp1: distance from the 1st interface to the front principle plane
        - ppk: distance from the last interface to the rear principle plane
        - pp_sep: distance from the front principle plane to the rear 
                  principle plane
        - ffl: front focal length, distance from the 1st interface to the 
               front focal point
        - bfl: back focal length, distance from the last interface to the back
               focal point
    """
    n_0 = n_k = 1.
    p_ray_0 = [1., 0.]
    q_ray_0 = [0., 1.]
    p_ray, q_ray = paraxial_trace_ynu(path, n_0, p_ray_0, q_ray_0)

    img = -1
    ak1 = p_ray[img][mc.ht]
    bk1 = q_ray[img][mc.ht]
    ck1 = p_ray[img][mc.slp]
    dk1 = q_ray[img][mc.slp]

    # print(p_ray[-2][mc.ht], q_ray[-2][mc.ht], n_k*p_ray[-2][mc.slp], n_k*q_ray[-2][mc.slp])
    # print(ak1, bk1, ck1, dk1)

    if ck1 == 0.0:
        power = 0.0
        fl_obj = 0.0
        fl_img = 0.0
        efl = 0.0
        pp1 = 0.0
        ppk = 0.0
    else:
        power = -ck1
        fl_obj = n_0/power
        fl_img = n_k/power
        efl = fl_img
        pp1 = (1.0 - dk1)*(fl_obj)
        ppk = (ak1 - 1.0)*(fl_img)

    ffl = pp1 + (-fl_obj)
    bfl = ppk + fl_img
    pp_sep = oal - pp1 + ppk

    return p_ray, q_ray, (power, efl, fl_obj, fl_img, 
                          pp1, ppk, pp_sep, ffl, bfl)


def list_parax_trace_fotr(opt_model, reduced=False):
    """ list the paraxial axial and chief ray data 
    
    This function lists the paraxial data as calculated by :func:`paraxial_trace` for the seq_model.
    """
    seq_model = opt_model.seq_model
    ax_ray, pr_ray, fod = opt_model['analysis_results']['parax_data']
    list_parax_trace(ax_ray, pr_ray, seq_model, reduced)


def list_parax_trace(ax_ray, pr_ray, seq_model, reduced=False):
    """ list the ray data for the input paraxial rays. """
    num_gaps = len(seq_model.gaps)
    print("stop surface:", seq_model.stop_surface)
    print("           y           u           n*i         ybar         ubar"
          "        n*ibar")
    for i in range(len(seq_model.ifcs)):
        if i < num_gaps:
            idx = i
        else:
            idx = i - 1
        n = seq_model.central_rndx(idx)
        n = n if seq_model.z_dir[idx] > 0 else -n
        
        ax_slp = n*ax_ray[i][mc.slp] if reduced else ax_ray[i][mc.slp]
        pr_slp = n*pr_ray[i][mc.slp] if reduced else pr_ray[i][mc.slp]
        print("{:2} {:12.6g} {:12.6g} {:12.6g} {:12.6g} {:12.6g} {:12.6g}"
              .format(i, ax_ray[i][mc.ht], ax_slp, n*ax_ray[i][mc.aoi],
                      pr_ray[i][mc.ht], pr_slp, n*pr_ray[i][mc.aoi]))


def specsheet_from_parax_data(opt_model, specsheet):
    """ update specsheet to contents of opt_model, while preserving inputs """
    if opt_model is None:
        return None
    optical_spec = opt_model.optical_spec
    if opt_model['analysis_results']['parax_data'] is None:
        return specsheet
    fod = opt_model['analysis_results']['parax_data'].fod
    conj_type = optical_spec.conjugate_type('object')

    specsheet.conjugate_type = conj_type

    # specsheet.imager_inputs contains values of independent variables of
    # the optical system. Augment these as needed to get a defined imager.
    imager_inputs = dict(specsheet.imager_inputs)
    num_imager_inputs = len(imager_inputs)
    if num_imager_inputs == 0:
        # no user inputs, use model values
        if conj_type == 'finite':
            imager_inputs['m'] = fod.m
            imager_inputs['f'] = fod.efl
            specsheet.frozen_imager_inputs = [False]*5
        else:  # conj_type == 'infinite'
            imager_inputs['s'] = -math.inf
            if fod.efl != 0:
                imager_inputs['f'] = fod.efl
            specsheet.frozen_imager_inputs = [True, True, True, True, False]
    elif num_imager_inputs == 1:
        # some/partial user input specification
        if conj_type == 'finite':
            # make sure that m is specified
            if 'm' in imager_inputs:
                imager_inputs['f'] = fod.efl
            else:
                imager_inputs['m'] = fod.m
            specsheet.frozen_imager_inputs = [False]*5
        else:  # conj_type == 'infinite'
            imager_inputs['s'] = -math.inf
            if fod.efl != 0:
                imager_inputs['f'] = fod.efl
            specsheet.frozen_imager_inputs = [True, True, True, True, False]

    specsheet.imager = ideal_imager_setup(**imager_inputs)

    ape_key, ape_value = optical_spec.pupil.get_input_for_specsheet()
    fld_key, fld_value = optical_spec.field_of_view.get_input_for_specsheet()

    etendue_inputs = specsheet.etendue_inputs
    etendue_inputs['aperture'][ape_key[0]][ape_key[1]] = ape_value
    etendue_inputs['field'][fld_key[0]][fld_key[1]] = fld_value
    specsheet.generate_from_inputs(imager_inputs, etendue_inputs)

    return specsheet
