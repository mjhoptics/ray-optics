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
from rayoptics.optical.model_constants import ht, slp, aoi
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
        efl: effective focal length
        pp1: distance of front principle plane from 1st interface
        ppk: distance of rear principle plane from last interface
        ffl: front focal length
        bfl: back focal length
        fno: focal ratio at working conjugates, f/#
        red: reduction ratio
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
        self.opt_inv = None
        self.power = None
        self.efl = None
        self.pp1 = None
        self.ppk = None
        self.ffl = None
        self.bfl = None
        self.fno = None
        self.m = None
        self.red = None
        self.n_obj = None
        self.n_img = None
        self.obj_dist = None
        self.img_dist = None
        self.obj_ang = None
        self.img_ht = None
        self.enp_dist = None
        self.enp_radius = None
        self.exp_dist = None
        self.exp_radius = None
        self.obj_na = None
        self.img_na = None

    def listobj_str(self):
        """ list the first order properties """
        o_str = (f"efl        {self.efl:12.4g}\n"
                 f"ffl        {self.ffl:12.4g}\n"
                 f"pp1        {self.pp1:12.4g}\n"
                 f"bfl        {self.bfl:12.4g}\n"
                 f"ppk        {self.ppk:12.4g}\n"
                 f"m          {self.m:12.4g}\n"
                 f"red        {self.red:12.4g}\n"
                 f"obj_dist   {self.obj_dist:12.4g}\n"
                 f"obj_ang    {self.obj_ang:12.4g}\n"
                 f"enp_dist   {self.enp_dist:12.4g}\n"
                 f"enp_radius {self.enp_radius:12.4g}\n"
                 f"na obj     {self.obj_na:12.4g}\n"
                 f"n obj      {self.n_obj:12.4g}\n"
                 f"img_dist   {self.img_dist:12.4g}\n"
                 f"img_ht     {self.img_ht:12.4g}\n"
                 f"exp_dist   {self.exp_dist:12.4g}\n"
                 f"exp_radius {self.exp_radius:12.4g}\n"
                 f"f/# img    {self.fno:12.4g}\n"
                 f"na img     {self.img_na:12.4g}\n"
                 f"n img      {self.n_img:12.4g}\n"
                 f"optical invariant {self.opt_inv:12.4g}\n")
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
            obj_ht = start_yu[ht] - t0*start_yu[slp]
            obj_htb = start_yu_bar[ht] - t0*start_yu_bar[slp]
        b4_yui = [obj_ht, start_yu[slp]]
        b4_yui_bar = [obj_htb, start_yu_bar[slp]]

    cv = b4_ifc.profile_cv
    # calculate angle of incidence (aoi)
    aoi = b4_yui[slp] + b4_yui[ht] * cv
    aoi_bar = b4_yui_bar[slp] + b4_yui_bar[ht] * cv
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
            cur_ht = b4_yui[ht] + t * b4_yui[slp]
            cur_htb = b4_yui_bar[ht] + t * b4_yui_bar[slp]

            # Refraction/Reflection
            if ifc.interact_mode == 'dummy':
                cur_slp = b4_yui[slp]
                cur_slpb = b4_yui_bar[slp]
            else:
                n_after = rndx if z_dir_after > 0 else -rndx

                k = n_before/n_after
    
                # calculate slope after refraction/reflection
                pwr = ifc.optical_power
                cur_slp = k * b4_yui[slp] - cur_ht * pwr/n_after
                cur_slpb = k * b4_yui_bar[slp] - cur_htb * pwr/n_after

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


def compute_first_order(opt_model, stop, wvl, src_model=None):
    """ Returns paraxial axial and chief rays, plus first order data. """
    sm = opt_model['seq_model']
    osp = opt_model['optical_spec']
    start = 1
    n_0 = sm.z_dir[start-1]*sm.central_rndx(start-1)
    n_k = sm.z_dir[-1]*sm.central_rndx(-1)
    p_ray, q_ray, ff = compute_principle_points(sm.path(wl=wvl), 
                                                n_0, n_k)
    img = -2 if sm.get_num_surfaces() > 2 else -1
    ak1 = p_ray[img][ht]
    bk1 = q_ray[img][ht]
    ck1 = n_k*p_ray[img][slp]
    dk1 = n_k*q_ray[img][slp]

    # print(p_ray[-2][ht], q_ray[-2][ht], n_k*p_ray[-2][slp], n_k*q_ray[-2][slp])
    # print(ak1, bk1, ck1, dk1)

    # The code below computes the object yu and yu_bar values
    orig_stop = stop
    if stop is None:
        # check for previously computed paraxial data and
        # use that to float the stop
        if (pm := opt_model['parax_model']) == src_model:
            if pm.pr[0][slp] == 0:
                enp_dist = misc_math.infinity_guard(np.inf)
            else:
                enp_dist = -pm.pr[1][ht]/pm.pr[0][slp]
        elif (opt_model['analysis_results'] is not None
            and 'parax_data' in opt_model['analysis_results'] 
            and opt_model['analysis_results']['parax_data'] is not None):
                pr = opt_model['analysis_results']['parax_data'].pr_ray
                enp_dist = -pr[1][ht]/(n_0*pr[0][slp])
        else:  # nothing pre-computed, assume 1st surface
            stop = 1

    if stop is not None:
        n_s = sm.z_dir[stop]*sm.central_rndx(stop)
        as1 = p_ray[stop][ht]
        bs1 = q_ray[stop][ht]
        cs1 = n_s*p_ray[stop][slp]
        ds1 = n_s*q_ray[stop][slp]

        # find entrance pupil location w.r.t. first surface
        ybar1 = -bs1
        ubar1 = as1
        n_0 = sm.gaps[0].medium.rindex(wvl)
        enp_dist = -ybar1/(n_0*ubar1)

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
            slp0 = n_0*math.tan(math.asin(pupil.value/n_0))
    elif pupil_oi_key == 'image':
        if pupil_key == 'height':
            slpk = pupil_value/obj2enp_dist
        elif pupil_key == 'slope':
            slpk = pupil_value
        elif pupil_key == 'f/#':
            slpk = -1./(2.0*pupil.value)
        elif pupil_key == 'NA':
            slpk = n_k*math.tan(math.asin(pupil.value/n_k))
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
        if field_key == 'height':
            ybar0 = red*field_value
            slpbar0 = -ybar0/obj2enp_dist
    yu_bar = [ybar0, slpbar0]

    stop = orig_stop
    idx = 0

    # We have the starting coordinates, now trace the rays
    ax_ray, pr_ray = paraxial_trace(sm.path(wl=wvl), idx, yu, yu_bar)

    # Calculate the optical invariant
    opt_inv = n_0*(ax_ray[1][ht]*pr_ray[0][slp] - pr_ray[1][ht]*ax_ray[0][slp])

    # Fill in the contents of the FirstOrderData struct
    fod = FirstOrderData()
    fod.opt_inv = opt_inv
    fod.obj_dist = obj_dist = sm.gaps[0].thi
    if ck1 == 0.0:
        fod.img_dist = img_dist = 1e10
        fod.power = 0.0
        fod.efl = 0.0
        fod.pp1 = 0.0
        fod.ppk = 0.0
    else:
        fod.img_dist = img_dist = -ax_ray[img][ht]/ax_ray[img][slp]
        fod.power = -ck1
        fod.efl = -n_k/ck1
        fod.pp1 = (dk1 - 1.0)*(n_0/ck1)
        fod.ppk = (p_ray[-2][ht] - 1.0)*(n_k/ck1)
    fod.ffl = fod.pp1 - fod.efl
    fod.bfl = fod.efl - fod.ppk
    fod.fno = -1.0/(2.0*n_k*ax_ray[-1][slp])

    fod.m = ak1 + ck1*img_dist/n_k
    fod.red = dk1 + ck1*obj_dist
    fod.n_obj = n_0
    fod.n_img = n_k
    fod.img_ht = -fod.opt_inv/(n_k*ax_ray[-1][slp])
    fod.obj_ang = math.degrees(math.atan(pr_ray[0][slp]))
    if pr_ray[0][slp] != 0:
        nu_pr0 = n_0*pr_ray[0][slp]
        fod.enp_dist = -pr_ray[1][ht]/nu_pr0
        fod.enp_radius = abs(fod.opt_inv/nu_pr0)
    else:
        fod.enp_dist = -1e10
        fod.enp_radius = 1e10

    if pr_ray[-1][slp] != 0:
        fod.exp_dist = -(pr_ray[-1][ht]/pr_ray[-1][slp] - fod.img_dist)
        fod.exp_radius = abs(fod.opt_inv/(n_k*pr_ray[-1][slp]))
    else:
        fod.exp_dist = -1e10
        fod.exp_radius = 1e10

    # compute object and image space numerical apertures
    fod.obj_na = n_0*math.sin(math.atan(sm.z_dir[0]*ax_ray[0][slp]))
    fod.img_na = n_k*math.sin(math.atan(sm.z_dir[-1]*ax_ray[-1][slp]))

    return ParaxData(ax_ray, pr_ray, fod)


def compute_principle_points(path, n_0=1.0, n_k=1.0):
    """ Returns paraxial p and q rays, plus partial first order data.

    Args:
        path: an iterator containing interfaces and gaps to be traced.
              for each iteration, the sequence or generator should return a
              list containing: **Intfc, Gap, Trfm, Index, Z_Dir**
        n_0: refractive index preceding the first interface
        n_k: refractive index following last interface

    Returns:
        (p_ray, q_ray, (power, efl, pp1, ppk, ffl, bfl))

        - p_ray: [ht, slp, aoi], [1, 0, -]
        - q_ray: [ht, slp, aoi], [0, 1, -]
        - power: optical power
        - efl: effective focal length
        - pp1: distance of front principle plane from 1st interface
        - ppk: distance of rear principle plane from last interface
        - ffl: front focal length
        - bfl: back focal length
    """
    uq0 = 1/n_0
    p_ray, q_ray = paraxial_trace(path, 1, [1., 0.], [0., uq0])

    img = -1
    ak1 = p_ray[img][ht]
    bk1 = q_ray[img][ht]
    ck1 = n_k*p_ray[img][slp]
    dk1 = n_k*q_ray[img][slp]

    # print(p_ray[-2][ht], q_ray[-2][ht], n_k*p_ray[-2][slp], n_k*q_ray[-2][slp])
    # print(ak1, bk1, ck1, dk1)

    if ck1 == 0.0:
        power = 0.0
        efl = 0.0
        pp1 = 0.0
        ppk = 0.0
    else:
        power = -ck1
        efl = 1/power
        pp1 = (dk1 - 1.0)*(n_0/ck1)
        ppk = (ak1 - 1.0)*(n_k/ck1)
    ffl = pp1 - efl
    bfl = efl - ppk

    return p_ray, q_ray, (power, efl, pp1, ppk, ffl, bfl)


def list_parax_trace(opt_model, reduced=False):
    """ list the paraxial axial and chief ray data """
    seq_model = opt_model.seq_model
    ax_ray, pr_ray, fod = opt_model['analysis_results']['parax_data']
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
        
        ax_slp = n*ax_ray[i][slp] if reduced else ax_ray[i][slp]
        pr_slp = n*pr_ray[i][slp] if reduced else pr_ray[i][slp]
        print("{:2} {:12.6g} {:12.6g} {:12.6g} {:12.6g} {:12.6g} {:12.6g}"
              .format(i, ax_ray[i][ht], ax_slp, n*ax_ray[i][aoi],
                      pr_ray[i][ht], pr_slp, n*pr_ray[i][aoi]))


def specsheet_from_parax_data(opt_model, specsheet):
    """ update specsheet to contents of opt_model, while preserving inputs """
    if opt_model is None:
        return None
    seq_model = opt_model.seq_model
    optical_spec = opt_model.optical_spec
    if opt_model['analysis_results']['parax_data'] is None:
        return specsheet
    fod = opt_model['analysis_results']['parax_data'].fod
    conj_type = 'finite'
    if seq_model.gaps[0].thi > 10e8:
        conj_type = 'infinite'

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
    etendue_inputs[ape_key[0]][ape_key[1]][ape_key[2]] = ape_value
    etendue_inputs[fld_key[0]][fld_key[1]][fld_key[2]] = fld_value
    specsheet.generate_from_inputs(imager_inputs, etendue_inputs)

    return specsheet
