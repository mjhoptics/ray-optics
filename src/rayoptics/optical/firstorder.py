#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Functions to support paraxial ray tracing a sequential optical model

.. Created on Tue Feb 13 10:48:19 2018

.. codeauthor: Michael J. Hayford
"""
import math
from collections import namedtuple
from rayoptics.optical.model_constants import Intfc, Gap, Tfrm, Indx, Zdir
from rayoptics.optical.model_constants import ht, slp, aoi
from rayoptics.optical.idealimager import ideal_imager_setup

ParaxData = namedtuple('ParaxData', ['ax_ray', 'pr_ray', 'fod'])
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
        img_dist: image distance
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

    def list_first_order_data(self):
        """ list the first order properties """
        print("efl        {:12.4g}".format(self.efl))
        print("ffl        {:12.4g}".format(self.ffl))
        print("pp1        {:12.4g}".format(self.pp1))
        print("bfl        {:12.4g}".format(self.bfl))
        print("ppk        {:12.4g}".format(self.ppk))
        print("f/#        {:12.4g}".format(self.fno))
        print("m          {:12.4g}".format(self.m))
        print("red        {:12.4g}".format(self.red))
        print("obj_dist   {:12.4g}".format(self.obj_dist))
        print("obj_ang    {:12.4g}".format(self.obj_ang))
        print("enp_dist   {:12.4g}".format(self.enp_dist))
        print("enp_radius {:12.4g}".format(self.enp_radius))
        print("na obj     {:12.4g}".format(self.obj_na))
        print("n obj      {:12.4g}".format(self.n_obj))
        print("img_dist   {:12.4g}".format(self.img_dist))
        print("img_ht     {:12.4g}".format(self.img_ht))
        print("exp_dist   {:12.4g}".format(self.exp_dist))
        print("exp_radius {:12.4g}".format(self.exp_radius))
        print("na img     {:12.4g}".format(self.img_na))
        print("n img      {:12.4g}".format(self.n_img))
        print("optical invariant {:12.4g}".format(self.opt_inv))


# paraxial_trace() - This routine performs a paraxial raytrace from object
#                    (surface 0) to image.
def paraxial_trace(path, start, start_yu, start_yu_bar):
    """ perform a paraxial raytrace of 2 linearly independent rays """
    p_ray = []
    p_ray_bar = []

    before = next(path)
    z_dir_before = before[Zdir]
    n_before = before[Indx] if z_dir_before > 0 else -before[Indx]

    b4_yui = start_yu
    b4_yui_bar = start_yu_bar
    if start == 1:
        # compute object coords from 1st surface data
        t0 = before[Gap].thi
        obj_ht = start_yu[ht] - t0*start_yu[slp]
        obj_htb = start_yu_bar[ht] - t0*start_yu_bar[slp]
        b4_yui = [obj_ht, start_yu[slp]]
        b4_yui_bar = [obj_htb, start_yu_bar[slp]]

    cv = before[Intfc].profile_cv
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
            after = next(path)

            z_dir_after = after[Zdir]
            n_after = after[Indx] if z_dir_after > 0 else -after[Indx]

            # Transfer
            t = before[Gap].thi
            cur_ht = b4_yui[ht] + t * b4_yui[slp]
            cur_htb = b4_yui_bar[ht] + t * b4_yui_bar[slp]

            # Refraction/Reflection
            srf = after[Intfc]
            k = n_before/n_after

            # calculate slope after refraction/reflection
            pwr = srf.optical_power
            cur_slp = k * b4_yui[slp] - cur_ht * pwr/n_after
            cur_slpb = k * b4_yui_bar[slp] - cur_htb * pwr/n_after

            # calculate angle of incidence (aoi)
            cv = srf.profile_cv
            aoi = cur_slp + cur_ht * cv
            aoi_bar = cur_slpb + cur_htb * cv

            yu = [cur_ht, cur_slp, aoi]
            yu_bar = [cur_htb, cur_slpb, aoi_bar]

            p_ray.append(yu)
            p_ray_bar.append(yu_bar)

            b4_yui = yu
            b4_yui_bar = yu_bar

            n_before = n_after
            z_dir_before = z_dir_after
            before = after

        except StopIteration:
            break

    return p_ray, p_ray_bar


def compute_first_order(opt_model, stop, wvl):
    """ Returns paraxial axial and chief rays, plus first order data. """
    seq_model = opt_model.seq_model
    start = 1
    n_0 = seq_model.z_dir[start-1]*seq_model.central_rndx(start-1)
    uq0 = 1/n_0
    p_ray, q_ray = paraxial_trace(seq_model.path(wl=wvl), start,
                                  [1., 0.], [0., uq0])

    n_k = seq_model.z_dir[-1]*seq_model.central_rndx(-1)
    img = -2 if seq_model.get_num_surfaces() > 2 else -1
    ak1 = p_ray[img][ht]
    bk1 = q_ray[img][ht]
    ck1 = n_k*p_ray[img][slp]
    dk1 = n_k*q_ray[img][slp]

    # print(p_ray[-2][ht], q_ray[-2][ht], n_k*p_ray[-2][slp], n_k*q_ray[-2][slp])
    # print(ak1, bk1, ck1, dk1)

    if stop is not None:
        n_s = seq_model.z_dir[stop]*seq_model.central_rndx(stop)
        as1 = p_ray[stop][ht]
        bs1 = q_ray[stop][ht]
        cs1 = n_s*p_ray[stop][slp]
        ds1 = n_s*q_ray[stop][slp]

        # find entrance pupil location w.r.t. first surface
        ybar1 = -bs1
        ubar1 = as1
        n_0 = seq_model.gaps[0].medium.rindex(wvl)
        enp_dist = -ybar1/(n_0*ubar1)

        thi0 = seq_model.gaps[0].thi

        # calculate reduction ratio for given object distance
        red = dk1 + thi0*ck1
        obj2enp_dist = thi0 + enp_dist

        yu = [0., 1.]
        pupil = opt_model.optical_spec.pupil
        aperture, obj_img_key, value_key = pupil.key
        if obj_img_key == 'object':
            if value_key == 'pupil':
                slp0 = 0.5*pupil.value/obj2enp_dist
            elif value_key == 'NA':
                slp0 = n_0*math.tan(math.asin(pupil.value/n_0))
        elif obj_img_key == 'image':
            if value_key == 'f/#':
                slpk = -1./(2.0*pupil.value)
                slp0 = slpk/red
            elif value_key == 'NA':
                slpk = n_k*math.tan(math.asin(pupil.value/n_k))
                slp0 = slpk/red
        yu = [0., slp0]

        yu_bar = [1., 0.]
        fov = opt_model.optical_spec.field_of_view
        field, obj_img_key, value_key = fov.key
        max_fld, fn = fov.max_field()
        if max_fld == 0.0:
            max_fld = 1.0
        if obj_img_key == 'object':
            if value_key == 'angle':
                ang = math.radians(max_fld)
                slpbar0 = math.tan(ang)
                ybar0 = -slpbar0*obj2enp_dist
            elif value_key == 'height':
                ybar0 = -max_fld
                slpbar0 = -ybar0/obj2enp_dist
        elif obj_img_key == 'image':
            if value_key == 'height':
                ybar0 = red*max_fld
                slpbar0 = -ybar0/obj2enp_dist
        yu_bar = [ybar0, slpbar0]

    else:  # floating stop surface - use parax_model for starting data
        ax = opt_model.parax_model.ax
        pr = opt_model.parax_model.pr
        yu = [0., ax[0][slp]/n_0]
        yu_bar = [pr[0][ht], pr[0][slp]/n_0]

    ax_ray, pr_ray = paraxial_trace(seq_model.path(wl=wvl), 0, yu, yu_bar)

    n_0 = seq_model.central_rndx(0)
    opt_inv = n_0*(ax_ray[1][ht]*pr_ray[0][slp] - pr_ray[1][ht]*ax_ray[0][slp])

    fod = FirstOrderData()
    fod.opt_inv = opt_inv
    fod.obj_dist = obj_dist = seq_model.gaps[0].thi
    fod.img_dist = img_dist = seq_model.gaps[-1].thi
    if ck1 == 0.0:
        fod.power = 0.0
        fod.efl = 0.0
        fod.pp1 = 0.0
        fod.ppk = 0.0
    else:
        fod.power = -ck1
        fod.efl = -1.0/ck1
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
    fod.obj_na = n_0*math.sin(math.atan(seq_model.z_dir[0]*ax_ray[0][slp]))
    fod.img_na = n_k*math.sin(math.atan(seq_model.z_dir[-1]*ax_ray[-1][slp]))

    return ParaxData(ax_ray, pr_ray, fod)


def compute_principle_points(path, n_0=1.0, n_k=1.0):
    """ Returns paraxial p and q rays, plus partial first order data. """
    uq0 = 1/n_0
    p_ray, q_ray = paraxial_trace(path, 0, [1., 0.], [0., uq0])

    img = -1
    ak1 = p_ray[img][ht]
    bk1 = q_ray[img][ht]
    ck1 = n_k*p_ray[img][slp]
    dk1 = n_k*q_ray[img][slp]

    # print(p_ray[-2][ht], q_ray[-2][ht], n_k*p_ray[-2][slp], n_k*q_ray[-2][slp])
    # print(ak1, bk1, ck1, dk1)

    if ck1 == 0.0:
        efl = 0.0
        pp1 = 0.0
        ppk = 0.0
    else:
        efl = -1.0/ck1
        pp1 = (dk1 - 1.0)*(n_0/ck1)
        ppk = (p_ray[-1][ht] - 1.0)*(n_k/ck1)
    ffl = pp1 - efl
    bfl = efl - ppk

    return p_ray, q_ray, (efl, pp1, ppk, ffl, bfl)


def list_parax_trace(opt_model):
    """ list the paraxial axial and chief ray data """
    seq_model = opt_model.seq_model
    ax_ray, pr_ray, fod = opt_model.optical_spec.parax_data
    print("stop surface:", seq_model.stop_surface)
    print("           y           u           n*i         ybar         ubar"
          "        n*ibar")
    for i, ax in enumerate(ax_ray):
        n = seq_model.central_rndx(i)
        n = n if seq_model.z_dir[i] > 0 else -n
        print("{:2} {:12.6g} {:12.6g} {:12.6g} {:12.6g} {:12.6g} {:12.6g}"
              .format(i, ax_ray[i][ht], ax_ray[i][slp], n*ax_ray[i][aoi],
                      pr_ray[i][ht], pr_ray[i][slp], n*pr_ray[i][aoi]))


def specsheet_from_parax_data(opt_model, specsheet):
    """ update specsheet to contents of opt_model, while preserving inputs """
    if opt_model is None:
        return None
    seq_model = opt_model.seq_model
    optical_spec = opt_model.optical_spec
    fod = optical_spec.parax_data.fod
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
