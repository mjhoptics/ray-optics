#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Functions to support paraxial ray tracing a sequential optical model

.. Created on Tue Feb 13 10:48:19 2018

.. codeauthor: Michael J. Hayford
"""
import math
from collections import namedtuple
from rayoptics.optical.model_constants import Surf, Gap
from rayoptics.optical.model_constants import ht, slp, aoi
from rayoptics.optical.model_enums import PupilType, FieldType

ParaxData = namedtuple('ParaxData', ['ax_ray', 'pr_ray', 'fod'])


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
        self.efl = None
        self.pp1 = None
        self.ppk = None
        self.ffl = None
        self.bfl = None
        self.fno = None
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
def paraxial_trace(path, start, start_yu, start_yu_bar, wl):
    """ perform a paraxial raytrace of 2 linearly independent rays """
    p_ray = []
    p_ray_bar = []

    before = next(path)
    n_before = before[Gap].medium.rindex(wl)

    b4_yui = start_yu
    b4_yui_bar = start_yu_bar
    if start is 1:
        # compute object coords from 1st surface data
        t0 = before[Gap].thi
        obj_ht = start_yu[ht] - t0*start_yu[slp]
        obj_htb = start_yu_bar[ht] - t0*start_yu_bar[slp]
        b4_yui = [obj_ht, start_yu[slp]]
        b4_yui_bar = [obj_htb, start_yu_bar[slp]]

    cv = before[Surf].profile.cv
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
            if after[Gap]:
                n_after = math.copysign(after[Gap].medium.rindex(wl), n_before)
            else:
                n_after = n_before

            # Transfer
            t = before[Gap].thi
            cur_ht = b4_yui[ht] + t * b4_yui[slp]
            cur_htb = b4_yui_bar[ht] + t * b4_yui_bar[slp]

            # Refraction/Reflection
            srf = after[Surf]
            if srf.refract_mode == 'REFL':
                k = -1.0
                n_after = -n_after
            else:
                k = n_before/n_after

#            cv = srf.profile.cv
##            print(cv, t, srf.refract_mode, n_before, n_after)
#
#            # calculate angle of incidence (aoi)
#            aoi = b4_yui[slp] + cur_ht * cv
#            aoi_bar = b4_yui_bar[slp] + cur_htb * cv
#            # calculate slope after refraction/reflection
#            cur_slp = b4_yui[slp] + (k - 1.0)*aoi
#            cur_slpb = b4_yui_bar[slp] + (k - 1.0)*aoi_bar

            # calculate slope after refraction/reflection
            pwr = srf.optical_power
            cur_slp = k * b4_yui[slp] - cur_ht * pwr/n_after
            cur_slpb = k * b4_yui_bar[slp] - cur_htb * pwr/n_after

            # calculate angle of incidence (aoi)
            delta_n = n_after - n_before
            if delta_n == 0.0:
                aoi = cur_slp
                aoi_bar = cur_slpb
            else:
                aoi = cur_ht * pwr/delta_n + cur_slp
                aoi_bar = cur_htb * pwr/delta_n + cur_slpb

            yu = [cur_ht, cur_slp, aoi]
            yu_bar = [cur_htb, cur_slpb, aoi_bar]

            p_ray.append(yu)
            p_ray_bar.append(yu_bar)

            b4_yui = yu
            b4_yui_bar = yu_bar

            before = after
            n_before = n_after

        except StopIteration:
            break

    return p_ray, p_ray_bar


def compute_first_order(opt_model, stop, wl):
    """ Returns paraxial axial and chief rays, plus first order data. """
    seq_model = opt_model.seq_model
    p_ray, q_ray = paraxial_trace(seq_model.path(), 1, [1., 0.], [0., 1.], wl)

    n_k = seq_model.central_rndx(-1)
    ak1 = p_ray[-1][ht]
    bk1 = q_ray[-1][ht]
    ck1 = n_k*p_ray[-1][slp]
    dk1 = n_k*q_ray[-1][slp]

#    print(p_ray[-2][ht], q_ray[-2][ht], n_k*p_ray[-2][slp], n_k*q_ray[-2][slp])
#    print(ak1, bk1, ck1, dk1)

    n_s = seq_model.central_rndx(stop)
    as1 = p_ray[stop][ht]
    bs1 = q_ray[stop][ht]
    cs1 = n_s*p_ray[stop][slp]
    ds1 = n_s*q_ray[stop][slp]

    # find entrance pupil location w.r.t. first surface
    ybar1 = -bs1
    ubar1 = as1
    n_0 = seq_model.gaps[0].medium.rindex(wl)
    enp_dist = -ybar1/(n_0*ubar1)

    thi0 = seq_model.gaps[0].thi

    # calculate reduction ratio for given object distance
    red = dk1 + thi0*ck1
    obj2enp_dist = thi0 + enp_dist

    yu = [0., 1.]
    pupil = opt_model.optical_spec.pupil
    if pupil.pupil_type == PupilType.EPD:
        slp0 = 0.5*pupil.value/obj2enp_dist
    if pupil.pupil_type == PupilType.NAO:
        slp0 = n_0*math.tan(math.asin(pupil.value/n_0))
    if pupil.pupil_type == PupilType.FNO:
        slpk = -1./(2.0*pupil.value)
        slp0 = slpk/red
    if pupil.pupil_type == PupilType.NA:
        slpk = n_k*math.tan(math.asin(pupil.value/n_k))
        slp0 = slpk/red
    yu = [0., slp0]

    yu_bar = [1., 0.]
    fov = opt_model.optical_spec.field_of_view
    max_fld, fn = fov.max_field()
    if max_fld == 0.0:
        max_fld = 1.0
    if fov.field_type == FieldType.OBJ_ANG:
        ang = math.radians(max_fld)
        slpbar0 = math.tan(ang)
        ybar0 = -slpbar0*obj2enp_dist
    elif fov.field_type == FieldType.IMG_HT:
        ybar0 = red*max_fld
        slpbar0 = -ybar0/obj2enp_dist
    else:
        ybar0 = -max_fld
        slpbar0 = -ybar0/obj2enp_dist
    yu_bar = [ybar0, slpbar0]

    ax_ray, pr_ray = paraxial_trace(seq_model.path(), 0, yu, yu_bar, wl)

    n_0 = seq_model.central_rndx(0)
    opt_inv = n_0*(ax_ray[1][ht]*pr_ray[0][slp] - pr_ray[1][ht]*ax_ray[0][slp])

    fod = FirstOrderData()
    fod.opt_inv = opt_inv
    fod.obj_dist = seq_model.gaps[0].thi
    fod.img_dist = seq_model.gaps[-1].thi
    if ck1 == 0.0:
        fod.efl = 0.0
        fod.pp1 = 0.0
        fod.ppk = 0.0
    else:
        fod.efl = -1.0/ck1
        fod.pp1 = (dk1 - 1.0)*(n_0/ck1)
        fod.ppk = (p_ray[-2][ht] - 1.0)*(n_k/ck1)
    fod.ffl = fod.pp1 - fod.efl
    fod.bfl = fod.efl - fod.ppk
    fod.fno = -1.0/(2.0*n_k*ax_ray[-1][slp])
    fod.red = dk1 + ck1*fod.obj_dist
    fod.n_obj = n_0
    fod.n_img = n_k
    fod.img_ht = -fod.opt_inv/(n_k*ax_ray[-1][slp])
    fod.obj_ang = math.degrees(math.atan(pr_ray[0][slp]))
    fod.enp_dist = -pr_ray[1][ht]/(n_0*pr_ray[0][slp])
    fod.enp_radius = abs(fod.opt_inv/(n_0*pr_ray[0][slp]))
    fod.exp_dist = -(pr_ray[-1][ht]/pr_ray[-1][slp] - fod.img_dist)
    fod.exp_radius = abs(fod.opt_inv/(n_k*pr_ray[-1][slp]))
    # compute object and image space numerical apertures
    fod.obj_na = n_0*math.sin(math.atan(seq_model.z_dir[0]*ax_ray[0][slp]))
    fod.img_na = n_k*math.sin(math.atan(seq_model.z_dir[-1]*ax_ray[-1][slp]))

    return ParaxData(ax_ray, pr_ray, fod)


def list_parax_trace(opt_model):
    """ list the paraxial axial and chief ray data """
    seq_model = opt_model.seq_model
    ax_ray, pr_ray, fod = opt_model.optical_spec.parax_data
    print("stop surface:", seq_model.stop_surface)
    print("           y           u           n*i         ybar         ubar"
          "        n*ibar")
    for i, ax in enumerate(ax_ray):
        n = seq_model.central_rndx(i)
        print("{:2} {:12.6g} {:12.6g} {:12.6g} {:12.6g} {:12.6g} {:12.6g}"
              .format(i, ax_ray[i][ht], ax_ray[i][slp], n*ax_ray[i][aoi],
                      pr_ray[i][ht], pr_ray[i][slp], n*pr_ray[i][aoi]))
