#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Functions to support paraxial ray tracing a sequential optical model

Created on Tue Feb 13 10:48:19 2018

@author: Michael J. Hayford
"""
import itertools
import math

Surf, Gap = range(2)
ht, slp = range(2)


class FirstOrderData:
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
                n_after = after[Gap].medium.rindex(wl)
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
            else:
                k = n_before/n_after

            cv = srf.profile.cv
#            print(cv, t, srf.refract_mode, n_before, n_after)

            # calculate angle of incidence (aoi)
            aoi = b4_yui[slp] + cur_ht * cv
            aoi_bar = b4_yui_bar[slp] + cur_htb * cv
            # calculate slope after refraction/reflection
            cur_slp = b4_yui[slp] + (k - 1.0)*aoi
            cur_slpb = b4_yui_bar[slp] + (k - 1.0)*aoi_bar

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


def compute_first_order(ldm, stop, wl):
    """ Returns paraxial axial and chief rays, plus first order data. """
    path = itertools.zip_longest(ldm.surfs, ldm.gaps)
    p_ray, q_ray = paraxial_trace(path, 1, [1., 0.], [0., 1.], wl)

    n_k = ldm.gaps[-1].medium.rindex(wl)
    ak1 = p_ray[-1][ht]
    bk1 = q_ray[-1][ht]
    ck1 = n_k*p_ray[-1][slp]
    dk1 = n_k*q_ray[-1][slp]

#    print(p_ray[-2][ht], q_ray[-2][ht], n_k*p_ray[-2][slp], n_k*q_ray[-2][slp])
#    print(ak1, bk1, ck1, dk1)

    n_s = ldm.gaps[stop].medium.rindex(wl)
    as1 = p_ray[stop][ht]
    bs1 = q_ray[stop][ht]
    cs1 = n_s*p_ray[stop][slp]
    ds1 = n_s*q_ray[stop][slp]

    # find entrance pupil location w.r.t. first surface
    ybar1 = -bs1
    ubar1 = as1
    n_0 = ldm.gaps[0].medium.rindex(wl)
    enp_dist = -ybar1/(n_0*ubar1)

    thi0 = ldm.gaps[0].thi

    red = dk1 + thi0*ck1
    obj2enp_dist = thi0 + enp_dist

    yu = [0., 1.]
    pupil = ldm.global_spec.pupil
    if pupil.type == 'EPD':
        slp0 = 0.5*pupil.value/obj2enp_dist
    if pupil.type == 'NAO':
        slp0 = n_0*math.tan(math.asin(pupil.value/n_0))
    if pupil.type == 'FNO':
        slpk = 0.5*pupil.value/obj2enp_dist
        slp0 = red*slpk
    if pupil.type == 'NA':
        slpk = n_k*math.tan(math.asin(pupil.value/n_k))
        slp0 = red*slpk
    yu = [0., slp0]

    yu_bar = [1., 0.]
    fov = ldm.global_spec.field_of_view
    max_fld, fn = fov.max_field()
    if fov.type == 'OBJ_ANG':
        ang = math.radians(max_fld)
        slpbar0 = math.tan(ang)
        ybar0 = -slpbar0*obj2enp_dist
    elif fov.type == 'IMG_HT':
        ybar0 = -red*max_fld
        slpbar0 = ybar0/obj2enp_dist
    else:
        ybar0 = -max_fld
        slpbar0 = ybar0/obj2enp_dist
    yu_bar = [ybar0, slpbar0]

    path = itertools.zip_longest(ldm.surfs, ldm.gaps)
    ax_ray, pr_ray = paraxial_trace(path, 0, yu, yu_bar, wl)

    n_1 = ldm.gaps[1].medium.rindex(wl)
    opt_inv = n_1*(ax_ray[1][ht]*pr_ray[1][slp] - pr_ray[1][ht]*ax_ray[1][slp])

    fod = FirstOrderData()
    fod.opt_inv = opt_inv
    fod.obj_dist = ldm.gaps[0].thi
    fod.img_dist = ldm.gaps[-1].thi
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
    fod.fno = -1.0/(2.0*ax_ray[-1][slp])
    fod.red = dk1 + ck1*fod.obj_dist
    fod.n_obj = n_0
    fod.n_img = n_k
    fod.img_ht = -fod.opt_inv/(n_k*ax_ray[-1][slp])
    fod.obj_ang = math.degrees(math.atan(pr_ray[0][slp]))
    fod.enp_dist = -pr_ray[1][ht]/(n_0*pr_ray[0][slp])
    fod.enp_radius = abs(fod.opt_inv/(n_0*pr_ray[0][slp]))
    fod.exp_dist = -(pr_ray[-1][ht]/pr_ray[-1][slp] - fod.img_dist)
    fod.exp_radius = abs(fod.opt_inv/(n_k*pr_ray[-1][slp]))

    return ax_ray, pr_ray, fod
