#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Module for thin lens interface type

Created on Wed May 16 14:05:38 2018

@author: Michael J. Hayford
"""


from math import sqrt
import numpy as np
import transforms3d as t3d
from util.misc_math import normalize
from optical.surface import Interface


class HolographicElement:
    def __init__(self, lbl=''):
        self.label = lbl
        self.ref_pt = np.array([0., 0., -1e10])
        self.ref_virtual = False
        self.obj_pt = np.array([0., 0., -1e10])
        self.obj_virtual = False
        self.ref_wl = 550.0

    def phase(self, pt, in_dir, srf_nrml, wl=None):
        normal = normalize(srf_nrml)
        ref_dir = normalize(pt - self.ref_pt)
        if self.ref_virtual:
            ref_dir = -ref_dir
        ref_cosI = np.dot(ref_dir, normal)
        obj_dir = normalize(pt - self.obj_pt)
        if self.obj_virtual:
            obj_dir = -obj_dir
        obj_cosI = np.dot(obj_dir, normal)
        in_cosI = np.dot(in_dir, normal)
        mu = 1.0 if wl is None else wl/self.ref_wl
        b = in_cosI + mu*(obj_cosI - ref_cosI)
        refp_cosI = np.dot(ref_dir, in_dir)
        objp_cosI = np.dot(obj_dir, in_dir)
        ro_cosI = np.dot(ref_dir, obj_dir)
        c = mu*(mu*(1.0 - ro_cosI) + (objp_cosI - refp_cosI))
        Q = -b + sqrt(b*b - 2*c)
        out_dir = in_dir + mu*(obj_dir - ref_dir) + Q*normal
        dW = 0.
        return out_dir, dW


class ThinLens(Interface):
    def __init__(self, lbl='', **kwargs):
        super(ThinLens, self).__init__(refract_mode='PHASE', **kwargs)
        self.label = lbl
        self.power = 0.0
        self.ref_index = 1.5
        self.bending = 0.0
        self.od = 1.0
        self.phase_mapper = HolographicElement()

    def __repr__(self):
        if len(self.label) > 0:
            return "ThinLens(%r: %r)" % (self.label, self.power)
        else:
            return "ThinLens(%r)" % (self.power)

    def update(self):
        super(ThinLens, self).update()

    def full_profile(self, sd, flat_id=None, dir=1, steps=6):
        prf = []
        prf.append([0, -dir*sd])
        prf.append([0, dir*sd])
        return prf

    def surface_od(self):
        return self.od

    def set_max_aperture(self, max_ap):
        self.od = 2.0*max_ap

    @property
    def optical_power(self):
        return self.power

    @optical_power.setter
    def optical_power(self, pwr):
        self.power = pwr
        self.phase_mapper.obj_pt[2] = 1./pwr
        self.phase_mapper.obj_virtual = True if pwr > 0. else False

    def set_optical_power(self, pwr, n_before, n_after):
        self.delta_n = n_after - n_before
        self.optical_power = pwr

    def from_first_order(self, nu_before, nu_after, y):
        # nu_before used for reference point
        ref = y/nu_before
        self.phase_mapper.ref_pt[2] = ref
        self.phase_mapper.ref_virtual = True if ref > 0. else False
        obj = y/nu_after
        self.phase_mapper.obj_pt[2] = obj
        self.phase_mapper.obj_virtual = True if obj > 0. else False
        self.optical_power = (nu_before - nu_after)/y

    def normal(self, p):
        return np.array([0., 0., 1.])

    def intersect(self, p0, d, eps):
        s1 = -p0[2]/d[2]
        p = p0 + s1*d
        return s1, p

    def phase(self, pt, d_in, normal, wl):
        return self.phase_mapper.phase(pt, d_in, normal)
