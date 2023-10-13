#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Module for thin lens interface type

.. Created on Wed May 16 14:05:38 2018

.. codeauthor: Michael J. Hayford
"""


import numpy as np
from rayoptics.seq.interface import Interface
from rayoptics.oprops.doe import HolographicElement


class ThinLens(Interface):
    def __init__(self, lbl='', power=0.0, ref_index=1.5, **kwargs):
        super().__init__(interact_mode='transmit',
                         phase_element=HolographicElement(),
                         **kwargs)
        self.label = lbl
        self.optical_power = power
        self.ref_index = ref_index
        self.bending = 0.0

    def list_thinlens(self):
        if len(self.label) > 0:
            print(self.label)
        print("power: {:12.6g}".format(self.optical_power))
        self.phase_element.list_hoe()

    def __repr__(self):
        if len(self.label) > 0:
            return "{!s}(lbl={!r}, power={!r}, ref_index={!r})" \
                   .format(type(self).__name__,
                           self.label, self.optical_power, self.ref_index)
        else:
            return "{!s}(power={!r}, ref_index={!r})" \
                   .format(type(self).__name__,
                           self.optical_power, self.ref_index)

    def listobj_str(self):
        o_str = f"{self.label}: " if self.label != "" else ""
        o_str += f"{self.interact_mode}\n"
        o_str += f"thinlens: power={self.optical_power}, ref_index={self.ref_index}\n"
        o_str += f"surface_od={self.surface_od()}\n"

        return o_str

    def ifc_token(self):
        return 'l'

    def update(self):
        super().update()

    def full_profile(self, sd, flat_id=None, dir=1, steps=6):
        prf = []
        if len(sd) == 1:
            sd_lwr = -sd[0]
            sd_upr = sd[0]
        else:
            sd_lwr = sd[0]
            sd_upr = sd[1]

        prf.append([0, dir*sd_lwr])
        prf.append([0, dir*sd_upr])

        return prf

    @property
    def profile_cv(self):
        return self._power

    @profile_cv.setter
    def profile_cv(self, cv):
        self._power = cv

    def surface_od(self):
        return self.max_aperture

    def set_max_aperture(self, max_ap):
        super().set_max_aperture(max_ap)

    @property
    def optical_power(self):
        return self._power

    @optical_power.setter
    def optical_power(self, pwr):
#        print("optical_power {}: pwr={}, {} obj={}, {}".format(self.label,
#              self._power, pwr, self.phase_element.obj_pt[2], 1./pwr))
        self._power = pwr
        try:
            self.phase_element.obj_pt[2] = 1./pwr
        except ZeroDivisionError:
            self.phase_element.obj_pt[2] = 1e+10
        finally:
            self.phase_element.obj_virtual = True if pwr > 0. else False

    def set_optical_power(self, pwr, n_before, n_after):
        self.delta_n = n_after - n_before
        self.optical_power = pwr

    def apply_scale_factor(self, scale_factor):
        super().apply_scale_factor(scale_factor)
        self.optical_power = self.optical_power/scale_factor

    def update_following_reflection(self):
        super().apply_scale_factor(-1)

    def from_first_order(self, nu_before, nu_after, y):
        # nu_before used for reference point
        ref = -y/nu_before if nu_before != 0.0 else 1e+10
        obj = -y/nu_after if nu_after != 0.0 else 1e+10
#        pm = self.phase_element
#        print("from_first_order {}:\n pwr={}, {}\n ref={} ({}), {} ({})"
#              "\n obj={} ({}), {} ({})"
#              .format(self.label,
#                      self._power, (nu_before - nu_after)/y,
#                      pm.ref_pt[2], pm.ref_virtual, ref, ref > 0.,
#                      pm.obj_pt[2], pm.obj_virtual, obj, obj > 0.))
        self.phase_element.ref_pt[2] = ref
        self.phase_element.ref_virtual = True if ref > 0. else False
        self.phase_element.obj_pt[2] = obj
        self.phase_element.obj_virtual = True if obj > 0. else False
        self._power = (nu_before - nu_after)/y

    def normal(self, p):
        return np.array([0., 0., 1.])

    def intersect(self, p0, d, **kwargs):
        s1 = -p0[2]/d[2]
        p = p0 + s1*d
        return s1, p

    def phase(self, pt, d_in, normal, ifc_cntxt):
        return self.phase_element.phase(pt, d_in, normal, ifc_cntxt)
