#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 11:01:04 2018

@author: Mike
"""

import math
import numpy as np
from numpy.linalg import norm
import itertools
from . import firstorder as fo
from . import raytrace as rt


class GlobalData:
    def __init__(self):
        self.spectral_region = WvlSpec()
        self.pupil = PupilSpec()
        self.field_of_view = FieldSpec()
        self.specs = SystemSpec()
        self.parax_data = None

    def update_model(self, ldm):
        stop = ldm.stop_surface
        wl = self.spectral_region.central_wvl()
        self.parax_data = fo.compute_first_order(ldm, stop, wl)

    def trace(self, sm, pupil, fld, wl=None, eps=1.0e-12):
        if wl is None:
            wvl = self.spectral_region.central_wvl()
        else:
            wvl = self.spectral_region.wavelengths[wl]

        f = self.field_of_view.fields[fld]
        f.apply_vignetting(pupil)
        fod = self.parax_data[2]
        eprad = fod.enp_radius
        pt1 = np.array([eprad*pupil[0], eprad*pupil[1],
                        fod.obj_dist+fod.enp_dist])
        pt0 = self.obj_coords(f)
        dir0 = pt1 - pt0
        length = norm(dir0)
        dir0 = dir0/length
        path = itertools.zip_longest(sm.surfs, sm.gaps)
        return rt.trace(path, pt0, dir0, wvl, eps)

    def trace_fan(self, sm, fan_rng, fld, img_only=True, wl=None, eps=1.0e-12):
        if wl is None:
            wvl = self.spectral_region.central_wvl()
        else:
            wvl = self.spectral_region.wavelengths[wl]

        f = self.field_of_view.fields[fld]
        pt0 = self.obj_coords(f)
        fod = self.parax_data[2]
        eprad = fod.enp_radius
        start = np.array(fan_rng[0])
        stop = fan_rng[1]
        num = fan_rng[2]
        step = (stop - start)/(num - 1)
        fan = []
        for r in range(num):
            pupil = np.array(start)
            f.apply_vignetting(pupil)
            pt1 = np.array([eprad*pupil[0], eprad*pupil[1],
                            fod.obj_dist+fod.enp_dist])
            dir0 = pt1 - pt0
            length = norm(dir0)
            dir0 = dir0/length
            path = itertools.zip_longest(sm.surfs, sm.gaps)
            ray, op = rt.trace(path, pt0, dir0, wvl, eps)
            if img_only:
                fan.append([pupil, ray[-1][0]])
            else:
                fan.append([pupil, ray])

            start += step
        return fan

    def obj_coords(self, fld):
        fov = self.field_of_view
        fod = self.parax_data[2]
        if fov.type == 'OBJ_ANG':
            ang_dg = np.array([fld.x, fld.y, 0.0])
            dir_tan = np.tan(np.deg2rad(ang_dg))
            obj_pt = -dir_tan*(fod.obj_dist+fod.enp_dist)
        elif fov.type == 'IMG_HT':
            img_pt = np.array([fld.x, fld.y, 0.0])
            obj_pt = -fod.red*img_pt
        else:
            obj_pt = np.array([fld.x, fld.y, 0.0])
        return obj_pt


class WvlSpec:
    def __init__(self):
        self.wavelengths = []
        self.spectral_wts = []
        self.reference_wvl = 0
        self.coating_wvl = 550.0

    def central_wvl(self):
        return self.wavelengths[self.reference_wvl]

    def add(self, wl, wt):
        self.spectrum.append([wl, wt])
        self.spectrum.sort(key=lambda w: w[0], reverse=True)


class PupilSpec:
    types = ('EPD', 'NA', 'NAO', 'FNO')

    def __init__(self):
        self.type = 'EPD'
        self.value = 1.0


class FieldSpec:
    types = ('OBJ_ANG', 'OBJ_HT', 'IMG_HT')

    def __init__(self):
        self.fields = []
        self.type = 'OBJ_ANG'
        self.wide_angle = False

    def update_fields_cv_input(self, tla, dlist):
        if tla == 'XOB' or tla == 'YOB':
            self.type = 'OBJ_HT'
        elif tla == 'XAN' or tla == 'YAN':
            self.type = 'OBJ_ANG'
        elif tla == 'XIM' or tla == 'YIM':
            self.type = 'IMG_HT'

        if len(self.fields) == 0:
            self.fields = [Field() for f in range(len(dlist))]

        if tla[0] == 'V':
            attr = tla.lower()
        elif tla[0] == 'X' or tla[0] == 'Y':
            attr = tla[0].lower()
        elif tla == 'WTF':
            attr = 'wt'

        for i, f in enumerate(self.fields):
            f.__setattr__(attr, dlist[i])

    def max_field(self):
        max_fld = None
        max_fld_sqrd = 0.0
        for i, f in enumerate(self.fields):
            fld_sqrd = f.x*f.x + f.y*f.y
            if fld_sqrd > max_fld_sqrd:
                max_fld_sqrd = fld_sqrd
                max_fld = i
        return math.sqrt(max_fld_sqrd), max_fld


class Field:
    def __init__(self):
        self.x = 0.0
        self.y = 0.0
        self.vux = 0.0
        self.vuy = 0.0
        self.vlx = 0.0
        self.vly = 0.0
        self.wt = 1.0

    def apply_vignetting(self, pupil):
        if pupil[0] < 0.0:
            if self.vlx != 0.0:
                pupil[0] *= (1.0 - self.vlx)
        else:
            if self.vux != 0.0:
                pupil[0] *= (1.0 - self.vux)
        if pupil[1] < 0.0:
            if self.vly != 0.0:
                pupil[1] *= (1.0 - self.vly)
        else:
            if self.vuy != 0.0:
                pupil[1] *= (1.0 - self.vuy)
        return pupil


class SystemSpec:
    dims = ('M', 'C', 'I')

    def __init__(self):
        self.title = ''
        self.initials = ''
        self.dimensions = 'M'
        self.aperture_override = ''
        self.temperature = 20.0
        self.pressure = 760.0
