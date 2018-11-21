#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Container class for optical usage information

Created on Thu Jan 25 11:01:04 2018

@author: Michael J. Hayford
"""

import math
import numpy as np

from rayoptics.optical.firstorder import compute_first_order
import rayoptics.util.colour_system as cs
srgb = cs.cs_srgb


class OpticalSpecs:
    """ Container class for optical usage information

    Contains optical usage information to specify the aperture, field of view,
    spectrum and focal position.

    It maintains a repository of paraxial data.
    """
    def __init__(self, opt_model):
        self.opt_model = opt_model
        self.spectral_region = WvlSpec()
        self.pupil = PupilSpec()
        self.field_of_view = FieldSpec()
        self.defocus = FocusRange(0.0)
        self.parax_data = None

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['opt_model']
        del attrs['parax_data']
        return attrs

    def set_from_list(self, dl):
        self.spectral_region = dl[0]
        self.pupil = dl[1]
        self.field_of_view = dl[2]

    def sync_to_restore(self, opt_model):
        self.opt_model = opt_model

    def update_model(self):
        self.pupil.update_model()
        self.field_of_view.update_model()
        stop = self.opt_model.seq_model.stop_surface
        wvl = self.spectral_region.central_wvl()
        if not hasattr(self, 'defocus'):
            self.defocus = FocusRange(0.0)

        self.parax_data = compute_first_order(self.opt_model, stop, wvl)

    def lookup_fld_wvl_focus(self, fi, wl=None, fr=0.0):
        if wl is None:
            wvl = self.spectral_region.central_wvl()
        else:
            wvl = self.spectral_region.wavelengths[wl]
        fld = self.field_of_view.fields[fi]
        foc = self.defocus.get_focus(fr)
        return fld, wvl, foc

    def obj_coords(self, fld):
        fov = self.field_of_view
        fod = self.parax_data.fod
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
    """ Class defining a spectral region

    A spectral region is a list of wavelengths (in nm) and corresponding
    weights. A reference wavelength index defines the "center" of the
    spectral region.

    """
    def __init__(self, wlwts=[(550., 1.)], ref_wl=0):
        self.set_from_list(wlwts)
        self.reference_wvl = ref_wl
        self.coating_wvl = 550.0

    def central_wvl(self):
        return self.wavelengths[self.reference_wvl]

    def set_from_list(self, wlwts):
        self.wavelengths = []
        self.spectral_wts = []
        for wlwt in wlwts:
            self.wavelengths.append(wlwt[0])
            self.spectral_wts.append(wlwt[1])
        self.calc_colors()

    def add(self, wl, wt):
        self.wavelengths.append(wl)
        self.spectral_wts.append(wt)
        self.spectrum.sort(key=lambda w: w[0], reverse=True)

    def calc_colors(self):
        self.render_colors = []
        num_wvls = len(self.wavelengths)
        if num_wvls == 1:
            self.render_colors.append('black')
        elif num_wvls == 2:
            self.render_colors.append('red')
            self.render_colors.append('blue')
        elif num_wvls == 3:
            self.render_colors.append('red')
            self.render_colors.append('green')
            self.render_colors.append('blue')
        else:
            for w in self.wavelengths:
                print("calc_colors", w)
                rgb = srgb.wvl_to_rgb(w)
                print("rgb", rgb)
                self.render_colors.append(rgb)


class PupilSpec:
    types = ('EPD', 'NA', 'NAO', 'FNO')
    default_pupil_rays = [[0., 0.], [1., 0.], [-1., 0.], [0., 1.], [0., -1.]]
    default_ray_labels = ['00', '+X', '-X', '+Y', '-Y']

    def __init__(self, type='EPD', value=1.0):
        self.type = type
        self.value = value
        self.pupil_rays = PupilSpec.default_pupil_rays
        self.ray_labels = PupilSpec.default_ray_labels

    def set_from_list(self, ppl_spec):
        self.type = ppl_spec[0]
        self.value = ppl_spec[1]

    def update_model(self):
        if not hasattr(self, 'pupil_rays'):
            self.pupil_rays = PupilSpec.default_pupil_rays
            self.ray_labels = PupilSpec.default_ray_labels


class FieldSpec:
    """ Field of view specification

    Attributes:
        type: string with field specification type
        fields: list of Field instances
    """
    types = ('OBJ_ANG', 'OBJ_HT', 'IMG_HT')

    def __init__(self, type='OBJ_ANG', flds=[0.], wide_angle=False):
        self.type = type
        self.fields = [Field() for f in range(len(flds))]
        for i, f in enumerate(self.fields):
            f.y = flds[i]
        self.wide_angle = wide_angle

    def set_from_list(self, flds):
        self.fields = [Field() for f in range(len(flds))]
        for i, f in enumerate(self.fields):
            f.y = flds[i]

    def update_model(self):
        for f in self.fields:
            f.update()

        max_field, fi = self.max_field()
        field_norm = 1.0 if max_field == 0 else 1.0/max_field
        self.index_labels = [str(field_norm*f.y)+'F' for f in self.fields]
        self.index_labels[0] = 'axis'
        if len(self.index_labels) > 1:
            self.index_labels[-1] = 'edge'
        return self

    def update_fields_cv_input(self, tla, dlist):
        if tla == 'XOB' or tla == 'YOB':
            self.type = 'OBJ_HT'
        elif tla == 'XAN' or tla == 'YAN':
            self.type = 'OBJ_ANG'
        elif tla == 'XIM' or tla == 'YIM':
            self.type = 'IMG_HT'

        if len(self.fields) != len(dlist):
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
        """ calculates the maximum field of view

        Returns:
            magnitude of maximum field, maximum Field instance
        """
        max_fld = None
        max_fld_sqrd = 0.0
        for i, f in enumerate(self.fields):
            fld_sqrd = f.x*f.x + f.y*f.y
            if fld_sqrd > max_fld_sqrd:
                max_fld_sqrd = fld_sqrd
                max_fld = i
        return math.sqrt(max_fld_sqrd), max_fld


class Field:
    def __init__(self, x=0., y=0., wt=1.):
        self.x = x
        self.y = y
        self.vux = 0.0
        self.vuy = 0.0
        self.vlx = 0.0
        self.vly = 0.0
        self.wt = wt
        self.chief_ray = None
        self.ref_sphere = None

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['chief_ray']
        del attrs['ref_sphere']
        return attrs

    def __str__(self):
        return "{}, {}".format(self.x, self.y)

    def update(self):
        self.chief_ray = None
        self.ref_sphere = None

    def apply_vignetting(self, pupil):
        vig_pupil = pupil[:]
        if pupil[0] < 0.0:
            if self.vlx != 0.0:
                vig_pupil[0] *= (1.0 - self.vlx)
        else:
            if self.vux != 0.0:
                vig_pupil[0] *= (1.0 - self.vux)
        if pupil[1] < 0.0:
            if self.vly != 0.0:
                vig_pupil[1] *= (1.0 - self.vly)
        else:
            if self.vuy != 0.0:
                vig_pupil[1] *= (1.0 - self.vuy)
        return vig_pupil


class FocusRange:
    def __init__(self, defocus, infocus=0.0):
        self.infocus = infocus
        self.defocus = defocus

    def update(self):
        self.chief_ray = None
        self.ref_sphere = None

    def get_focus(self, fr):
        """ return focus position for input focus range parameter

        Args:
            fr: focus range parameter, -1.0 to 1.0

        Returns:
            focus position for input focus range parameter
        """
        return self.infocus + fr*self.defocus
