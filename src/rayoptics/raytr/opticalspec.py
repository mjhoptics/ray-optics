#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Container class for optical usage information

.. Created on Thu Jan 25 11:01:04 2018

.. codeauthor: Michael J. Hayford
"""

import math
import numpy as np

from rayoptics.parax.firstorder import compute_first_order, list_parax_trace
from rayoptics.raytr.trace import aim_chief_ray
from rayoptics.optical import model_enums
import rayoptics.optical.model_constants as mc
from opticalglass.spectral_lines import get_wavelength
import rayoptics.util.colour_system as cs
from rayoptics.util import colors
srgb = cs.cs_srgb


class OpticalSpecs:
    """ Container class for optical usage information

    Contains optical usage information to specify the aperture, field of view,
    spectrum and focal position. These can be accessed via the mapping
    interface:

        - self['wvls']: instance of :class:`~.WvlSpec`
        - self['pupil']: instance of :class:`~.PupilSpec`
        - self['fov']: instance of :class:`~.FieldSpec`
        - self['focus']: instance of :class:`~.FocusRange`

    It also maintains a repository of paraxial data.

    Attributes:
        parax_data: tuple of :obj:`~.firstorder.ParaxData`

    """

    do_aiming_default = True

    def __init__(self, opt_model, specsheet=None, **kwargs):
        self.opt_model = opt_model
        self._submodels = {}
        self['wvls'] = WvlSpec(**kwargs)
        self['pupil'] = PupilSpec(self)
        self['fov'] = FieldSpec(self)
        self['focus'] = FocusRange(0.0)
        self.parax_data = None
        self.do_aiming = OpticalSpecs.do_aiming_default
        if specsheet:
            self.set_from_specsheet(specsheet)
        
    def __getitem__(self, key):
        """ Provide mapping interface to submodels. """
        return self._submodels[key]
        
    def __setitem__(self, key, value):
        """ Provide mapping interface to submodels. """
        self._submodels[key] = value

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['opt_model']
        del attrs['_submodels']
        del attrs['parax_data']
        del attrs['do_aiming']

        attrs['spectral_region'] = self['wvls']
        attrs['pupil'] = self['pupil']
        attrs['field_of_view'] = self['fov']
        attrs['defocus'] = self['focus']

        return attrs

    def __json_decode__(self, **attrs):
        submodels = {}
        submodels['wvls'] = attrs['spectral_region']
        submodels['pupil'] = attrs['pupil']
        submodels['fov'] = attrs['field_of_view']
        submodels['focus'] = (attrs['defocus'] if 'defocus' in attrs
                              else FocusRange(0.0))

        self._submodels = submodels

    @property
    def spectral_region(self):
        return self._submodels['wvls']

    @spectral_region.setter
    def spectral_region(self, sr):
        self._submodels['wvls'] = sr

    @property
    def pupil(self):
        return self._submodels['pupil']

    @pupil.setter
    def pupil(self, pup):
        self._submodels['pupil'] = pup

    @property
    def field_of_view(self):
        return self._submodels['fov']

    @field_of_view.setter
    def field_of_view(self, fov):
        self._submodels['fov'] = fov

    @property
    def defocus(self):
        return self._submodels['focus']

    @defocus.setter
    def defocus(self, foc):
        self._submodels['focus'] = foc
    
    def set_from_list(self, dl):
        self.spectral_region = dl[0]
        self.pupil = dl[1]
        self.field_of_view = dl[2]

    def set_from_specsheet(self, ss):
        self.spectral_region.set_from_specsheet(ss)
        self.pupil.set_from_specsheet(ss)
        self.field_of_view.set_from_specsheet(ss)
        self.defocus.set_from_specsheet(ss)

    def sync_to_restore(self, opt_model):
        self.opt_model = opt_model
        if not hasattr(self, 'do_aiming'):
            self.do_aiming = OpticalSpecs.do_aiming_default

        self['wvls'].sync_to_restore(self)
        self['pupil'].sync_to_restore(self)
        self['fov'].sync_to_restore(self)

    def update_model(self, **kwargs):
        self.spectral_region.update_model(**kwargs)
        self.pupil.update_model(**kwargs)
        self.field_of_view.update_model(**kwargs)
        stop = self.opt_model.seq_model.stop_surface
        wvl = self.spectral_region.central_wvl

        if self.opt_model.seq_model.get_num_surfaces() > 2:
            self.parax_data = compute_first_order(self.opt_model, stop, wvl)
            if self.do_aiming:
                for i, fld in enumerate(self.field_of_view.fields):
                    aim_pt = aim_chief_ray(self.opt_model, fld, wvl)
                    fld.aim_pt = aim_pt

    def lookup_fld_wvl_focus(self, fi, wl=None, fr=0.0):
        """ returns field, wavelength and defocus data

        Args:
            fi (int): index into the field_of_view list of Fields
            wl (int): index into the spectral_region list of wavelengths
            fr (float): focus range parameter, -1.0 to 1.0

        Returns:
            (**fld**, **wvl**, **foc**)

            - **fld** - :class:`Field` instance for field_of_view[fi]
            - **wvl** - wavelength in nm
            - **foc** - focus shift from image interface
        """
        if wl is None:
            wvl = self.spectral_region.central_wvl
        else:
            wvl = self.spectral_region.wavelengths[wl]
        fld = self.field_of_view.fields[fi]
        foc = self.defocus.get_focus(fr)
        return fld, wvl, foc

    def obj_coords(self, fld):
        return self.field_of_view.obj_coords(fld)

    def list_first_order_data(self):
        self.parax_data.fod.list_first_order_data()

    def list_parax_trace(self, **kwargs):
        list_parax_trace(self.opt_model, **kwargs)


class WvlSpec:
    """ Class defining a spectral region

    A spectral region is a list of wavelengths (in nm) and corresponding
    weights. The central wavelength of the spectral region is central_wvl.
    The index into the wavelength list for central_wvl is reference_wvl.

    """

    def __init__(self, wlwts=[('d', 1.)], ref_wl=0, do_init=True, **kwargs):
        if do_init:
            self.set_from_list(wlwts)
        else:
            self.wavelengths = []
            self.spectral_wts = []
        self.reference_wvl = ref_wl
        self.coating_wvl = 550.0

    @property
    def central_wvl(self):
        return self.wavelengths[self.reference_wvl]

    @central_wvl.setter
    def central_wvl(self, wvl):
        self.wavelengths[self.reference_wvl] = wvl

    def set_from_list(self, wlwts):
        self.wavelengths = []
        self.spectral_wts = []
        for wlwt in wlwts:
            self.wavelengths.append(get_wavelength(wlwt[0]))
            self.spectral_wts.append(wlwt[1])
        self.calc_colors()

    def sync_to_restore(self, optical_spec):
        self.calc_colors()

    def set_from_specsheet(self, ss):
        pass

    def update_model(self, **kwargs):
        self.calc_colors()

    def add(self, wl, wt):
        self.wavelengths.append(get_wavelength(wl))
        self.spectral_wts.append(wt)
        self.spectrum.sort(key=lambda w: w[0], reverse=True)

    def calc_colors(self):
        accent = colors.accent_colors()
        self.render_colors = []
        num_wvls = len(self.wavelengths)
        if num_wvls == 1:
            self.render_colors.append(accent['green'])
        elif num_wvls > 1:
            step = 1 if self.wavelengths[0] < self.wavelengths[-1] else -1
            if num_wvls == 2:
                c = ['blue', 'red']
            elif num_wvls == 3:
                c = ['blue', 'green', 'red']
            elif num_wvls == 4:
                c = ['blue', 'green', 'yellow', 'red']
            elif num_wvls == 5:
                c = ['violet', 'cyan', 'green', 'yellow', 'red']
            elif num_wvls == 6:
                c = ['violet', 'cyan', 'green', 'yellow', 'red', 'magenta']
            else:
                c = ['violet', 'blue', 'cyan', 'green', 'yellow',
                     'red', 'magenta']
            self.render_colors = [accent[clr] for clr in c[::step]]
        # else:
        #     for w in self.wavelengths:
        #         print("calc_colors", w)
        #         rgb = srgb.wvl_to_rgb(w)
        #         print("rgb", rgb)
        #         self.render_colors.append(rgb)


class PupilSpec:
    """ Aperture specification

    Attributes:
        key: 'aperture', 'object'|'image', 'pupil'|'NA'|'f/#'
        value: size of the pupil
        pupil_rays: list of relative pupil coordinates for pupil limiting rays
        ray_labels: list of string labels for pupil_rays
    """
    default_pupil_rays = [[0., 0.], [1., 0.], [-1., 0.], [0., 1.], [0., -1.]]
    default_ray_labels = ['00', '+X', '-X', '+Y', '-Y']

    def __init__(self, parent, key=('object', 'pupil'), value=1.0):
        self.optical_spec = parent
        self.key = 'aperture', key[0], key[1]
        self.value = value
        self.pupil_rays = PupilSpec.default_pupil_rays
        self.ray_labels = PupilSpec.default_ray_labels

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['optical_spec']
        return attrs

    def sync_to_restore(self, optical_spec):
        if hasattr(self, 'pupil_type'):
            self.key = model_enums.get_ape_key_for_type(self.pupil_type)
            del self.pupil_type
        self.optical_spec = optical_spec

    def set_from_specsheet(self, ss):
        self.key, self.value = ss.get_etendue_inputs('aperture')

    def get_input_for_specsheet(self):
        return self.key, self.value

    def update_model(self, **kwargs):
        if not hasattr(self, 'pupil_rays'):
            self.pupil_rays = PupilSpec.default_pupil_rays
            self.ray_labels = PupilSpec.default_ray_labels

    def get_pupil_type(self):
        return model_enums.get_ape_type_for_key(self.key).value

    def mutate_pupil_type(self, new_pupil_type):
        ape_key = model_enums.get_ape_key_for_type(new_pupil_type)
        aperture, obj_img_key, value_key = ape_key
        if self.optical_spec is not None:
            if self.optical_spec.parax_data is not None:
                fod = self.optical_spec.parax_data.fod
                if obj_img_key == 'object':
                    if value_key == 'pupil':
                        self.value = 2*fod.enp_radius
                    elif value_key == 'NA':
                        self.value = fod.obj_na
                elif obj_img_key == 'image':
                    if value_key == 'f/#':
                        self.value = fod.fno
                    elif value_key == 'NA':
                        self.value = fod.img_na

        self.key = ape_key


class FieldSpec:
    """ Field of view specification

    Attributes:
        key: 'field', 'object'|'image', 'height'|'angle'
        value: maximum field, per the key
        fields: list of Field instances
        is_relative: if True, `fields` are relative to max field

    """

    def __init__(self, parent, key=('object', 'angle'), value=0., flds=[0.],
                 is_relative=False, do_init=True, **kwargs):
        self.optical_spec = parent
        self.key = 'field', key[0], key[1]
        self.value = value
        self.is_relative = is_relative
        if do_init:
            self.set_from_list(flds)
        else:
            self.fields = []

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['optical_spec']
        return attrs

    def sync_to_restore(self, optical_spec):
        if hasattr(self, 'field_type'):
            self.key = model_enums.get_fld_key_for_type(self.field_type)
            del self.field_type
        if not hasattr(self, 'is_relative'):
            self.is_relative = False
        if not hasattr(self, 'value'):
            self.value, _ = self.max_field()
        self.optical_spec = optical_spec

    def __str__(self):
        return "key={}, max field={}".format(self.key, self.max_field()[0])

    def set_from_list(self, flds):
        self.fields = [Field() for f in range(len(flds))]
        for i, f in enumerate(self.fields):
            f.y = flds[i]
        self.value, _ = self.max_field()

    def set_from_specsheet(self, ss):
        key, value = ss.get_etendue_inputs('field')
        if value != 0 and len(self.fields) == 1:
            # just one field, add a second one for max value
            self.is_relative = True
            self.fields.append(Field(x=0, y=1))

        if not self.is_relative:
            fld_scale = 1 if self.value == 0 else value/self.value

            for i, f in enumerate(self.fields):
                f.x *= fld_scale
                f.y *= fld_scale

        self.key, self.value = key, value

    def get_input_for_specsheet(self):
        return self.key, self.value

    def update_model(self, **kwargs):
        for f in self.fields:
            f.update()

        # recalculate max_field and relabel fields.
        #  relabeling really assumes the fields are radial, specifically,
        #  y axis only
        if self.is_relative:
            field_norm = 1
        else:
            field_norm = 1 if self.value == 0 else 1.0/self.value

        self.index_labels = []
        for f in self.fields:
            if f.x != 0.0:
                fldx = '{:5.2f}x'.format(field_norm*f.x)
            else:
                fldx = ''
            if f.y != 0.0:
                fldy = '{:5.2f}y'.format(field_norm*f.y)
            else:
                fldy = ''
            self.index_labels.append(fldx + fldy)
        self.index_labels[0] = 'axis'
        if len(self.index_labels) > 1:
            self.index_labels[-1] = 'edge'
        return self

    def get_field_type(self):
        return model_enums.get_fld_type_for_key(self.key).value

    def mutate_field_type(self, new_field_type):
        osp = self.optical_spec
        fld_key = model_enums.get_fld_key_for_type(new_field_type)
        field, obj_img_key, value_key = fld_key
        if self.optical_spec is not None:
            if osp.parax_data is not None:
                fod = self.optical_spec.parax_data.fod
                if obj_img_key == 'object':
                    if value_key == 'height':
                        self.value = osp.parax_data.pr_ray[0][mc.ht]
                    elif value_key == 'angle':
                        self.value = fod.obj_ang
                elif obj_img_key == 'image':
                    if value_key == 'height':
                        self.value = fod.img_ht
        self.key = fld_key

    def obj_coords(self, fld):
        fld_coord = np.array([fld.x, fld.y, 0.0])
        if self.is_relative:
            fld_coord *= self.value

        field, obj_img_key, value_key = self.key

        fod = self.optical_spec.parax_data.fod
        if obj_img_key == 'object':
            if value_key == 'angle':
                dir_tan = np.tan(np.deg2rad(fld_coord))
                obj_pt = -dir_tan*(fod.obj_dist+fod.enp_dist)
            elif value_key == 'height':
                obj_pt = fld_coord
        elif obj_img_key == 'image':
            if value_key == 'height':
                img_pt = fld_coord
                obj_pt = fod.red*img_pt
        return obj_pt

    def max_field(self):
        """ calculates the maximum field of view

        Returns:
            magnitude of maximum field, maximum Field instance
        """
        max_fld = None
        max_fld_sqrd = -1.0
        for i, f in enumerate(self.fields):
            fld_sqrd = f.x*f.x + f.y*f.y
            if fld_sqrd > max_fld_sqrd:
                max_fld_sqrd = fld_sqrd
                max_fld = i
        max_fld_value = math.sqrt(max_fld_sqrd)
        if self.is_relative:
            max_fld_value *= self.value
        return max_fld_value, max_fld


class Field:
    """ a single field point, largely a data container

    Attributes:
        x: x field component
        y: y field component
        vux: +x vignetting factor
        vuy: +y vignetting factor
        vlx: -x vignetting factor
        vly: -y vignetting factor
        wt: field weight
        aim_pt: x, y chief ray coords on the paraxial entrance pupil plane
        chief_ray: ray package for the ray from the field point throught the
                   center of the aperture stop, traced in the central
                   wavelength
        ref_sphere: a tuple containing (image_pt, ref_dir, ref_sphere_radius)

    """

    def __init__(self, x=0., y=0., wt=1.):
        self.x = x
        self.y = y
        self.vux = 0.0
        self.vuy = 0.0
        self.vlx = 0.0
        self.vly = 0.0
        self.wt = wt
        self.aim_pt = None
        self.chief_ray = None
        self.ref_sphere = None

    def __json_encode__(self):
        attrs = dict(vars(self))
        items = ['chief_ray', 'ref_sphere', 'pupil_rays']
        for item in items:
            if item in attrs:
                del attrs[item]
        return attrs

    def __str__(self):
        return "{}, {}".format(self.x, self.y)

    def __repr__(self):
        return "Field(x={}, y={}, wt={})".format(self.x, self.y, self.wt)

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
    """ Focus range specification

    Attributes:
        focus_shift: focus shift (z displacement) from nominal image interface
        defocus_range: +/- half the total focal range, from the focus_shift
                       position
    """

    def __init__(self, focus_shift=0.0, defocus_range=0.0):
        self.focus_shift = focus_shift
        self.defocus_range = defocus_range

    def __repr__(self):
        return ("FocusRange(focus_shift={}, defocus_range={})"
                .format(self.focus_shift, self.defocus_range))

    def set_from_specsheet(self, ss):
        pass

    def update(self):
        pass

    def get_focus(self, fr=0.0):
        """ return focus position for input focus range parameter

        Args:
            fr (float): focus range parameter, -1.0 to 1.0

        Returns:
            focus position for input focus range parameter
        """
        return self.focus_shift + fr*self.defocus_range
