#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Top level model classes

.. Created on Wed Mar 14 11:08:28 2018

.. codeauthor: Michael J. Hayford
"""

import os.path
import json_tricks

import rayoptics.codev.cmdproc as cvp

from rayoptics.optical.model_enums import DimensionType as dt

from rayoptics.optical.elements import ElementModel
from rayoptics.optical.paraxialdesign import ParaxialModel
from rayoptics.optical.sequential import SequentialModel
from rayoptics.optical.opticalspec import OpticalSpecs


def open_model(file_name):
    """ open a file and populate an optical model with the data

    Args:
        file_name (str): a filename of a supported file type

            - .roa - a rayoptics JSON encoded file
            - .seq - a CODE V (TM) sequence file

    Returns:
        if successful, an OpticalModel instance, otherwise, None
    """
    file_extension = os.path.splitext(file_name)[1]
    opm = None
    if file_extension == '.seq':
        opm = OpticalModel()
        cvp.read_lens(opm, file_name)
    elif file_extension == '.roa':
        with open(file_name, 'r') as f:
            obj_dict = json_tricks.load(f)
            if 'optical_model' in obj_dict:
                opm = obj_dict['optical_model']
                opm.sync_to_restore()
    return opm


class SystemSpec:
    """ Container for units and other system level constants

    Attributes:
        title (str): a short description of the model
        initials (str): user initials or other id
        dimensions: a DimensionType enum of the model linear units
        temperature (double): model temperature in degrees Celsius
        pressure (double): model pressure in mm/Hg
    """
    def __init__(self):
        self.title = ''
        self.initials = ''
        self.dimensions = dt.MM
        self.temperature = 20.0
        self.pressure = 760.0

    def nm_to_sys_units(self, nm):
        """ convert nm to system units

        Args:
            nm (double): value in nm

        Returns:
            double: value converted to system units
        """
        if self.dimensions == dt.M:
            return 1e-9 * nm
        elif self.dimensions == dt.CM:
            return 1e-7 * nm
        elif self.dimensions == dt.MM:
            return 1e-6 * nm
        elif self.dimensions == dt.IN:
            return 1e-6 * nm/25.4
        elif self.dimensions == dt.FT:
            return 1e-6 * nm/304.8
        else:
            return nm


class OpticalModel:
    """ Top level container for optical model.

    The OpticalModel serves as a top level container of model properties.
    Key aspects are built-in element and surface based repesentations of the
    optical surfaces.
    A sequential optical model is a sequence of surfaces and gaps.
    Additionally, it includes optical usage information to specify the
    aperture, field of view, spectrum and focus.
    """

    def __init__(self):
        self.radius_mode = False
        self.system_spec = SystemSpec()
        self.seq_model = SequentialModel(self)
        self.optical_spec = OpticalSpecs(self)
        self.parax_model = ParaxialModel(self)
        self.ele_model = ElementModel(self)

    def name(self):
        return self.system_spec.title

    def reset(self):
        rdm = self.radius_mode
        self.__init__()
        self.radius_mode = rdm

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['parax_model']
        return attrs

    def save_model(self, file_name):
        file_extension = os.path.splitext(file_name)[1]
        filename = file_name if len(file_extension) > 0 else file_name+'.roa'
        fs_dict = {}
        fs_dict['optical_model'] = self
        with open(filename, 'w') as f:
            json_tricks.dump(fs_dict, f, indent=1,
                             separators=(',', ':'))

    def sync_to_restore(self):
        self.seq_model.sync_to_restore(self)
        self.ele_model.sync_to_restore(self)
        self.optical_spec.sync_to_restore(self)
        if not hasattr(self, 'parax_model'):
            self.parax_model = ParaxialModel(self)
        else:
            self.parax_model.sync_to_restore(self)
        self.update_model()

    def update_model(self):
        self.seq_model.update_model()
        self.optical_spec.update_model()
        self.parax_model.update_model()
        self.ele_model.update_model()

    def nm_to_sys_units(self, nm):
        """ convert nm to system units

        Args:
            nm (double): value in nm

        Returns:
            double: value converted to system units
        """
        return self.system_spec.nm_to_sys_units(nm)
