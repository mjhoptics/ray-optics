#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Top level model classes

.. Created on Wed Mar 14 11:08:28 2018

.. codeauthor: Michael J. Hayford
"""
import os.path
import json_tricks

import rayoptics

from rayoptics.optical import elements
import rayoptics.optical.model_constants as mc

from rayoptics.optical.paraxialdesign import ParaxialModel
from rayoptics.optical.sequential import SequentialModel
from rayoptics.optical.opticalspec import OpticalSpecs
from rayoptics.optical.specsheet import create_specsheet_from_model
from rayoptics.optical.model_enums import DimensionType as dt


class SystemSpec:
    """ Container for units and other system level constants

    Attributes:
        title (str): a short description of the model
        initials (str): user initials or other id
        dimensions: a DimensionType enum of the model linear units
        temperature (float): model temperature in degrees Celsius
        pressure (float): model pressure in mm/Hg
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
            nm (float): value in nm

        Returns:
            float: value converted to system units
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

    def __init__(self, radius_mode=False, specsheet=None):
        self.ro_version = rayoptics.__version__
        self.radius_mode = radius_mode
        self.specsheet = specsheet
        self.system_spec = SystemSpec()
        self.seq_model = SequentialModel(self)
        self.optical_spec = OpticalSpecs(self)
        self.parax_model = ParaxialModel(self)
        self.ele_model = elements.ElementModel(self)

        if self.specsheet:
            self.set_from_specsheet()

    def name(self):
        return self.system_spec.title

    def reset(self):
        rdm = self.radius_mode
        self.__init__()
        self.radius_mode = rdm

    def __json_encode__(self):
        attrs = dict(vars(self))
        if hasattr(self, 'app_manager'):
            del attrs['app_manager']
        return attrs

    def set_from_specsheet(self, specsheet=None):
        if specsheet:
            self.specsheet = specsheet
        else:
            specsheet = self.specsheet
        self.seq_model.set_from_specsheet(specsheet)
        self.optical_spec.set_from_specsheet(specsheet)

    def save_model(self, file_name):
        file_extension = os.path.splitext(file_name)[1]
        filename = file_name if len(file_extension) > 0 else file_name+'.roa'
        fs_dict = {}
        fs_dict['optical_model'] = self
        with open(filename, 'w') as f:
            json_tricks.dump(fs_dict, f, indent=1,
                             separators=(',', ':'), allow_nan=True)

    def sync_to_restore(self):
        if not hasattr(self, 'ro_version'):
            self.ro_version = rayoptics.__version__

        self.seq_model.sync_to_restore(self)
        self.ele_model.sync_to_restore(self)
        self.optical_spec.sync_to_restore(self)

        if hasattr(self, 'parax_model'):
            self.parax_model.sync_to_restore(self)
        else:
            self.parax_model = ParaxialModel(self)

        if hasattr(self, 'specsheet'):
            self.specsheet.sync_to_restore(self)
        else:
            self.specsheet = None

        self.update_model()

    def update_model(self):
        self.seq_model.update_model()
        self.optical_spec.update_model()
        self.parax_model.update_model()
        self.ele_model.update_model()
        if self.specsheet is None:
            self.specsheet = create_specsheet_from_model(self)

    def nm_to_sys_units(self, nm):
        """ convert nm to system units

        Args:
            nm (float): value in nm

        Returns:
            float: value converted to system units
        """
        return self.system_spec.nm_to_sys_units(nm)

    def add_lens(self, **kwargs):
        seq, ele = elements.create_lens(**kwargs)
        kwargs['insert'] = True
        self.insert_ifc_gp_ele(seq, ele, **kwargs)

    def add_mirror(self, **kwargs):
        seq, ele = elements.create_mirror(**kwargs)
        kwargs['insert'] = True
        self.insert_ifc_gp_ele(seq, ele, **kwargs)

    def add_thinlens(self, **kwargs):
        seq, ele = elements.create_thinlens(**kwargs)
        kwargs['insert'] = True
        self.insert_ifc_gp_ele(seq, ele, **kwargs)

    def add_dummy_plane(self, **kwargs):
        seq, ele = elements.create_dummy_plane(**kwargs)
        kwargs['insert'] = True
        self.insert_ifc_gp_ele(seq, ele, **kwargs)

    def add_from_file(self, filename, **kwargs):
        seq, ele = elements.create_from_file(filename, **kwargs)
        kwargs['insert'] = True
        self.insert_ifc_gp_ele(seq, ele, **kwargs)

    def insert_ifc_gp_ele(self, *args, **kwargs):
        """ insert interfaces and gaps into seq_model and eles into ele_model
        """
        seq, ele = args
        if 'idx' in kwargs:
            self.seq_model.cur_surface = kwargs['idx']

        # distinguish between adding a new chunk, which requires splitting a
        #  gap in two, and replacing a node, which uses the existing gaps.
        if 'insert' in kwargs:
            t = kwargs['t'] if 't' in kwargs else 0.
            g, ag = elements.create_air_gap(t=t, ref_ifc=seq[-1][mc.Intfc])
            seq[-1][mc.Gap] = g
            ele.append(ag)
        else:
            # replacing an existing node. need to hook new chunk final
            # interface to the existing gap and following (air gap) element
            g = self.seq_model.gaps[self.seq_model.cur_surface+1]
            seq[-1][mc.Gap] = g
            ag = self.ele_model.gap_dict[g]
            ag.ref_ifc = seq[-1][mc.Intfc]  # tacit assumption is ag == AirGap

        for sg in seq:
            self.seq_model.insert(sg[mc.Intfc], sg[mc.Gap])

        for e in ele:
            self.ele_model.add_element(e)
        self.ele_model.sequence_elements()

    def remove_ifc_gp_ele(self, *args, **kwargs):
        """ remove interfaces and gaps from seq_model and eles from ele_model
        """
        seq, ele = args
        sg = seq[0]
        idx = self.seq_model.ifcs.index(sg[mc.Intfc])

        # verify that the sequences match
        seq_match = True
        for i, sg in enumerate(seq):
            if sg[0] is not self.seq_model.ifcs[idx+i]:
                seq_match = False
                break

        if seq_match:
            # remove interfaces in reverse
            for i in range(idx+len(seq)-1, idx-1, -1):
                self.seq_model.remove(i)

        for e in ele:
            self.ele_model.remove_element(e)
