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

import rayoptics.elem.elements as ele
import rayoptics.optical.model_constants as mc

from rayoptics.elem.elements import ElementModel
from rayoptics.elem.parttree import (PartTree, elements_from_sequence)
from rayoptics.parax.paraxialdesign import ParaxialModel
from rayoptics.seq.sequential import SequentialModel
from rayoptics.raytr.opticalspec import OpticalSpecs
from rayoptics.parax.specsheet import create_specsheet_from_model
from rayoptics.optical.model_enums import get_dimension_for_type


class SystemSpec:
    """ Container for units and other system level constants

    Attributes:
        title (str): a short description of the model
        initials (str): user initials or other id
        temperature (float): model temperature in degrees Celsius
        pressure (float): model pressure in mm/Hg
    """

    def __init__(self):
        self.title = ''
        self.initials = ''
        self.dimensions = 'mm'
        self.temperature = 20.0
        self.pressure = 760.0

    def __json_decode__(self, **attrs):
        for a_key, a_val in attrs.items():
            if a_key == 'dimensions':
                self._dimensions = (a_val if isinstance(a_val, str)
                                    else get_dimension_for_type(a_val))
            else:
                setattr(self, a_key, a_val)

    @property
    def dimensions(self):
        """ the model linear units (str). """
        return self._dimensions

    @dimensions.setter
    def dimensions(self, value):
        self._dimensions = (value if isinstance(value, str)
                            else get_dimension_for_type(value))

    def nm_to_sys_units(self, nm):
        """ convert nm to system units

        Args:
            nm (float): value in nm

        Returns:
            float: value converted to system units
        """
        if self.dimensions == 'm':
            return 1e-9 * nm
        elif self.dimensions == 'cm':
            return 1e-7 * nm
        elif self.dimensions == 'mm':
            return 1e-6 * nm
        elif self.dimensions == 'in':
            return 1e-6 * nm/25.4
        elif self.dimensions == 'ft':
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

    Attributes:
        ro_version: current version of rayoptics
        radius_mode: if True output radius, else output curvature
        specsheet: :class:`~rayoptics.parax.specsheet.SpecSheet`
        system_spec: :class:`.SystemSpec`
        seq_model: :class:`~rayoptics.seq.sequential.SequentialModel`
        optical_spec: :class:`~rayoptics.raytr.opticalspec.OpticalSpecs`
        parax_model: :class:`~rayoptics.parax.paraxialdesign.ParaxialModel`
        ele_model: :class:`~rayoptics.elem.elements.ElementModel`
    """

    def __init__(self, radius_mode=False, specsheet=None, **kwargs):
        self.ro_version = rayoptics.__version__
        self.radius_mode = radius_mode
        
        self.specsheet = specsheet
        self.system_spec = SystemSpec()
        self.seq_model = SequentialModel(self, **kwargs)
        self.optical_spec = OpticalSpecs(self, specsheet=specsheet, **kwargs)
        self.parax_model = ParaxialModel(self, **kwargs)
        self.ele_model = ElementModel(self, **kwargs)
        self.part_tree = PartTree(self, **kwargs)

        self.map_submodels()

        if self.specsheet:
            self.set_from_specsheet()

        if kwargs.get('do_init', True):
            # need to do this after OpticalSpec is initialized
            self.seq_model.update_model()
            elements_from_sequence(self.ele_model,
                                   self.seq_model,
                                   self.part_tree)

    def map_submodels(self):
        """Setup machinery for model mapping api. """
        submodels = {}
        submodels['specsheet'] = self.specsheet
        submodels['system_spec'] = self.system_spec
        submodels['seq_model'] = self.seq_model
        submodels['optical_spec'] = self.optical_spec
        submodels['parax_model'] = self.parax_model
        submodels['ele_model'] = self.ele_model
        submodels['part_tree'] = self.part_tree
        # Add a level of indirection to allow short and long aliases
        submodel_aliases = {
            'ss': 'specsheet', 'specsheet': 'specsheet',
            'sys': 'system_spec', 'system_spec': 'system_spec',
            'sm': 'seq_model', 'seq_model': 'seq_model',
            'osp': 'optical_spec', 'optical_spec': 'optical_spec',
            'pm': 'parax_model', 'parax_model': 'parax_model',
            'em': 'ele_model', 'ele_model': 'ele_model',
            'pt': 'part_tree', 'part_tree': 'part_tree',
            }
        self._submodels = submodels, submodel_aliases

    def __getitem__(self, key):
        """ Provide mapping interface to submodels. """
        submodels, submodel_aliases = self._submodels
        return submodels[submodel_aliases[key]]

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
        del attrs['_submodels']
        return attrs

    def set_from_specsheet(self, specsheet=None):
        if specsheet:
            self.specsheet = specsheet
        else:
            specsheet = self.specsheet
        self.optical_spec.set_from_specsheet(specsheet)
        self.seq_model.set_from_specsheet(specsheet)

    def save_model(self, file_name, version=None):
        """Save the optical_model in a ray-optics JSON file.
        
        Args:
            file_name: str or Path
            version: optional override for rayoptics version number
        """
        file_extension = os.path.splitext(file_name)[1]
        filename = file_name if len(file_extension) > 0 else file_name+'.roa'
        # update version number prior to writing file.
        self.ro_version = rayoptics.__version__ if version is None else version

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

        if hasattr(self, 'part_tree'):
            self.part_tree.sync_to_restore(self)
        else:
            self.part_tree = PartTree(self)
            self.part_tree.add_element_model_to_tree(self.ele_model)

        self.map_submodels()

        self.update_model()

    def update_model(self, **kwargs):
        self.seq_model.update_model(**kwargs)
        self.optical_spec.update_model(**kwargs)
        self.parax_model.update_model(**kwargs)
        self.ele_model.update_model(**kwargs)
        self.part_tree.update_model(**kwargs)
        if self.specsheet is None:
            self.specsheet = create_specsheet_from_model(self)
            self.map_submodels()

    def nm_to_sys_units(self, nm):
        """ convert nm to system units

        Args:
            nm (float): value in nm

        Returns:
            float: value converted to system units
        """
        return self.system_spec.nm_to_sys_units(nm)

    def add_lens(self, **kwargs):
        descriptor = ele.create_lens(**kwargs)
        kwargs['insert'] = True
        self.insert_ifc_gp_ele(*descriptor, **kwargs)

    def add_mirror(self, **kwargs):
        descriptor = ele.create_mirror(**kwargs)
        kwargs['insert'] = True
        self.insert_ifc_gp_ele(*descriptor, **kwargs)

    def add_thinlens(self, **kwargs):
        descriptor = ele.create_thinlens(**kwargs)
        kwargs['insert'] = True
        self.insert_ifc_gp_ele(*descriptor, **kwargs)

    def add_dummy_plane(self, **kwargs):
        descriptor = ele.create_dummy_plane(**kwargs)
        kwargs['insert'] = True
        self.insert_ifc_gp_ele(*descriptor, **kwargs)

    def add_from_file(self, filename, **kwargs):
        descriptor = ele.create_from_file(filename, **kwargs)
        kwargs['insert'] = True
        self.insert_ifc_gp_ele(*descriptor, **kwargs)

    def insert_ifc_gp_ele(self, *descriptor, **kwargs):
        """ insert interfaces and gaps into seq_model and eles into ele_model

        Args:
            descriptor: a tuple of additions for the sequential, element and
                        part tree models
            kwargs: keyword arguments including
                idx: insertion point in the sequential model
                insert: if True, insert the chunk, otherwise replace it
                t: the thickness following a chuck when inserting
        """
        sm = self['seq_model']
        seq, elm, e_node = descriptor
        if 'idx' in kwargs:
            sm.cur_surface = kwargs['idx']
        idx = sm.cur_surface

        e_node.parent = self.part_tree.root_node

        # distinguish between adding a new chunk, which requires splitting a
        #  gap in two, and replacing a node, which uses the existing gaps.
        ins_prev_gap = False
        if 'insert' in kwargs:
            t_after = kwargs['t'] if 't' in kwargs else 0.
            if sm.get_num_surfaces() == 2:
                # only object space gap, add image space gap following this
                gap_label = "Image space"
                ins_prev_gap = False
            else:
                # we have both object and image space gaps; retain the image
                # space gap by splitting and inserting the new gap before the
                # inserted chunk, unless we're inserting before idx=1.
                gap_label = None
                if idx > 0:
                    ins_prev_gap = True

            if ins_prev_gap:
                t_air, sm.gaps[idx].thi = sm.gaps[idx].thi, t_after
            else:
                t_air = t_after
            g, ag, ag_node = ele.create_air_gap(t=t_air, label=gap_label)
            if not ins_prev_gap:
                seq[-1][mc.Gap] = g
            elm.append(ag)
            ag_node.parent = self.part_tree.root_node
        else:
            # replacing an existing node. need to hook new chunk final
            # interface to the existing gap and following (air gap) element
            g = sm.gaps[sm.cur_surface+1]
            seq[-1][mc.Gap] = g
            ag, ag_node = self.part_tree.parent_object(g, '#airgap')            # ag.idx = seq[-1][mc.Intfc]

        for sg in seq:
            if ins_prev_gap:
                gap, g = g, sg[mc.Gap]
            else:
                gap = sg[mc.Gap]
            sm.insert(sg[mc.Intfc], gap, prev=ins_prev_gap)

        for e in elm:
            self.ele_model.add_element(e)
        self.ele_model.sequence_elements()

    def remove_ifc_gp_ele(self, *descriptor, **kwargs):
        """ remove interfaces and gaps from seq_model and eles from ele_model
        """
        seq, elm, e_node = descriptor
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

        for e in elm:
            self.ele_model.remove_element(e)

        e_node.parent = None

    def remove_node(self, e_node):
        # remove interfaces from seq_model
        self.seq_model.remove_node(e_node)
        # remove elements from ele_model
        self.ele_model.remove_node(e_node)
        # unhook node
        e_node.parent = None

    def rebuild_from_seq(self):
        """ Rebuild ele_model and part_tree from seq_model. """
        self['em'].elements = []
        self['pt'].root_node.children = []
        elements_from_sequence(self['em'], self['sm'], self['pt'])
