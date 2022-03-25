#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Top level model classes

.. Created on Wed Mar 14 11:08:28 2018

.. codeauthor: Michael J. Hayford
"""
import os.path
import json_tricks
from collections.abc import Sequence

import rayoptics

import rayoptics.elem.elements as ele
import rayoptics.optical.model_constants as mc

from rayoptics.elem.elements import (ElementModel, do_flip_with_part_list)
from rayoptics.elem import parttree
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

    def listobj_str(self):
        vs = vars(self)
        o_str = f"{type(self).__name__}:\n"
        for k, v in vs.items():
             o_str += f"{k}: {v}\n"
        return o_str

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
        self.analysis_results = {'parax_data': None}

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
        submodels['analysis_results'] = self.analysis_results
        # Add a level of indirection to allow short and long aliases
        submodel_aliases = {
            'ss': 'specsheet', 'specsheet': 'specsheet',
            'sys': 'system_spec', 'system_spec': 'system_spec',
            'sm': 'seq_model', 'seq_model': 'seq_model',
            'osp': 'optical_spec', 'optical_spec': 'optical_spec',
            'pm': 'parax_model', 'parax_model': 'parax_model',
            'em': 'ele_model', 'ele_model': 'ele_model',
            'pt': 'part_tree', 'part_tree': 'part_tree',
            'ar': 'analysis_results', 'analysis_results': 'analysis_results',
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
        # not sure about saving analysis_results...
        if hasattr(self, 'analysis_results'):
            del attrs['analysis_results']
        del attrs['_submodels']
        return attrs

    def listobj_str(self):
        vs = vars(self)
        o_str = f"{type(self).__name__}:\n"
        for k, v in vs.items():
             o_str += f"{k}: {v}\n"
        return o_str

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

        self.profile_dict = self._build_profile_dict()
        self.parts_dict = {id(p):p for p in self.ele_model.elements}

        fs_dict = {}
        fs_dict['optical_model'] = self
        with open(filename, 'w') as f:
            json_tricks.dump(fs_dict, f, indent=1,
                             separators=(',', ':'), allow_nan=True)
        delattr(self, 'profile_dict')
        delattr(self, 'parts_dict')

    def _build_profile_dict(self):
        """ build a profile dict for the union of the seq_model and part_tree. """
        profile_dict = {}
        for ifc in self.seq_model.ifcs:
            if ifc.profile is not None:
                profile_dict[str(id(ifc.profile))] = ifc.profile
        profile_nodes = self.part_tree.nodes_with_tag(tag='#profile')
        for profile_node in profile_nodes:
            profile = profile_node.id
            profile_id = str(id(profile))
            if profile_id not in profile_dict:
                profile_dict[profile_id] = profile
                print(f"found new profile in part_tree: "
                      f"{profile_node.parent.name}.{profile_node.name}")
        return profile_dict

    def sync_to_restore(self):
        if not hasattr(self, 'ro_version'):
            self.ro_version = rayoptics.__version__

        self.profile_dict = (self.profile_dict if hasattr(self, 'profile_dict')
                             else {})
        self.parts_dict = (self.parts_dict if hasattr(self, 'parts_dict')
                             else {})
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

        if not hasattr(self, 'analysis_results'):
            self.analysis_results = {'parax_data': None}

        self.map_submodels()

        self.update_model()

        # delete the profile_dict used for save/restore fidelity
        # the idea is one instance of an object is saved to a dictionary,
        #  keyed to the original object's id. each class that uses the original
        #  object save the object id, instead of the object. During the restore
        #  process, sync_to_restore can use the saved profile id to lookup
        #  the restored profile in the profile_dict.
        delattr(self, 'profile_dict')
        delattr(self, 'parts_dict')

    def update_model(self, **kwargs):
        """Model and its constituents are updated.

        Args:
            kwargs: possible keyword arguments including:

                - build:
    
                    - 'rebuild': rebuild the model "from scratch", e.g number of nodes changes
                    - 'update': number of nodes unchanged, just the parameters
    
                - src_model: model that originated the modification

        """
        self['seq_model'].update_model(**kwargs)
        self['optical_spec'].update_model(**kwargs)
        self.update_optical_properties(**kwargs)
        if len(self['part_tree'].nodes_with_tag(tag='#element')) == 0:
            elements_from_sequence(self['ele_model'], self['seq_model'], 
                                   self['part_tree'])

        self['ele_model'].update_model(**kwargs)
        self['part_tree'].update_model(**kwargs)
        if self.specsheet is None:
            self.specsheet = create_specsheet_from_model(self)
            self.map_submodels()

    def update_optical_properties(self, **kwargs):
        """Compute first order and other optical properties. """

        # OpticalSpec maintains first order and ray aiming for fields
        self['optical_spec'].update_optical_properties(**kwargs)
        # Update the ParaxialModel as needed
        self['parax_model'].update_model(**kwargs)
        # Update surface apertures, if requested (do_apertures=True)
        self['seq_model'].update_optical_properties(**kwargs)

    def nm_to_sys_units(self, nm):
        """ convert nm to system units

        Args:
            nm (float): value in nm

        Returns:
            float: value converted to system units
        """
        return self.system_spec.nm_to_sys_units(nm)

    def add_lens(self, **kwargs):
        """ Add a lens into the optical model

        Args:
            kwargs: keyword arguments including:

                - idx: insertion point in the sequential model
                - t: the thickness following a chunk when inserting
                - lens: tuple of `cv1, cv2, th, glass_name_catalog, sd` where:

                    - cv1: front curvature
                    - cv2: rear curvature
                    - th: lens thickness
                    - glass_name_catalog: a str, e.g. 'N-BK7, Schott'
                    - sd: lens semi-diameter

        """
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

    def add_assembly_from_seq(self, idx1, idx2, **kwargs):
        """ Create an Assembly from the elements in the sequence range. """
        pt = self['part_tree']
        ifc1 = self['seq_model'].ifcs[idx1]
        # record where to insert asm_node in root_node.children
        start_idx = pt.root_node.children.index(pt.node(ifc1).ancestors[1])
        asm, asm_node = ele.create_assembly_from_seq(self, idx1, idx2, **kwargs)
        self['ele_model'].add_element(asm)
        # use anytree mechanism to add asm_node to tree
        asm_node.parent = pt.root_node
        root_children = list(pt.root_node.children)
        # move asm_node to desired spot and update root's children
        root_children.remove(asm_node)
        root_children.insert(start_idx, asm_node)
        pt.root_node.children = root_children

    def rebuild_from_seq(self):
        """ Rebuild ele_model and part_tree from seq_model. 
        
        When in doubt about whether there is a problem with bad data in an 
        OpticalModel, this function can be used to rebuild everything from 
        the sequential model.
        
        """
        self['em'].elements = []
        self['pt'].root_node.children = []
        elements_from_sequence(self['em'], self['sm'], self['pt'])

    def flip(self, *args, **kwargs):
        """ Flip a `Part` or an `Interface` range in the optical model. 
        
        The flip operation supports several different ways of specifying what is
        to be flipped.

        Args:

        - None: This flips the model from `ifc 1` to `ifc image-1`
        - idx1, idx2: This flips the model between interfaces idx1 and idx2. Flipping from object to image is disallowed.
        - part: This flips the corresponding part in the model
        - list: This flips a list of parts in the model

        """
        import numpy as np
        rot_around_x = np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]])
        rot_around_y = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
        sm = self['seq_model']
        em = self['ele_model']
        pt = self['part_tree']

        if len(args) == 0:
            # default behavior: flip 1st to image-1 interfaces
            args = 1, len(sm.gaps)-1

        if isinstance(args[0], int):
            # flip a range of interfaces/gaps
            # disallow flipping from object to image
            idx1 = args[0] if args[0] > 0 else 1
            idx2 = args[1] if args[1] < len(sm.gaps) else len(sm.gaps)-1
            part_list, node_list = parttree.part_list_from_seq(self,
                                                               idx1, idx2)
            flip_pt = 0.5*(sm.gbl_tfrms[idx2][1] + sm.gbl_tfrms[idx1][1])
            flip_pt_tfrm = sm.gbl_tfrms[idx1][0], flip_pt
            do_flip_with_part_list(part_list, flip_pt_tfrm)

        elif isinstance(args[0], list):
            # flip a list of parts
            part_list = args[0]
            idxs = []
            for p in part_list:
                idxs += p.idx_list()
            idx1, idx2 = idxs[0], idxs[-1]

            flip_pt = 0.5*(sm.gbl_tfrms[idx2][1] + sm.gbl_tfrms[idx1][1])
            flip_pt_tfrm = sm.gbl_tfrms[idx1][0], flip_pt
            do_flip_with_part_list(part_list, flip_pt_tfrm)

        elif isinstance(args[0], ele.Part):
            # flip a Part, any subtype
            p = args[0]
            p.flip()
            idx_list = p.idx_list()
            idx1 = idx_list[0]
            idx2 = idx_list[-1]

        # flip the range in the sequential model
        sm.flip(idx1, idx2)
        self.update_optical_properties()

        pt.update_model()

    def insert_ifc_gp_ele(self, *descriptor, **kwargs):
        """ insert interfaces and gaps into seq_model and eles into ele_model

        Args:
            descriptor: a tuple of additions for the sequential, element and
                        part tree models
            kwargs: keyword arguments including
                idx: insertion point in the sequential model
                insert: if True, insert the chunk, otherwise replace it
                t: the thickness following a chunk when inserting
        """
        sm = self['seq_model']
        em = self['ele_model']
        pt = self['part_tree']
        seq, elm, e_nodez = descriptor
        if 'idx' in kwargs:
            sm.cur_surface = kwargs['idx']
        idx = sm.cur_surface

        if isinstance(e_nodez, Sequence):
            for node in e_nodez:
                node.parent = pt.root_node
        else:
            e_nodez.parent = pt.root_node

        # distinguish between adding a new chunk, which requires splitting a
        #  gap in two, and replacing a node, which uses the existing gaps.
        ins_prev_gap = False
        if 'insert' in kwargs:
            t_after = kwargs['t'] if 't' in kwargs else 0.
            if sm.get_num_surfaces() == 2:
                # only object space gap, add image space gap following this
                gap_label = "Image space"
                gap_tag = '#image'
                ins_prev_gap = False
            else:
                # we have both object and image space gaps; retain the image
                # space gap by splitting and inserting the new gap before the
                # inserted chunk, unless we're inserting before idx=1.
                gap_label = None
                gap_tag = ''
                if idx > 0:
                    ins_prev_gap = True

            if ins_prev_gap:
                t_air, sm.gaps[idx].thi = sm.gaps[idx].thi, t_after
                z_dir = sm.z_dir[idx]
            else:
                t_air = t_after
                z_dir = seq[-1][mc.Zdir]
            g, ag, ag_node = ele.create_air_gap(t=t_air, label=gap_label,
                                                z_dir=z_dir, tag=gap_tag)
            if not ins_prev_gap:
                seq[-1][mc.Gap] = g
            elm.append(ag)
            ag_node.parent = pt.root_node
        else:
            # replacing an existing node. need to hook new chunk final
            # interface to the existing gap and following (air gap) element
            g = sm.gaps[sm.cur_surface+1]
            seq[-1][mc.Gap] = g
            ag, ag_node = pt.parent_object(g, '#airgap')

        # insert the new seq into the seq_model
        for sg in seq:
            if ins_prev_gap:
                gap, g = g, sg[mc.Gap]
            else:
                gap = sg[mc.Gap]
            sm.insert(sg[mc.Intfc], gap, z_dir=sg[mc.Zdir], prev=ins_prev_gap)
        sm.update_model()
        self.update_optical_properties()

        # add new elements into the ele_model and 
        #  re-sync them with the seq_model
        for e in elm:
            em.add_element(e)
        em.sync_to_seq(sm)

        # re-sort the part_tree to incorporate the new seq
        pt.update_model()
        # re-sort the ele_model by position on Z axis
        em.sequence_elements()

    def remove_ifc_gp_ele(self, *descriptor, **kwargs):
        """ remove interfaces and gaps from seq_model and eles from ele_model
        """
        sm = self['seq_model']
        em = self['ele_model']
        pt = self['part_tree']
        seq, elm, e_nodez = descriptor
        sg = seq[0]
        idx = sm.ifcs.index(sg[mc.Intfc])

        # verify that the sequences match
        seq_match = True
        for i, sg in enumerate(seq):
            if sg[0] is not sm.ifcs[idx+i]:
                seq_match = False
                break

        if seq_match:
            # remove interfaces in reverse
            for i in range(idx+len(seq)-1, idx-1, -1):
                sm.remove(i)
            sm.update_model()

        for e in elm:
            em.remove_element(e)
        em.sync_to_seq(sm)

        if isinstance(e_nodez, Sequence):
            for node in e_nodez:
                node.parent = None
        else:
            e_nodez.parent = None
        pt.update_model()
        # re-sort the ele_model by position on Z axis
        em.sequence_elements()

    def remove_node(self, e_node):
        # remove interfaces from seq_model
        self.seq_model.remove_node(e_node)
        # remove elements from ele_model
        self.ele_model.remove_node(e_node)
        # unhook node
        e_node.parent = None
