#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Manager class for a sequential optical model

.. codeauthor: Michael J. Hayford
"""

import itertools
import logging

from typing import Optional

from anytree import Node

import rayoptics.optical.model_constants as mc
from rayoptics.typing import SeqPath, Z_DIR
from rayoptics.coord_geometry_types import Tfm3d

from rayoptics.elem import surface
from . import gap
from rayoptics.seq.interface import Interface
from . import medium
from rayoptics.raytr import raytrace as rt
from rayoptics.raytr import trace as trace
from rayoptics.raytr import waveabr
from rayoptics.elem import transform as trns
from opticalglass import glassfactory as gfact
from opticalglass import glasserror as ge
from opticalglass import modelglass as mg
from opticalglass import opticalmedium as om

import numpy as np
from math import copysign, sqrt
from rayoptics.util.misc_math import isanumber

logger = logging.getLogger(__name__)

class SequentialModel:
    """ Manager class for a sequential optical model

    A sequential optical model is a sequence of surfaces and gaps.

    The sequential model has this structure
    ::

        IfcObj  Ifc1  Ifc2  Ifc3 ... Ifci-1   IfcImg
             \  /  \  /  \  /             \   /
             GObj   G1    G2              Gi-1

    where

        - Ifc is a :class:`~rayoptics.seq.interface.Interface` instance
        - G   is a :class:`~rayoptics.seq.gap.Gap` instance

    There are N interfaces and N-1 gaps. The initial configuration has an
    object and image Surface and an object gap.

    The Interface API supports implementation of an optical action, such as
    refraction, reflection, scatter, diffraction, etc. The Interface may be
    realized as a physical profile separating the adjacent gaps or an idealized
    object, such as a thin lens or 2 point HOE.

    The Gap class maintains a simple separation (z translation) and the medium
    filling the gap. More complex coordinate transformations are handled
    through the Interface API.

    Attributes:
        opt_model: parent optical model
        ifcs: list of :class:`~rayoptics.seq.interface.Interface`
        gaps: list of :class:`~rayoptics.seq.gap.Gap`
        lcl_tfrms: forward transform, interface to interface
        rndx: a list with refractive indices for all **wvls**
        z_dir: -1 if gap follows an odd number of reflections, otherwise +1
        gbl_tfrms: global coordinates of each interface wrt the 1st interface
        stop_surface (int): index of stop interface
        cur_surface (int): insertion index for next interface
    """

    def __init__(self, opt_model, do_init=True, **kwargs):
        self.opt_model = opt_model

        self.ifcs: list[Interface] = []
        self.gaps: list[gap.Gap] = []
        self.z_dir: list[Z_DIR] = []

        self.do_apertures: bool = True
        
        self.stop_surface: Optional[int] = None
        self.cur_surface: int = 0

        # derived attributes
        self.gbl_tfrms: list[Tfm3d] = []
        self.lcl_tfrms: list[Tfm3d] = []

        self.seq_def = SequentialStr(self)

        # data for a wavelength vs index vs gap data arrays
        self.wvlns: list[float] = []  # sampling wavelengths in nm
        self.rndx: list[list[float]] = []  # refractive index vs wv and gap

        if do_init:
            self._initialize_arrays()

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['opt_model']
        del attrs['gbl_tfrms']
        del attrs['lcl_tfrms']
        del attrs['wvlns']
        del attrs['rndx']
        del attrs['seq_def']
        return attrs

    def _initialize_arrays(self):
        """ initialize object and image interfaces and intervening gap """
        # add object interface
        self.ifcs.append(surface.Surface('Obj', interact_mode='dummy'))

        tfrm = np.identity(3), np.array([0., 0., 0.])
        self.gbl_tfrms.append(tfrm)
        self.lcl_tfrms.append(tfrm)

        # add object gap
        self.gaps.append(gap.Gap())
        self.z_dir.append(1)
        self.rndx.append([1.0])

        # interfaces are inserted after cur_surface
        self.cur_surface = 0

        # add image interface
        self.ifcs.append(surface.Surface('Img', interact_mode='dummy'))
        self.gbl_tfrms.append(tfrm)
        self.lcl_tfrms.append(tfrm)

        self.seq_def.seq_str = "dad"

    def reset(self):
        self.__init__(self.opt_model)

    def get_num_surfaces(self) -> int:
        return len(self.ifcs)

    def path(self, wl=None, start=None, stop=None, step=1):
        """ returns an iterable path tuple for a range in the sequential model

        Args:
            wl: wavelength in nm for path, defaults to central wavelength
            start: start of range
            stop: first value beyond the end of the range
            step: increment or stride of range

        Returns:
            (**ifcs, gaps, lcl_tfrms, rndx, z_dir**)
        """
        if wl is None:
            wl = self.central_wavelength()

        if step < 0:
            gap_start = start - 1 if start is not None else start
        else:
            gap_start = start

        wl_idx = self.index_for_wavelength(wl)
        try:
            rndx = [n[wl_idx] for n in self.rndx[start:stop:step]]
        except IndexError:
            self.wvlns = self.opt_model['osp']['wvls'].wavelengths
            self.rndx = self.calc_ref_indices_for_spectrum(self.wvlns)
            rndx = [n[wl_idx] for n in self.rndx[start:stop:step]]

        path = itertools.zip_longest(self.ifcs[start:stop:step],
                                     self.gaps[gap_start:stop:step],
                                     self.lcl_tfrms[start:stop:step],
                                     rndx,
                                     self.z_dir[start:stop:step])
        return path

    def reverse_path(self, start: int, stop=None, step=-1, wl=None):
        """ returns an iterable path tuple for a range in the sequential model
    
        Args:
            wl: wavelength in nm for path, defaults to central wavelength
            start: start of range
            stop: first value beyond the end of the range
            step: increment or stride of range
    
        Returns:
            (**ifcs, gaps, lcl_tfrms, rndx, z_dir**)
        """
        if wl is None:
            wl = self.central_wavelength()
    
        if step < 0:
            gap_start = start - 1
            rndx_start = start - 1
        else:
            gap_start = start
    
        tfrms = self.compute_local_transforms(step=-1)
        wl_idx = self.index_for_wavelength(wl)
        rndx = [n[wl_idx] for n in self.rndx[rndx_start:stop:step]]
        z_dir = [-z_dir for z_dir in self.z_dir[start:stop:step]]
        path = itertools.zip_longest(self.ifcs[start:stop:step],
                                     self.gaps[gap_start:stop:step],
                                     tfrms[-(start+1)::+1],
                                     rndx,
                                     z_dir)
        return path

    def seq_str(self):
        """ return a character encoding of `ifcs` and `gaps` """
        return self.seq_def.seq_str

    def calc_ref_indices_for_spectrum(self, wvls):
        """ returns a list with refractive indices for all **wvls**

        Args:
            wvls: list of wavelengths in nm
        """
        indices = []
        for g in self.gaps:
            ri = []
            mat = g.medium
            for w in wvls:
                rndx = mat.rindex(w)
                ri.append(rndx)
            indices.append(ri)

        return indices

    def central_wavelength(self):
        """ returns the central wavelength in nm of the model's `WvlSpec` """
        spectral_region = self.opt_model['optical_spec'].spectral_region
        return spectral_region.central_wvl

    def index_for_wavelength(self, wvl):
        """ returns index into rndx array for wavelength `wvl` in nm """
        spectral_region = self.opt_model['optical_spec'].spectral_region
        self.wvlns = spectral_region.wavelengths
        return self.wvlns.index(wvl)

    def central_rndx(self, i):
        """ returns the central refractive index of the model's `WvlSpec` """
        spectral_region = self.opt_model['optical_spec'].spectral_region
        central_wvl = spectral_region.reference_wvl
        return self.rndx[i][central_wvl]

    def get_surface_and_gap(self, srf: Optional[int]=None):
        if srf is None:
            srf = self.cur_surface
        s = self.ifcs[srf]
        if srf == len(self.gaps):
            g = None
        else:
            g = self.gaps[srf]
        return s, g

    def set_cur_surface(self, s: int):
        self.cur_surface = s

    def set_stop(self, cur_idx: Optional[int] = None) -> Optional[int]:
        """ sets the stop surface to the current surface """
        if cur_idx is None:
            cur_idx = self.cur_surface
        if cur_idx is not None and len(self.ifcs)>2:
            self.stop_surface = cur_idx if cur_idx>0 else 1
        else:
            self.stop_surface = None
        return self.stop_surface

    def insert(self, ifc, gap, z_dir=1, idx=None):
        """ insert ifc and gap *after* cur_surface in seq_model lists """

        if idx is None:
            num_ifcs = len(self.ifcs)
            if self.stop_surface is not None:
                if num_ifcs > 2:
                    if self.stop_surface > self.cur_surface and \
                    self.stop_surface < num_ifcs - 2:
                        self.stop_surface += 1
            
            idx = 0 if num_ifcs < 1 else self.cur_surface+1

        self.cur_surface = idx

        self.ifcs.insert(idx, ifc)
        self.seq_def.insert_token(idx, ifc)
        if gap is not None:
            self.gaps.insert(idx, gap)
            self.seq_def.insert_token(idx, gap)
            z_dir = 1 if z_dir is None else z_dir
            new_z_dir = z_dir*self.z_dir[idx-1] if idx > 1 else z_dir
            self.z_dir.insert(idx, new_z_dir)
        else:
            gap = self.gaps[idx]

        tfrm: Tfm3d = np.identity(3), np.array([0., 0., 0.])
        self.gbl_tfrms.insert(idx, tfrm)
        self.lcl_tfrms.insert(idx, tfrm)

        wvls = self.opt_model.optical_spec.spectral_region.wavelengths
        rindex = [gap.medium.rindex(w) for w in wvls]
        self.rndx.insert(idx, rindex)

        if ifc.interact_mode == 'reflect':
            self.update_reflections(start=idx)

    def remove(self, *args, prev=False):
        """Remove surf and gap at cur_surface or an input index argument.

        If no arguments are provided, the cur_surface is deleted. If a single 
        index is provided, that interface/gap pair is deleted, otherwise a 
        range of interface/gap pairs is deleted.

        To avoid invalid sequence states, both an interface and a gap must be
        removed at the same time. The ``prev`` argument, if True, removes the
        gap preceding the interface. The default behavior is to remove the
        following gap.

        Neither the object not the image interfaces may be removed.

        """
        num_args = len(args)
        if num_args == 0:
            idx_1 = idx_k = self.cur_surface
        elif num_args == 1:
            idx_1 = idx_k = args[0]
        else:
            idx_1 = args[0]
            idx_k = args[1]

        num_ifcs = len(self.ifcs)
        # don't allow object or image interfaces to be removed
        if prev:
            if idx_1 == 0 or idx_1 == 1 or idx_k == num_ifcs:
                raise IndexError
        else:
            if idx_1 == 0 or idx_k == -1 or idx_k == num_ifcs:
                raise IndexError

        self.seq_def.remove_tokens(idx_1, idx_k)

        # loop backward over the range to delete the intended objects
        for idx in range(idx_k, idx_1-1, -1):
            if self.ifcs[idx].interact_mode == 'reflect':
                self.update_reflections(start=idx)

            # decrement stop surface as needed
            if self.stop_surface is not None:
                if len(self.ifcs) > 2:
                    if self.stop_surface > idx and self.stop_surface > 1:
                        self.stop_surface -= 1
        
            # decrement cur_surface as needed
            if self.cur_surface is not None:
                if len(self.ifcs) > 2:
                    if self.cur_surface > idx and self.cur_surface > 0:
                        self.cur_surface -= 1

            # interface related attribute lists
            del self.ifcs[idx]
            del self.gbl_tfrms[idx]
            del self.lcl_tfrms[idx]

            # gap node and related attribute lists
            idx = idx-1 if prev else idx

            del self.gaps[idx]
            del self.z_dir[idx]
            del self.rndx[idx]

        if self.stop_surface is not None:
            if len(self.ifcs)-1 == self.stop_surface:
                self.stop_surface -= 1

    def replace_node_with_seq(self, e_node, seq, **kwargs):
        """ Replace a sub-sequence of e_node with seq. """
        if e_node is None:
            idx_1 = idx_k = kwargs.get('idx', self.cur_surface)
        else:
            pt = self.opt_model['part_tree']
            ifcs = [n.id for n in pt.nodes_with_tag(tag='#ifc', root=e_node)]
            if len(ifcs) > 0:
                idx_1 = self.ifcs.index(ifcs[0])
                idx_k = self.ifcs.index(ifcs[-1])
            else:
                idx_1 = idx_k = self.cur_surface

        self.remove_node(idx_1, idx_k, merge=False)

        # add seq into self.
        tfrm = np.identity(3), np.array([0., 0., 0.])
        for idx, sg in enumerate(seq, start=idx_1):
            self.ifcs.insert(idx, sg[mc.Intfc])
            self.lcl_tfrms.insert(idx, sg[mc.Tfrm])
            self.gbl_tfrms.insert(idx, tfrm)
            if sg[mc.Gap] is not None:
                self.gaps.insert(idx, sg[mc.Gap])
                self.z_dir.insert(idx, sg[mc.Zdir])
                self.rndx.insert(idx, sg[mc.Indx])

        # figure out where the stop belongs
        if self.stop_surface is not None:
            idx_stop = self.stop_surface
            if idx_stop <= idx_1:
                pass
            elif idx_stop <= idx_k:
                # stop was interior to replaced node. set to float because
                # entrance pupil should be well defined.
                self.stop_surface = None
            else:
                idx_delta = idx_stop - idx_k
                idx_stop_new = idx_1 + idx_delta
                self.stop_surface = idx_stop_new

        # handle inserted reflecting interfaces
        self.scan_for_reflections(start=idx_1)

    def remove_node(self, idx_1, idx_k, merge: bool = True, **kwargs):
        """ Remove a range of ifcs/gaps by indices. """
        
        if (idx_stop:=self.stop_surface) is not None:
            if idx_stop > idx_1 and idx_stop <= idx_k:
                idx_stop = idx_1
            elif idx_stop > idx_k:
                idx_stop -= idx_k - idx_1 + 1
            
            # make sure the stop and image surfs are separate;
            #  move stop to previous surf if in conflict
            if idx_stop == 0:
                idx_stop = 1
            elif idx_stop == len(self.ifcs)-1:
                idx_stop = len(self.ifcs)-2

            self.stop_surface = idx_stop

        idx_0 = idx_1-1 if idx_1 > 0 else 0 
        if merge:
            self.gaps[idx_0].thi += self.gaps[idx_k].thi
    
        # delete in reverse 
        for idx in range(idx_k, idx_0, -1):
            del self.ifcs[idx]
            del self.lcl_tfrms[idx]
            del self.gbl_tfrms[idx]
            if idx != idx_k or merge:
                del self.gaps[idx]
                del self.z_dir[idx]
                del self.rndx[idx]

        incr = 1 if merge else 0
        self.seq_def.remove_tokens(idx_1, idx_k, inc_k=incr)
    
    def scan_for_reflections(self, start=0):
        """ Rectify any inconsistent z_dir setting due to insertions. """
        b4_idx = start if start == 0 else start-1
        z_dir_before = self.z_dir[b4_idx]
        reflected = z_dir_before
        
        seq = itertools.zip_longest(self.ifcs[start:], 
                                    self.gaps[start:],
                                    self.z_dir[start:])

        for i, sgz in enumerate(seq, start=start):
            ifc, g, z_dir_after = sgz
            if ifc.interact_mode == 'reflect' or reflected != z_dir_after:
                if i != start:
                    ifc.update_following_reflection()
                if g:
                    g.apply_scale_factor(-1)
                    z_dir_after = -z_dir_after
                # update the reflected state (-1 if odd # of reflections)
                if ifc.interact_mode == 'reflect':
                    reflected = -reflected

            if g is not None:
                z_dir_before = z_dir_after
                self.z_dir[i] = z_dir_after

    def add_surface(self, surf_data, **kwargs):
        """ add a surface where `surf_data` is a list that contains:

        [curvature, thickness, refractive_index, v-number, semi-diameter]

        The `curvature` entry is interpreted as radius if `radius_mode` is **True**

        The `thickness` is the signed thickness

        The `refractive_index, v-number` entry can have several forms:

            - **refractive_index, v-number** (numeric)
            - **refractive_index** only -> constant index model
            - **glass_name, catalog_name** as 1 or 2 strings
            - an instance with a :meth:`~opticalglass.opticalmedium.OpticalMedium.rindex` attribute
            - **air**, str -> :class:`~opticalglass.opticalmedium.Air`
            - blank -> defaults to :class:`~opticalglass.opticalmedium.Air`
            - **'REFL'** -> set interact_mode to 'reflect'

        If no catalog_name is provided, the gfact._cat_names list is used.

        The `semi-diameter` entry is optional. It may also be entered using 
        the `sd` keyword argument.

        """
        radius_mode = self.opt_model.radius_mode
        mat = None
        if len(surf_data) > 2:
            if not isanumber(surf_data[2]):
                if (isinstance(surf_data[2], str) 
                    and 
                    surf_data[2].upper() == 'REFL'):
                    mat = self.gaps[self.cur_surface].medium
        s, g, z_dir, rn, tfrm = create_surface_and_gap(surf_data,
                                                       prev_medium=mat,
                                                       radius_mode=radius_mode,
                                                       **kwargs)
        self.insert(s, g, z_dir=z_dir)

    def add_coord_break(self, thi: float, decenter_data=None, **kwargs):
        """ add a phantom dummy surface to alter the following coordinate systems

        Args:
            thi:  the gap thickness following the break
            decenter_data: the coordinate transform to apply
        """
        s = surface.Surface(interact_mode='phantom', decenter=decenter_data)
        mat = self.gaps[self.cur_surface].medium
        g = gap.Gap(thi, mat)
        z_dir = self.z_dir[self.cur_surface]

        self.insert(s, g, z_dir=z_dir)

    def sync_to_restore(self, opt_model):
        self.opt_model = opt_model
        if hasattr(self, 'optical_spec'):
            opt_model.optical_spec = self.optical_spec
            delattr(self, 'optical_spec')
        init_z_dir = False
        if not hasattr(self, 'z_dir'):
            self.z_dir = []
            init_z_dir = True
            z_dir_work = 1
        for sg in itertools.zip_longest(self.ifcs, self.gaps):
            ifc, g = sg
            if hasattr(ifc, 'sync_to_restore'):
                ifc.sync_to_restore(opt_model)
            if g:
                if hasattr(g, 'sync_to_restore'):
                    g.sync_to_restore(self)
                if init_z_dir:
                    if ifc.interact_mode == 'reflect':
                        z_dir_work = -z_dir_work
                    self.z_dir.append(z_dir_work)

        self.ifcs[0].interact_mode = 'dummy'
        self.ifcs[-1].interact_mode = 'dummy'
                
        if not hasattr(self, 'seq_def'):
            self.seq_def = SequentialStr(self)

        if not hasattr(self, 'do_apertures'):
            self.do_apertures = True

    def update_model(self, **kwargs):
        # delta n across each surface interface must be set to some
        #  reasonable default value. use the index at the central wavelength
        spectral_region = self.opt_model['optical_spec'].spectral_region
        ref_wl = spectral_region.reference_wvl

        self.wvlns = spectral_region.wavelengths
        self.rndx = self.calc_ref_indices_for_spectrum(self.wvlns)

        num_ifcs = len(self.ifcs)
        if self.cur_surface is not None:
            if num_ifcs == 2:
                self.cur_surface = 0
            elif self.cur_surface >= num_ifcs:
                self.cur_surface = num_ifcs - 1
        else:  # if None set cur_surface to insert before image surface
            self.cur_surface = num_ifcs - 2
        
        start = kwargs.get('start', 0)
        b4_idx = start if start == 0 else start-1
        n_before = self.rndx[b4_idx][ref_wl]
        z_dir_before = self.z_dir[b4_idx]

        seq = itertools.zip_longest(self.ifcs[start:], self.gaps[start:])

        for i, sg in enumerate(seq, start=start):
            ifc, g = sg
            z_dir_after = int(copysign(1, z_dir_before))
            if ifc.interact_mode == 'reflect':
                z_dir_after = -z_dir_after

            # leave rndx data unsigned, track change of sign using z_dir
            if g is not None:
                n_after = self.rndx[i][ref_wl]
                if z_dir_after < 0:
                    n_after = -n_after
                ifc.delta_n = n_after - n_before
                n_before = n_after
    
                z_dir_before = z_dir_after
                self.z_dir[i] = z_dir_after

            # call update() on the surface interface
            ifc.update()

        self.gbl_tfrms = self.compute_global_coords()
        self.lcl_tfrms = self.compute_local_transforms()

        self.seq_def.update()

    def update_optical_properties(self, **kwargs):
        if self.do_apertures:
            if len(self.ifcs) > 2:
                self.set_clear_apertures()

    def apply_scale_factor(self, scale_factor):
        """ Apply the `scale_factor` to entire seq_model. """
        self.apply_scale_factor_over(scale_factor)

    def apply_scale_factor_over(self, scale_factor, *surfs):
        """ Apply the `scale_factor` to the `surfs` arg. 
        
        - If `surfs` isn't present, the `scale_factor` is applied to all interfaces and gaps.
        - If `surfs` contains a single value, it is applied to that interface and gap.
        - If `surfs` contains 2 values it is considered an interface range and the `scale_factor` is applied to the interface range and the gaps contained between the outer interfaces.

        """
        if len(surfs) == 0:
            surfs = 0, len(self.ifcs)

        if len(surfs) == 1:
            idx = surfs[0]
            logger.debug(f"{scale_factor=}, {idx=}")
            self.ifcs[idx].apply_scale_factor(scale_factor)
            if idx < len(self.gaps):
                self.gaps[idx].apply_scale_factor(scale_factor)

        elif len(surfs) == 2:
            idx1, idx2 = surfs
            logger.debug(f"{scale_factor=}, {idx1=}, {idx2=}")
            for i in range(idx1, idx2+1):
                try:
                    self.ifcs[i].apply_scale_factor(scale_factor)
                    logger.debug(f"{i}: ifc")
                    if i < idx2:
                        self.gaps[i].apply_scale_factor(scale_factor)
                        logger.debug(f"{i}: gap")
                except IndexError as ie:
                    break

        self.gbl_tfrms = self.compute_global_coords()
        self.lcl_tfrms = self.compute_local_transforms()

    def flip(self, idx1: int, idx2: int) -> None:
        """Flip interfaces and gaps from *idx1* thru *idx2*."""
        def partial_reverse(list_, idx1: int, idx2: int):
            for i in range(0, int((idx2 - idx1)/2)+1):
                a, b = idx1+i, idx2-i
                if a < b:
                    (list_[a], list_[b]) = (list_[b], list_[a])

        if idx2 < idx1:
            idx1, idx2 = idx2, idx1
        partial_reverse(self.ifcs, idx1, idx2)
        partial_reverse(self.gaps, idx1, idx2-1)

        for ifc in self.ifcs[idx1:idx2+1]:
            ifc.flip()

        if self.stop_surface is not None:
            # if the stop surface is in the flip range, flip it too
            stop_idx = self.stop_surface
            if stop_idx >= idx1 and stop_idx <= idx2:
                self.stop_surface = idx2 - (stop_idx - idx1)

        self.update_model()

    def set_from_specsheet(self, specsheet):
        if 'parax_data' not in self.opt_model['analysis_results']:
            return
        if self.opt_model['analysis_results']['parax_data'] is None:
            return
        if specsheet.imager_defined():
            fod = self.opt_model['analysis_results']['parax_data'].fod
            imager = specsheet.imager
            f_old = fod.efl
            f_new = imager.f
            scale_factor = f_new/f_old
            if scale_factor < 1.0-1e-5 or scale_factor > 1.0+1e-5:
                self.apply_scale_factor(scale_factor)

            if specsheet.conjugate_type == 'finite':
                self.gaps[0].thi = -(imager.s + scale_factor*fod.pp1)
                self.gaps[-1].thi = imager.sp - scale_factor*fod.ppk
            elif specsheet.conjugate_type == 'infinite':
                self.gaps[-1].thi = imager.sp - scale_factor*fod.ppk

    def insert_surface_and_gap(self):
        s = surface.Surface()
        g = gap.Gap()
        self.insert(s, g)
        return s, g

    def update_reflections(self, start):
        """ update interfaces and gaps following insertion of a mirror """

        for i, sg in enumerate(self.path(start=start), start=start):
            ifc, g, lcl_tfrm, rndx, z_dir = sg
            if i > start:
                ifc.update_following_reflection()
                if g:
                    g.apply_scale_factor(-1)
                    self.z_dir[i] = -z_dir

    def get_rndx_and_imode(self):
        """ get list of signed refractive index and interact mode for sequence. """
        central_wvl = self.opt_model['osp']['wvls'].reference_wvl
        rndx_and_imode = []
        for i in range(len(self.rndx)):
            rndx = self.rndx[i][central_wvl]
            n = rndx if self.z_dir[i] > 0 else -rndx
            imode = self.ifcs[i].interact_mode
            rndx_and_imode += [(n, imode)]
        rndx_and_imode += [(n, imode)]
        return rndx_and_imode

    def overall_length(self, os_idx=1, is_idx=-1):
        """ Sum gap thicknesses from `os_idx` to `is_idx` 
        
        The default arguments return the thickness sum between the 1st and last surfaces.

        To include the image surface, is_idx=len(sm.gaps)

        Args:
            os_idx: starting gap index
            is_idx: final gap index

        Returns:
            oal: float, overal length of gap range
        """
        oal = 0
        for g in self.gaps[os_idx:is_idx]:
            oal += g.thi
            # print(f"{oal}    +{g.thi}")
        return oal
    
    def total_track(self):
        """ Total track length, distance from object to image. """
        return self.overall_length(0, len(self.gaps))

    def surface_label_list(self):
        """ list of surface labels or surface number, if no label """
        labels = []
        for i, s in enumerate(self.ifcs):
            if len(s.label) == 0:
                if i == self.stop_surface:
                    labels.append('Stop')
                else:
                    labels.append(str(i))
            else:
                labels.append(s.label)
        return labels

    def list_model(self, path=None):
        cvr = 'r' if self.opt_model.radius_mode else 'c'
        print("              {}            t        medium     mode   zdr"
              "      sd".format(cvr))
        labels = self.surface_label_list()
        path = self.path() if path is None else path
        prev_z_dir = 1
        for i, sg in enumerate(path):
            ifc, gap, _, _, z_dir = sg
            s = self.list_surface_and_gap(ifc, gp=gap)
            if gap is not None:
                s.append(z_dir)
            else:
                s.append(prev_z_dir)
            fmt = "{0:>5s}: {1:12.6f} {2:#12.6g} {3:>9s} {4:>10s} {6:2n}"
            if s[4] is not None:  # if the sd is not None...
                fmt += "  {5:#10.5g}"
            print(fmt.format(labels[i], *s))
            prev_z_dir = z_dir

    def list_model_old(self):
        cvr = 'r' if self.opt_model.radius_mode else 'c'
        print("           {}            t        medium     mode         sd"
              .format(cvr))
        for i, sg in enumerate(self.path()):
            ifc, g, _, _, _ = sg
            s = self.list_surface_and_gap(ifc, gp=g)
            fmt = "{0:2n}: {1:12.6f} {2:#12.6g} {3:>9s} {4.name:>10s}"
            if s[4] is not None:  # if the sd is not None...
                fmt += " {5:#10.5g}"
            print(fmt.format(i, *s))

    def list_gaps(self):
        for i, gp in enumerate(self.gaps):
            print(i, gp)

    def list_surfaces(self):
        for i, s in enumerate(self.ifcs):
            print(i, s)

    def list_surface_and_gap(self, ifc, gp=None):
        """Returns cvr, thi, med, imode, sd for input ifc and gap."""
        cvr = ifc.profile_cv
        if self.opt_model.radius_mode:
            if cvr != 0.0:
                cvr = 1.0/cvr
        sd = ifc.surface_od()
        imode = "" if ifc.interact_mode == 'transmit' else ifc.interact_mode

        if gp is not None:
            thi = gp.thi
            med = gp.medium.name()
        else:
            thi = 0.
            med = ''
        return [cvr, thi, med, imode, sd]

    def list_decenters(self, full=False):
        """List decenter data and gap separations.

        Arguments:
            full: lists all values if True, else only y offset and alpha tilt
        """
        fmt0a = ("              thi    medium/mode          type          x"
                 "          y       alpha      beta       gamma")
        fmt0b = ("              thi    medium/mode          type          y"
                 "       alpha")
        fmt1a = ("{:6s}                {:>10s}  {:>14s} {:#10.5g} {:#10.5g}"
                 " {:#10.5g} {:#10.5g} {:#10.5g}")
        fmt1b = ("{:6s}                {:>10s}  {:>14s} {:#10.5g}"
                 " {:#10.5g}")
        fmt1c = "{:6s}                {:>10s}"
        fmt2 = "{:6s} {:#12.6g}    {:>9s}"

        # print header
        if full:
            print(fmt0a)
        else:
            print(fmt0b)

        for i, sg in enumerate(self.path()):
            ifc, gap, lcl_tfrm, rndx, z_dir = sg
            idx = f"{i:5n}:"
            imode = (ifc.interact_mode if ifc.interact_mode != 'transmit'
                     else "")

            if ifc.decenter is not None:
                d = ifc.decenter
                if full:
                    print(fmt1a.format(idx, imode, d.dtype,
                                       d.dec[0], d.dec[1],
                                       d.euler[0], d.euler[1], d.euler[2]))
                else:
                    print(fmt1b.format(idx, imode, d.dtype,
                                       d.dec[1], d.euler[0]))
                idx = f"{'':5s} "
            elif gap is None:  # final interface, just list interact_mode
                print(fmt1c.format(idx, imode))

            if gap:
                print(fmt2.format(idx, gap.thi, gap.medium.name()))

    def list_sg(self):
        """List decenter data and gap separations. """
        cvrd = 'r' if self.opt_model.radius_mode else 'c'
        fmt0a = ("               {}               mode              type"
                 "          y       alpha")
        fmt0b = ("                       t           medium")
        fmt1 = ("{:>5s}: {:#12.6g}       {:>10s}     {:>14s} {:#10.5g}"
                " {:#10.5g}")
        fmt2 = ("{:>5s}: {:#12.6g}       {:>10s}")
        fmt3 = ("                {:#12.6g}    {:>9s}")

        # print header
        print(fmt0a.format(cvrd))
        print(fmt0b)

        labels = self.surface_label_list()
        for i, sg in enumerate(self.path()):
            ifc, gap, lcl_tfrm, rndx, z_dir = sg
            s = self.list_surface_and_gap(ifc, gap)
            s.append(z_dir)
            cvr, thi, med, imode, sd, z_dir = s
            if ifc.decenter is not None:
                d = ifc.decenter
                print(fmt1.format(labels[i], cvr, imode, d.dtype,
                                  d.dec[1], d.euler[0]))
            else:
                print(fmt2.format(labels[i], cvr, imode))

            if gap:
                print(fmt3.format(gap.thi, gap.medium.name()))

    def list_sg_ele(self, part_tree):
        seq_str = ''
        ele_list = []
        for i, sgz in enumerate(
            itertools.zip_longest(self.ifcs, self.gaps, self.z_dir)):
            s, g, z_dir = sgz

            s_str = s.ifc_token()
            s_parent_node = part_tree.parent_node(s)
            if s_parent_node is not None:
                s_parent = s_parent_node.name
            else:
                s_parent = " "
            seq_str += s_str
            ele_list.append(s_parent)
            if g is not None:
                g_str = 'a' if g.medium.name() == 'air' else 't'
                g_parent_node = part_tree.parent_node((g, z_dir))
                if g_parent_node is not None:
                    g_parent = g_parent_node.name
                else:
                    g_parent = " "
                seq_str += g_str
                ele_list.append(g_parent)
            else:
                g_str = " "
                g_parent = " "
            print(f"{s_str}{g_str}      {s_parent:10s} {g_parent}")
        return seq_str, ele_list

    def listobj_str(self):
        o_str = ""
        stop_idx = self.stop_surface
        for i, sg in enumerate(self.path()):
            ifc, gap, lcl_tfrm, rndx, z_dir = sg
            if stop_idx is not None and i == stop_idx:
                o_str += f'{i} (stop): ' + ifc.listobj_str()
            else:
                o_str += f'{i}: ' + ifc.listobj_str()
            if gap is not None:
                gap_str = gap.listobj_str()
                semicolon_indx = gap_str.find(';')
                o_str += (gap_str[:semicolon_indx] +
                          f" ({int(z_dir):+})" +
                          gap_str[semicolon_indx:] + '\n')
        
        o_str += f'\ndo apertures: {self.do_apertures}'
        return o_str

    def trace_fan(self, fct, fi, xy, num_rays=21, **kwargs):
        """ xy determines whether x (=0) or y (=1) fan """
        osp = self.opt_model.optical_spec
        fld = osp.field_of_view.fields[fi]
        wvl = self.central_wavelength()
        foc = osp.defocus.get_focus()

        rs_pkg, cr_pkg = trace.setup_pupil_coords(self.opt_model,
                                                  fld, wvl, foc)
        fld.chief_ray = cr_pkg
        fld.ref_sphere = rs_pkg

        # Use the central wavelength reference image point for the wavefront error calculations
        ref_img_pt = rs_pkg[0]

        wvls = osp.spectral_region
        fans_x = []
        fans_y = []
        fan_start = np.array([0., 0.])
        fan_stop = np.array([0., 0.])
        fan_start[xy] = -1.0
        fan_stop[xy] = 1.0
        fan_def = [fan_start, fan_stop, num_rays]
        max_rho_val = 0.0
        max_y_val = 0.0
        rc = []
        for wi, wvl in enumerate(wvls.wavelengths):
            rc.append(wvls.render_colors[wi])

            rs_pkg, cr_pkg = trace.setup_pupil_coords(self.opt_model,
                                                      fld, wvl, foc,
                                                      image_pt=ref_img_pt)
            fld.chief_ray = cr_pkg
            fld.ref_sphere = rs_pkg
            fan = trace.trace_fan(self.opt_model, fan_def, fld, wvl, foc,
                                  img_filter=lambda p, ray_pkg:
                                  fct(p, xy, ray_pkg, fld, wvl, foc), **kwargs)
            f_x = []
            f_y = []
            for p, y_val in fan:
                f_x.append(p[xy])
                f_y.append(y_val)
                if abs(p[xy]) > max_rho_val:
                    max_rho_val = abs(p[xy])
                if abs(y_val) > max_y_val:
                    max_y_val = abs(y_val)
            fans_x.append(f_x)
            fans_y.append(f_y)
        fans_x = np.array(fans_x)
        fans_y = np.array(fans_y)
        return fans_x, fans_y, (max_rho_val, max_y_val), rc

    def trace_grid(self, fct, fi, wl=None, num_rays=21, form='grid',
                   append_if_none=True, **kwargs):
        """ fct is applied to the raw grid and returned as a grid  """
        osp = self.opt_model.optical_spec
        wvls = osp.spectral_region
        wvl = self.central_wavelength()
        wv_list = wvls.wavelengths if wl is None else [wvl]
        fld = osp.field_of_view.fields[fi]
        foc = osp.defocus.get_focus()

        rs_pkg, cr_pkg = trace.setup_pupil_coords(self.opt_model,
                                                  fld, wvl, foc)
        fld.chief_ray = cr_pkg
        fld.ref_sphere = rs_pkg

        grids = []
        grid_start = np.array([-1., -1.])
        grid_stop = np.array([1., 1.])
        grid_def = [grid_start, grid_stop, num_rays]
        for wi, wvl in enumerate(wv_list):
            grid = trace.trace_grid(self.opt_model, grid_def, fld, wvl, foc,
                                    form=form, append_if_none=append_if_none,
                                    img_filter=lambda p, ray_pkg:
                                    fct(p, wi, ray_pkg, fld, wvl, foc),
                                    **kwargs)
            grids.append(grid)
        rc = wvls.render_colors
        return grids, rc

    def trace_wavefront(self, fld, wvl, foc, num_rays=32):

        def wave(p, ray_pkg, fld, wvl, foc):
            x = p[0]
            y = p[1]
            if ray_pkg is not None:
                fod = self.opt_model['analysis_results']['parax_data'].fod
                opd = waveabr.wave_abr_full_calc(fod, fld, wvl, foc, ray_pkg,
                                                 fld.chief_ray,
                                                 fld.ref_sphere)
                opd = opd/self.opt_model.nm_to_sys_units(wvl)
            else:
                opd = 0.0
            return np.array([x, y, opd])

        rs_pkg, cr_pkg = trace.setup_pupil_coords(self.opt_model,
                                                  fld, wvl, foc)
        fld.chief_ray = cr_pkg
        fld.ref_sphere = rs_pkg

        grid_start = np.array([-1., -1.])
        grid_stop = np.array([1., 1.])
        grid_def = (grid_start, grid_stop, num_rays)

        grid = trace.trace_grid(self.opt_model, grid_def, fld, wvl, foc,
                                img_filter=lambda p, ray_pkg:
                                wave(p, ray_pkg, fld, wvl, foc), form='grid')
        return grid

    def set_clear_apertures_paraxial(self):
        ax_ray, pr_ray, _ = self.opt_model['analysis_results']['parax_data']
        for i, ifc in enumerate(self.ifcs):
            sd = abs(ax_ray[i][0]) + abs(pr_ray[i][0])
            ifc.set_max_aperture(sd)

    def set_clear_apertures(self, avoid_list: Optional[list[int]]=None):
        """ set each surface aperture using edge rays at each field

        Args:
        avoid_list: list of surfaces to skip when setting apertures.

        The avoid_list idea and implementation was contributed by Quentin Bécar
        """
        if avoid_list is None:
            avoid_list = []

        rayset = trace.trace_boundary_rays(self.opt_model, 
                                           use_named_tuples=True)

        for i, s in enumerate(self.ifcs):
            if i not in avoid_list:
                max_ap = -1.0e+10
                update = True
                for f in rayset:
                    for p in f:
                        ray = p.ray
                        if len(ray) > i:
                            ap = sqrt(ray[i].p[0]**2 + ray[i].p[1]**2)
                            if ap > max_ap:
                                max_ap = ap
                        else:  # ray failed before this interface, don't update
                            update = False
                if update:
                    s.set_max_aperture(max_ap)

    def trace(self, pt0, dir0, wvl, **kwargs):
        return rt.trace(self, pt0, dir0, wvl, **kwargs)

    def compute_global_coords(self, glo=1, origin=None, **kwargs):
        """ Return global surface coordinates (rot, t) wrt surface `glo`. 
        
        If origin isn't None, it should be a tuple (r, t) being the transform
        from the desired global origin to the specified global surface.
        """
        return trns.compute_global_coords(self, glo, origin)

    def compute_local_transforms(self, seq=None, step=1):
        """ Return forward surface coordinates (r.T, t) for each interface. """
        return trns.compute_local_transforms(self, seq, step)

    def list_lcl_tfrms(self, *args):
        trns.list_tfrms(self.lcl_tfrms, *args)

    def list_gbl_tfrms(self, *args):
        trns.list_tfrms(self.gbl_tfrms, *args)

    def list_tfrms(self, tfrms, sel: str='r+t', *args):
        return trns.list_tfrms(tfrms, sel, *args)

    def filter_out_phantoms(self, items: list) -> list:
        """ given items, filter out those corresponding to a phantom interface.

        len(items) == len(self.ifcs)
        """
        return [gt for ifc, gt in zip(self.ifcs, items) 
                    if not ifc.interact_mode == 'phantom']

    def find_matching_ifcs(self):
        rot_tols = dict(atol=1e-14, rtol=1e-8)
        tols = dict(atol=1e-14, rtol=1e-14)
        matches = []
        for i, gi in enumerate(self.gbl_tfrms):
            i1 = i+1
            for j, gj in enumerate(self.gbl_tfrms[i1:], start=i1):
                if (
                        np.allclose(gi[0], gj[0], **rot_tols) and
                        np.allclose(gi[1], gj[1], **tols)
                        ):
                    print(f'coincident surfs: {i} - {j}')
                    matches.append((i, j))
        return matches


def gen_sequence(surf_data_list, **kwargs):
    """ create a sequence iterator from the surf_data_list

    Args:
        surf_data_list: a list of lists containing:
                        [curvature, thickness, refractive_index, v-number]
        **kwargs: keyword arguments

    Returns:
        (**ifcs**, **gaps**, **rndx**, **lcl_tfrms**, **z_dir**)
    """
    ifcs = []
    gaps = []
    rndx = []
    lcl_tfrms = []
    z_dir = []

    for surf_data in surf_data_list:
        s, g, zdir, rn, tfrm = create_surface_and_gap(surf_data, **kwargs)
        ifcs.append(s)
        gaps.append(g)
        rndx.append(rn)
        lcl_tfrms.append(tfrm)
        z_dir.append(zdir)
    ifcs[-1].interact_mode = 'dummy'

    n_before = 1.0
    z_dir_before = 1
    for i, s in enumerate(ifcs):
        z_dir_after = int(copysign(1, z_dir_before))
        n_after = np.copysign(rndx[i], n_before)
        if s.interact_mode == 'reflect':
            n_after = -n_after
            z_dir_after = -z_dir_after

        n_before = n_after
        rndx[i] = n_after
        z_dir_before = z_dir_after
        z_dir[i] = z_dir_after

    seq = itertools.zip_longest(ifcs, gaps[:-2], lcl_tfrms, rndx, z_dir)
    return seq


def create_surface_and_gap(surf_data, radius_mode=False, prev_medium=None,
                           wvl=550.0, **kwargs):
    """ create a surface and gap where `surf_data` is a list that contains:

    [curvature, thickness, refractive_index, v-number, semi-diameter]
    
    The `curvature` entry is interpreted as radius if `radius_mode` is **True**

    The `thickness` is the signed thickness

    The `refractive_index, v-number` entry can have several forms:
        
        - **refractive_index, v-number** (numeric)
        - **refractive_index** only -> constant index model
        - **glass_name, catalog_name** as 1 or 2 strings
        - an instance with a :meth:`~opticalglass.opticalmedium.OpticalMedium.rindex` attribute
        - **air**, str -> :class:`~opticalglass.opticalmedium.Air`
        - blank -> defaults to :class:`~opticalglass.opticalmedium.Air`
        - **'REFL'** -> set interact_mode to 'reflect'

    If no catalog_name is provided, the gfact._cat_names list is used.

    The `semi-diameter` entry is optional. It may also be entered using the 
    `sd` keyword argument.
    """
    s = surface.Surface()

    if radius_mode:
        if surf_data[0] != 0.0:
            s.profile.cv = 1.0/surf_data[0]
        else:
            s.profile.cv = 0.0
    else:
        s.profile.cv = surf_data[0]

    z_dir = 1
    sd_indx = None
    num_inputs = len(surf_data)
    if num_inputs > 2:  # look for medium data, possibly followed by sd
        last_k = 3      # assume medium with 1 input
        if num_inputs >= 5:  # 2 input medium plus sd
            last_k = 4
            sd_indx = 4
        elif num_inputs == 4:  # 2 inputs left
            if type(surf_data[2]) == type(surf_data[3]):
                # if same type, assume 2 medium inputs, no sd
                last_k = 4
            else:
                # different types, assume 1 medium input and sd
                last_k = 3
                sd_indx = 3
        try:
            # Feed the right number of inputs into decode_medium
            if last_k == 3:
                mat = medium.decode_medium(surf_data[2])
            else:
                mat = medium.decode_medium(surf_data[2], surf_data[3])
        except ValueError:
            if isinstance(surf_data[2], str):  # string args
                if surf_data[2].upper() == 'REFL':
                    s.interact_mode = 'reflect'
                    mat = prev_medium
                    z_dir = -1

        if sd_indx:
            s.set_max_aperture(surf_data[sd_indx])

    else:  # only curvature and thickness entered, set material to air
        mat = om.Air()

    if kwargs.get('sd', None) is not None:
        s.set_max_aperture(kwargs.get('sd'))
    thi = surf_data[1]
    g = gap.Gap(thi, mat)
    rndx = mat.rindex(wvl)
    tfrm = np.identity(3), np.array([0., 0., thi])

    return s, g, z_dir, rndx, tfrm


def overall_length(path: SeqPath) -> float:
    """ Sum gap thicknesses for `path` 

    Returns:
        oal: float, overal length of gap range
    """
    oal = 0
    for sg in path:
        g = sg[mc.Gap]
        if g is not None:
            oal += g.thi
        # print(f"{oal}    +{g.thi}")
    return oal


class SequentialStr:
    """Manage the string representation of the sequential path. """
    def __init__(self, seq_model: SequentialModel):
        self.sm = seq_model
        self.seq_str = self.gen_seq_str()
    
    @property
    def seq_str(self) -> str:
        return ''.join(self._seq_str_list)

    @seq_str.setter
    def seq_str(self, s: str):
        self._seq_str_list = list(s)

    def update(self):
        """ Regenerate the seq_str from the current seq_model. """
        self.seq_str = self.gen_seq_str()

    def insert_token(self, idx: int, token: Interface|gap.Gap):
        if isinstance(token, Interface):
            self._seq_str_list.insert(2*idx, token.ifc_token())
        elif isinstance(token, gap.Gap):
            self._seq_str_list.insert(2*idx+1, token.gap_token())

    def remove_tokens(self, idx_1: int, idx_k: int, 
                      inc_1: int = 0, inc_k: int = 0):
        """ Remove tokens from ifcs[idx_1] to ifcs[idx_k] """
        self._seq_str_list[2*idx_1+inc_1:2*idx_k+1+inc_k] = []

    def gen_seq_str(self) -> str:
        """return a character encoding of `ifcs` and `gaps` """

        seq_str = ''
        for sg in itertools.zip_longest(self.sm.ifcs, self.sm.gaps):
            s, g = sg

            s_str = s.ifc_token()
            seq_str += s_str

            if g is not None:
                g_str = g.gap_token()
                seq_str += g_str

        return seq_str
