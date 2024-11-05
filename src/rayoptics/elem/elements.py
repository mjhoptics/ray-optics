#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Module for element modeling

.. Created on Sun Jan 28 16:27:01 2018

.. codeauthor: Michael J. Hayford
"""

from collections import namedtuple
from copy import deepcopy
import itertools
from itertools import zip_longest
from packaging import version

from abc import abstractmethod
from typing import Protocol, ClassVar, List, Dict, Any, runtime_checkable

from math import sqrt
import numpy as np

from anytree import Node  # type: ignore

import rayoptics.optical.model_constants as mc

import rayoptics.util.rgbtable as rgbt
from rayoptics.oprops import thinlens
from rayoptics.elem import parttree
from rayoptics.elem.profiles import SurfaceProfile, Spherical, Conic
from rayoptics.elem.surface import Surface
from rayoptics.elem import transform as trns
from rayoptics.seq.gap import Gap
from rayoptics.seq.medium import decode_medium

from rayoptics.seq.sequential import SequentialModel
from rayoptics.seq.interface import Interface

import rayoptics.gui.appcmds as cmds
from rayoptics.gui.actions import (Action, AttrAction, SagAction, BendAction,
                                   ReplaceGlassAction)
from rayoptics.gui.util import (calc_render_color_for_material, transform_poly)

import opticalglass.glassfactory as gfact  # type: ignore
from opticalglass.modelglass import ModelGlass  # type: ignore

GraphicsHandle = namedtuple('GraphicsHandle', ['polydata', 'tfrm', 'polytype',
                                               'color'], defaults=(None,))
GraphicsHandle.polydata.__doc__ = "poly data in local coordinates"
GraphicsHandle.tfrm.__doc__ = "global transformation for polydata"
GraphicsHandle.polytype.__doc__ = "'polygon' (for filled) or 'polyline'"
GraphicsHandle.color.__doc__ = "RGBA for the polydata or None for default"


""" tuple grouping together graphics rendering data

    Attributes:
        polydata: poly data in local coordinates
        tfrm: global transformation for polydata
        polytype: 'polygon' (for filled) or 'polyline'
"""

rot_around_x = np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]])
rot_around_y = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])


# --- Factory functions
def create_thinlens(power=0., indx=1.5, sd=None, **kwargs):
    tl = thinlens.ThinLens(power=power, ref_index=indx, max_ap=sd, **kwargs)
    tle = ThinElement(ifc=tl)
    tree = tle.tree()

    if 'prx' in kwargs:
        pm, node, type_sel = prx = kwargs['prx']
        dgm_pkg = [pm.get_pt(node)], [[1.0, 'transmit']]
        dgm = prx, dgm_pkg
    else:
        dgm = None

    descriptor = [[tl, None, None, 1, +1]], [tle], tree, dgm
    return descriptor


def create_mirror(c=0.0, r=None, cc=0.0, ec=None,
                  power=None, profile=None, sd=None, **kwargs):
    '''Create a sequence and element for a mirror.

    Args:
        c: vertex curvature
        r: vertex radius of curvature
        cc: conic constant
        ec: = 1 + cc
        power:  optical power of the mirror
        sd:  semi-diameter
        profile: Spherical or Conic type, or a profile instance
    '''
    delta_n = kwargs['delta_n'] if 'delta_n' in kwargs else -2
    if power:
        cv = power/delta_n
    elif r:
        cv = 1.0/r
    else:
        cv = c

    if ec:
        k = ec - 1.0
    else:
        k = cc

    if profile is Spherical:
        prf = Spherical(c=cv)
    elif profile is Conic:
        prf = Conic(c=cv, cc=k)
    elif profile is not None:
        prf = profile
    else:
        if k == 0.0:
            prf = Spherical(c=cv)
        else:
            prf = Conic(c=cv, cc=k)

    sd = sd if sd is not None else 1

    m = Surface(profile=prf, interact_mode='reflect', max_ap=sd,
                delta_n=delta_n, **kwargs)
    ele_kwargs = {'label': kwargs['label']} if 'label' in kwargs else {}
    me = Mirror(ifc=m, sd=sd, **ele_kwargs)

    tree = me.tree()

    if 'prx' in kwargs:
        pm, node, type_sel = prx = kwargs['prx']
        dgm_pkg = [pm.get_pt(node)], [[-1, 'reflect']]
        dgm = prx, dgm_pkg
    else:
        dgm = None

    return [[m, None, None, 1, -1]], [me], tree, dgm


def lens_from_power(power=0., bending=0., th=None, sd=1.,
                    med=None, nom_wvl='d'):
    if med is None:
        med = ModelGlass(1.517, 64.2, '517642')
    else:
        med = decode_medium(med)
    rndx = med.rindex(nom_wvl)

    if th is None:
        th = sd/5
    
    if power == 0:
        cv1 = cv2 = 0
    else:
        abs_bending = abs(bending)
        if abs_bending == 1:
            cv1 = power/(rndx - 1)
            cv2 = 0
        else:
            B = (abs_bending - 1)/(abs_bending + 1)
            a = (rndx - 1)*(th/rndx)*B
            b = 1 - B
            c = -power/(rndx - 1)
            cv1 = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
            cv2 = cv1*B

        if bending < 0:
            cv1, cv2 = -cv2, -cv1

    return cv1, cv2, th, rndx, sd


def _create_lens(power=0., bending=0., th=None, sd=1., med=None, 
                lens=None, **kwargs):
    """ Create a lens element chunk of sm, em, and pt tree 
    
    Args:
        kwargs: keyword arguments including:

            - idx: insertion point in the sequential model
            - t: the thickness following a chunk when inserting
            - lens: tuple of `cv1, cv2, th, glass_name_catalog, sd` where:

                - cv1: front curvature
                - cv2: rear curvature
                - th: lens thickness
                - glass_input: a str, e.g. 'N-BK7, Schott' or index (+V-number)
                - sd: lens semi-diameter

        """
    if med is None:
        mat = ModelGlass(1.517, 64.2, '517642')
    else:
        mat = decode_medium(med)
        
    if lens is None:
        lens = lens_from_power(power=power, bending=bending, th=th, sd=sd,
                               med=mat)
        cv1, cv2, th, rndx, sd = lens
    else:
        cv1, cv2, th, glass, sd = lens
        mat = decode_medium(glass)

    rndx = mat.rindex('d')
    lens = cv1, cv2, th, rndx, sd

    s1 = Surface(profile=Spherical(c=cv1), max_ap=sd, delta_n=(rndx - 1))
    s2 = Surface(profile=Spherical(c=cv2), max_ap=sd, delta_n=(1 - rndx))
    g = Gap(t=th, med=mat)
    le = Element(sg_def=(s1, s2, g), sd=sd)
    tree = le.tree()

    return [[s1, g, None, rndx, 1], [s2, None, None, 1, 1]], [le], tree


def create_lens(power=0., bending=0., th=None, sd=1., med=None, 
                lens=None, **kwargs):
    """ Create a lens element chunk of sm, em, pt tree, and |ybar| entry
    
    Args:
        kwargs: keyword arguments including:

            - idx: insertion point in the sequential model
            - t: the thickness following a chunk when inserting
            - lens: tuple of `cv1, cv2, th, glass_name_catalog, sd` where:

                - cv1: front curvature
                - cv2: rear curvature
                - th: lens thickness
                - glass_input: a str, e.g. 'N-BK7, Schott' or index (+V-number)
                - sd: lens semi-diameter

        """
    descriptor = _create_lens(power, bending, th, sd, med, lens, **kwargs)
    descriptor += (None,)
    return descriptor


def create_lens_from_dgm(prx=None, **kwargs):
    """ Use diagram points to create a lens. 
    
    Adds a |ybar| component to the descriptor tuple.
    dgm = prx, dgm_pkg
    
    prx = parax_model, node_idx, type_sel
    dgm_pkg = node_list, sys_data
    
    sys_data = list([rndx, 'transmit'|'reflect'])
    """
    pm, node, type_sel = prx
    dgm_pkg, lens_from_dgm = pm.lens_from_dgm(node, **kwargs)
    descriptor = _create_lens(lens=lens_from_dgm, **kwargs)
    descriptor += ((prx, dgm_pkg),)
    return descriptor


def achromat(power, Va, Vb):
    """Compute lens powers for a thin doublet achromat, given their V-numbers."""
    power_a = (Va/(Va - Vb))*power
    power_b = (Vb/(Vb - Va))*power
    return power_a, power_b


def create_cemented_doublet(power=0., bending=0., th=None, sd=1.,
                            glasses=('N-BK7,Schott', 'N-F2,Schott'),
                            **kwargs):
    from opticalglass.spectral_lines import get_wavelength  # type: ignore
    from opticalglass import util
    wvls = np.array([get_wavelength(w) for w in ['d', 'F', 'C']])
    gla_a = gfact.create_glass(glasses[0])
    rndx_a = gla_a.calc_rindex(wvls)
    Va, PcDa = util.calc_glass_constants(*rndx_a)
    gla_b = gfact.create_glass(glasses[1])
    rndx_b = gla_b.calc_rindex(wvls)
    Vb, PcDb = util.calc_glass_constants(*rndx_b)

    power_a, power_b = achromat(power, Va, Vb)

    if th is None:
        th = sd/4
    t1 = 3*th/4
    t2 = th/4
    if power_a < 0:
        t1, t2 = t2, t1

    lens_a = lens_from_power(power=power_a, bending=bending, th=t1, sd=sd,
                              med=gla_a)
    cv1, cv2, t1, indx_a, sd = lens_a

    # cv1 = power_a/(rndx_a[0] - 1)
    # delta_cv = -cv1/2
    # cv1 += delta_cv
    # cv2 = delta_cv
    # cv3 = power_b/(1 - rndx_b[0]) + delta_cv
    indx_b = rndx_b[0]
    cv3 = (power_b/(indx_b-1) - cv2)/((t2*cv2*(indx_b-1)/indx_b) - 1)

    s1 = Surface(profile=Spherical(c=cv1), max_ap=sd,
                 delta_n=(rndx_a[0] - 1))
    s2 = Surface(profile=Spherical(c=cv2), max_ap=sd,
                 delta_n=(rndx_b[0] - rndx_a[0]))
    s3 = Surface(profile=Spherical(c=cv3), max_ap=sd,
                 delta_n=(1 - rndx_b[0]))

    g1 = Gap(t=t1, med=gla_a)
    g2 = Gap(t=t2, med=gla_b)

    g_tfrm = np.identity(3), np.array([0., 0., 0.])

    ifc_list = []
    ifc_list.append([0, s1, g1, 1, g_tfrm])
    ifc_list.append([1, s2, g2, 1, g_tfrm])
    ifc_list.append([2, s3, None, 1, g_tfrm])
    ce = CementedElement(ifc_list=ifc_list)
    tree = ce.tree()

    return [[s1, g1, None, rndx_a, 1],
            [s2, g2, None, rndx_b, 1],
            [s3, None, None, 1, 1]], [ce], tree, None


def create_dummy_plane(sd=1., **kwargs):
    s = Surface(interact_mode='dummy', max_ap=sd, **kwargs)
    se = DummyInterface(ifc=s, sd=sd)
    tree = se.tree()

    if 'prx' in kwargs:
        pm, node, type_sel = prx = kwargs['prx']
        dgm_pkg = [pm.get_pt(node)], [[1, 'transmit']]
        dgm = prx, dgm_pkg
    else:
        dgm = None

    descriptor = [[s, None, None, 1, +1]], [se], tree, dgm
    return descriptor


def create_air_gap(t=0., **kwargs):
    g = Gap(t=t)
    ag = AirGap(g=g, **kwargs)
    kwargs.pop('label', None)
    tree = ag.tree(**kwargs)
    return g, ag, tree, None


def create_from_file(filename, **kwargs):
    opm_file = cmds.open_model(filename, post_process_imports=False)
    sm_file = opm_file['seq_model']
    osp_file = opm_file['optical_spec']
    pm_file = opm_file['parax_model']
    em_file = opm_file['ele_model']
    pt_file = opm_file['part_tree']
    ar_file = opm_file['analysis_results']
    if len(pt_file.nodes_with_tag(tag='#element')) == 0:
        parttree.sequence_to_elements(sm_file, em_file, pt_file)

    if 'power' in kwargs:
        desired_power = kwargs['power']
        cur_power = ar_file['parax_data'].fod.power
        # scale_factor is linear, power is 1/linear
        #  so use reciprocal of power to compute scale_factor
        scale_factor = cur_power/desired_power
        opm_file.apply_scale_factor(scale_factor)

    # extract the system definition, minus object and image
    seq = [list(node) for node in sm_file.path(start=1, stop=-1)]
    seq[-1][1] = None
    
    if 'prx' in kwargs:
        dgm = pm_file.match_pupil_and_conj(kwargs['prx'])
    else:
        dgm = None

    # get the top level nodes of the input system, minus object and image
    part_nodes = pt_file.nodes_with_tag(tag='#element#airgap#assembly',
                                        not_tag='#object#image',
                                        node_list=pt_file.root_node.children)
    parts = [part_node.id for part_node in part_nodes]

    if (len(part_nodes) == 1 and '#assembly' in part_nodes[0].tag):
        asm_node = part_nodes[0]
        print("found root assembly node")
    else:
        # create an Assembly from the top level part list
        label = kwargs.get('label', None)
        tfrm = kwargs.get('tfrm', opm_file['seq_model'].gbl_tfrms[1])
        asm = Assembly(parts, idx=1, label=label, tfrm=tfrm)
        asm_node = asm.tree(part_tree=opm_file['part_tree'], tag='#file')
    asm_node.parent = None

    return seq, parts, part_nodes, dgm


def create_assembly_from_seq(opt_model, idx1, idx2, **kwargs):
    part_list, node_list = parttree.part_list_from_seq(opt_model, idx1, idx2)
    label = kwargs.get('label', None)
    tfrm = kwargs.get('tfrm', opt_model['seq_model'].gbl_tfrms[idx1])
    asm = Assembly(part_list, idx=idx1, label=label, tfrm=tfrm)
    asm_node = asm.tree(part_tree=opt_model['part_tree'])

    return asm, asm_node


def render_lens_shape(s1, profile1, s2, profile2, thi, extent, sd, 
                      is_flipped, hole_sd=None, apply_tfrm=True,
                      flat1_pkg=None, flat2_pkg=None):
    is_concave_s1 = s1.profile_cv < 0.0
    is_concave_s2 = s2.profile_cv > 0.0

    profile_polys = []

    flat = None
    if flat1_pkg is not None:
        do_flat1, flat1 = flat1_pkg
        if use_flat(do_flat1, is_concave_s1):
            if flat1 is None:
                flat = flat1 = compute_flat(s1, sd)
            else:
                flat = flat1

    poly1 = full_profile(profile1, is_flipped, extent, 
                         flat, hole_id=hole_sd)
    profile_polys.append(poly1)

    flat = None
    if flat2_pkg is not None:
        do_flat2, flat2 = flat2_pkg
        if use_flat(do_flat2, is_concave_s2):
            if flat2 is None:
                flat = flat2 = compute_flat(s2, sd)
            else:
                flat = flat2

    poly2 = full_profile(profile2, is_flipped, extent,
                         flat, hole_id=hole_sd, dir=-1)

    if apply_tfrm:
        # get the full transform between lens surfaces
        r_new, t_new = trns.forward_transform(s1, thi, s2)
    else:
        r_new, t_new = np.identity(3), np.array([0., 0., thi])

    # apply transformation to poly2 profile
    poly2_tfrmd = []
    for polyline in poly2:
        poly = np.array(polyline)
        poly_tfrmd = transform_poly((r_new, t_new), poly)
        poly2_tfrmd.append(poly_tfrmd.tolist())
    profile_polys.append(poly2_tfrmd)

    # assemble the tuple for handling holes, 1 or 2 polylines
    if hole_sd is None:
        poly = []
        poly += poly1[0]
        orig_pt = deepcopy(poly1[0][0])
        poly += poly2_tfrmd[0]
        poly.append(orig_pt)
        poly = poly,
    else:
        poly_list = []
        for p1, p2 in zip(poly1, poly2_tfrmd):
            poly = []
            poly += p1
            orig_pt = deepcopy(p1[0])
            poly += p2
            poly.append(orig_pt)
            poly_list.append(poly)
        poly = tuple(poly_list)

    return poly, profile_polys


def render_surf_shape(srf, profile, extent, sd, is_flipped, 
                      hole_sd=None, flat_pkg=None):
    is_concave_srf = srf.profile_cv < 0.0

    flat = None
    if flat_pkg is not None:
        do_flat1, flat1 = flat_pkg
        if use_flat(do_flat1, is_concave_srf):
            if flat1 is None:
                flat = flat1 = compute_flat(srf, sd)
            else:
                flat = flat1

    poly_list = full_profile(profile, is_flipped, extent, 
                                flat, hole_id=hole_sd)

    return poly_list


def full_profile(profile, is_flipped, edge_extent,
                 flat_id=None, hole_id=None, dir=1, steps=6):
    """Produce a 2d segmented approximation to the *profile*. 

    In the case of a hole, 2 polyline segments are returned. Otherwise, a single polyline is returned (as a tuple of len=1)

    Args:

        profile: optical profile to be sampled
        is_flipped: the flipped state of the profile
        edge_extent: tuple with symmetric or asymetric bounds
        flat_id: if not None, inside diameter of flat zone
        hole_id: if not None, inside diameter of centered surface hole
        dir: sampling direction, +1 for up, -1 for down
        steps: number of profile curve samples

    Returns:
        a tuple of 2d coord lists. a tuple is returned even in the case of a single coord list

    """
    from rayoptics.raytr.traceerror import TraceError
    def flip_profile(prf):
        return [[-pt[0], pt[1]] for pt in prf]

    if len(edge_extent) == 1:
        sd_upr = edge_extent[0]
        sd_lwr = -edge_extent[0]
    else:
        sd_upr = edge_extent[1]
        sd_lwr = edge_extent[0]

    if flat_id is None:
        if hole_id is None:
            prf = profile.profile((sd_lwr, sd_upr), dir, steps),
        else:
            prf_lwr = profile.profile((sd_lwr, -hole_id), dir, steps)
            prf_upr = profile.profile((hole_id, sd_upr), dir, steps)
            prf = prf_lwr, prf_upr

    else:
        prf = []
        # compute top part of flat
        try:
            sag = profile.sag(0, flat_id)
        except TraceError:
            sag = None
        else:
            sd_lwr_pt = [sag, sd_lwr]
            sd_upr_pt = [sag, sd_upr]

        if hole_id is None:
            if dir > 0:
                prf_full = [sd_lwr_pt] if sag is not None else []
                prf_full += profile.profile((flat_id,), dir, steps)
                if sag is not None:
                    prf_full.append(sd_upr_pt)
            else:
                prf_full = [sd_upr_pt] if sag is not None else []
                prf_full += profile.profile((flat_id,), dir, steps)
                if sag is not None:
                    prf_full.append(sd_lwr_pt)
            prf = prf_full,
        else:
            if dir > 0:
                prf_lwr = [sd_lwr_pt] if sag is not None else []
                prf_lwr += profile.profile((-flat_id, -hole_id), dir, steps)
                prf_upr = profile.profile((hole_id, flat_id), dir, steps)
                if sag is not None:
                    prf_upr.append(sd_upr_pt)
            else:
                prf_upr = [sd_upr_pt] if sag is not None else []
                prf_upr += profile.profile((hole_id, flat_id), dir, steps)
                prf_lwr = profile.profile((-flat_id, -hole_id), dir, steps)
                if sag is not None:
                    prf_lwr.append(sd_lwr_pt)
            prf = prf_lwr, prf_upr

    if is_flipped:
        if hole_id is None:
            prf = flip_profile(prf[0]),
        else:
            prf_lwr, prf_upr = prf
            prf = flip_profile(prf_lwr), flip_profile(prf_upr)

    return prf


def use_flat(do_flat, is_concave):
    if do_flat == 'always':
        return True
    elif do_flat == 'if concave' and is_concave:
        return True
    elif do_flat == 'if convex' and not is_concave:
        return True
    return False


def compute_flat(ifc, sd, under_fract=0.05):
    ca = ifc.surface_od()
    if (1.0 - ca/sd) >= under_fract:
        flat = ca
    else:
        flat = None
    return flat


def encode_obj_reference(obj, obj_attr_str, attrs):
    attrs[obj_attr_str+'_id'] = str(id(getattr(obj, obj_attr_str)))
    del attrs[obj_attr_str]


def sync_obj_reference(obj, obj_attr_str, obj_dict, alt_attr_value):
    if hasattr(obj, obj_attr_str+'_id'):
        obj_id = getattr(obj, obj_attr_str+'_id')
        setattr(obj, obj_attr_str, obj_dict[obj_id])
        delattr(obj, obj_attr_str+'_id')
    else:
        setattr(obj, obj_attr_str, alt_attr_value)


# --- Element definitions
@runtime_checkable
class Part(Protocol):
    """Abstract base class for all types of elements. """
    label_format: ClassVar[str]
    label: str
    parent: Any
    is_flipped: bool = False
    ele_token: str

    def flip(self):
        """Called by opt_model.flip when a Part is flipped. """
        self.do_flip()
        self.is_flipped = not self.is_flipped

    @abstractmethod
    def do_flip(self):
        """Subclass action when it is flipped. """
        raise NotImplementedError

    @abstractmethod
    def sync_to_restore(self, ele_model, surfs, gaps, tfrms, 
                        profile_dict, parts_dict):
        raise NotImplementedError

    @abstractmethod
    def sync_to_seq(self, seq_model: SequentialModel):
        raise NotImplementedError

    @abstractmethod
    def sync_to_ele_def(self, seq_model, ele_def):
        """ Update idx_list and gap_list according to ele_def.
        
        ele_def: (ele_type, idx_list, gap_list)
        """
        raise NotImplementedError

    @abstractmethod
    def tree(self, **kwargs) -> Node:
        raise NotImplementedError

    @abstractmethod
    def idx_list(self) -> List[int]:
        raise NotImplementedError

    @abstractmethod
    def reference_idx(self) -> int:
        raise NotImplementedError

    @abstractmethod
    def reference_interface(self) -> Interface:
        raise NotImplementedError

    @abstractmethod
    def profile_list(self) -> List[SurfaceProfile]:
        raise NotImplementedError

    @abstractmethod
    def gap_list(self) -> List[Gap]:
        raise NotImplementedError

    @abstractmethod
    def update_size(self) -> None:
        raise NotImplementedError

    @abstractmethod
    def render_shape(self) -> List[GraphicsHandle]:
        '''return a polyline that is representative of the cemented element. '''
        raise NotImplementedError

    @abstractmethod
    def render_handles(self, opt_model) -> Dict[str, GraphicsHandle]:
        raise NotImplementedError

    @abstractmethod
    def handle_actions(self) -> Dict[str, Any]:
        raise NotImplementedError


def do_flip_with_part_list(part_list: List[Part], flip_pt_tfrm) -> None:
    """Flip a list of parts around a flip_pt. """

    r_asm, flip_pt = flip_pt_tfrm
            
    r_asm_new = np.matmul(r_asm, rot_around_y)
    
    for p in part_list:
        r, t = p.tfrm
        # get part into flip_pt_tfrm coordinate system
        r_part = np.matmul(r_asm.T, r)
        t_part = r_asm.T.dot(t - flip_pt)

        # apply the 180 degree rotation to the part
        r_part_new = np.matmul(r_asm_new, r_part)
        t_part_new = flip_pt + r_asm_new.dot(t_part)

        # set the new transform and flip the is_flipped bit
        p.tfrm = r_part_new, t_part_new
        p.is_flipped = not p.is_flipped


class Element(Part):
    """Lens element domain model. Manage rendering and selection/editing.

    An Element consists of 2 Surfaces, 1 Gap, and edge_extent information.

    Attributes:
        s1: first/origin :class:`~rayoptics.seq.interface.Interface`
        s2: second/last :class:`~rayoptics.seq.interface.Interface`
        gap: element thickness and material :class:`~rayoptics.seq.gap.Gap`
        tfrm: global transform to element origin, (Rot3, trans3)
        medium_name: the material filling the gap
        flat1, flat2: semi-diameter of flat or None. Setting to None will 
                      result in re-evaluation of flat ID
        do_flat1, do_flat2: 'if concave', 'always', 'never', 'if convex'
        handles: dict of graphical entities
        actions: dict of actions associated with the graphical handles
    """
    clut = rgbt.RGBTable(filename='red_blue64.csv',
                         data_range=[10.0, 100.])

    label_format = 'E{}'
    serial_number = 0
    default_ele_token = 'lens'

    def __init__(self, sg_def=None, ele_def_pkg=None, tfrm=None, 
                 idx=0, idx2=1, sd=1., label=None):
        if label is None:
            Element.serial_number += 1
            self.label = Element.label_format.format(Element.serial_number)
        else:
            self.label = label
        if tfrm is not None:
            self.tfrm = tfrm
        else:
            self.tfrm = (np.identity(3), np.array([0., 0., 0.]))

        if sg_def is not None:
            s1, s2, g = sg_def
            self.s1 = s1
            self.profile1 = s1.profile
            self.s1_indx = idx
            self.s2 = s2
            self.s2_indx = idx2
            self.profile2 = s2.profile
            self.gap = g
            self.medium_name = self.gap.medium.name()
            self.ele_token = Element.default_ele_token
        elif ele_def_pkg is not None:
            seq_model, ele_def = ele_def_pkg
            self.sync_to_ele_def(seq_model, ele_def)

        self._sd = sd
        self.hole_sd = None
        self.flat1 = None
        self.flat2 = None
        self.do_flat1 = 'if concave'  # alternatives are 'never', 'always',
        self.do_flat2 = 'if concave'  # or 'if convex'
        self.handles = {}
        self.actions = {}

    @property
    def sd(self):
        """Semi-diameter """
        return self._sd

    @sd.setter
    def sd(self, semidiam):
        self._sd = semidiam
        self.edge_extent = (-semidiam, semidiam)

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['parent']
        del attrs['s1']
        del attrs['s2']
        del attrs['gap']
        del attrs['handles']
        del attrs['actions']
        encode_obj_reference(self, 'profile1', attrs)
        encode_obj_reference(self, 'profile2', attrs)
        if hasattr(self, 'profile_polys'):
            del attrs['profile_polys']
        return attrs

    def __str__(self):
        fmt = 'Element: {!r}, {!r}, t={:.4f}, sd={:.4f}, glass: {}'
        return fmt.format(self.s1.profile, self.s2.profile, self.gap.thi,
                          self.sd, self.gap.medium.name())

    def listobj_str(self):
        ele_type = self.ele_token, type(self).__module__, type(self).__name__
        idx_list = tuple(i for i in self.idx_list())
        gap_list = tuple(idx_list[i] for i, g in enumerate(self.gap_list()))

        o_str = f"{ele_type[0]}: {ele_type[2]}\n"
        o_str += f"idx={idx_list},   gaps={gap_list}   conic cnst={self.cc}\n"
        o_str += f"coefficients: {self.coefs}\n"
        return o_str
        fmt = f"Element: {self.s1.profile!r}, {self.s2.profile!r}, t={self.gap.thi:.4f}, sd={self.sd:.4f}, glass: {self.gap.medium.name()}"

    def sync_to_restore(self, ele_model, surfs, gaps, tfrms, 
                        profile_dict, parts_dict):
        # when restoring, we want to use the stored indices to look up the
        # new object instances
        self.parent = ele_model
        if not hasattr(self, 'tfrm'):
            self.tfrm = tfrms[self.s1_indx]
        self.s1 = surfs[self.s1_indx]
        sync_obj_reference(self, 'profile1', profile_dict, self.s1.profile)

        if self.is_flipped:
            self.gap = gaps[self.s2_indx]
        else:
            self.gap = gaps[self.s1_indx]
        self.s2 = surfs[self.s2_indx]
        sync_obj_reference(self, 'profile2', profile_dict, self.s2.profile)

        if not hasattr(self, 'medium_name'):
            self.medium_name = self.gap.medium.name()
        if not hasattr(self, 'ele_token'):
            self.ele_token = Element.default_ele_token

        if not hasattr(self, 'do_flat1'):
            self.do_flat1 = 'if concave'
        if not hasattr(self, 'do_flat2'):
            self.do_flat2 = 'if concave'
        if not hasattr(self, 'hole_sd'):
            self.hole_sd = None
        self.handles = {}
        self.actions = {}

    def sync_to_seq(self, seq_model):
        # when updating, we want to use the stored object instances to get the
        # current indices into the interface list (e.g. to handle insertion and
        # deletion of interfaces)
        self.s1_indx = seq_model.ifcs.index(self.s1)
        self.s2_indx = seq_model.ifcs.index(self.s2)
        self.profile1 = self.s1.profile
        self.profile2 = self.s2.profile
        self.medium_name = self.gap.medium.name()

    def sync_to_ele_def(self, seq_model, ele_def):
        ele_type, idx_list, gap_list = ele_def
        ele_token, ele_module, ele_class = ele_type
        self.ele_token = ele_token
        self.s1_indx = idx_list[0]
        self.s2_indx = idx_list[1]
        self.s1 = seq_model.ifcs[self.s1_indx]
        self.s2 = seq_model.ifcs[self.s2_indx]
        self.profile1 = self.s1.profile
        self.profile2 = self.s2.profile
        self.gap = seq_model.gaps[gap_list[0]]
        self.medium_name = self.gap.medium.name()

    def tree(self, **kwargs):
        """Build tree linking sequence to element model. """

        default_tag = '#element#lens'
        tag = default_tag + kwargs.get('tag', '')
        zdir = kwargs.get('z_dir', 1)

        # Interface branch 1
        e = Node(self.label, id=self, tag=tag)
        p1 = Node('p1', id=self.profile1, tag='#profile', parent=e)
        Node(f'i{self.s1_indx}', id=self.s1, tag='#ifc', parent=p1)

        # Gap branch
        t = Node('t', id=self.gap, tag='#thic', parent=e)
        Node(f'g{self.s1_indx}', id=(self.gap, zdir), tag='#gap', parent=t)

        # Interface branch 2
        p2 = Node('p2', id=self.profile2, tag='#profile', parent=e)
        Node(f'i{self.s2_indx}', id=self.s2, tag='#ifc', parent=p2)

        return e

    def idx_list(self):
        if hasattr(self, 'parent') and self.parent is not None:
            seq_model = self.parent.opt_model['seq_model']
            try:
                self.s1_indx = seq_model.ifcs.index(self.s1)
            except ValueError:
                self.s1_indx = str(self.s1_indx)
            try:
                self.s2_indx = seq_model.ifcs.index(self.s2)
            except ValueError:
                self.s2_indx = str(self.s2_indx)
            return [self.s1_indx, self.s2_indx]
        else:
            print(f"idx_list: {self.label}")
            return []

    def reference_idx(self):
        return self.s1_indx

    def reference_interface(self):
        return self.s1

    def profile_list(self):
        return [self.profile1, self.profile2]

    def gap_list(self):
        return [self.gap]

    def get_power(self, nom_wvl='d'):
        cv1 = self.s1.profile_cv
        cv2 = self.s2.profile_cv
        rndx = self.gap.medium.rindex(nom_wvl)
        th = self.gap.thi
        power = (rndx - 1)*(cv1 - cv2 + th*cv1*cv2*(rndx - 1)/rndx)
        return power

    def get_bending(self):
        cv1 = self.s1.profile_cv
        cv2 = self.s2.profile_cv
        delta_cv = cv1 - cv2
        bending = 0.
        if delta_cv != 0.0:
            bending = (cv1 + cv2)/delta_cv
        return bending

    def set_bending(self, bending):
        power = self.get_power()
        lens = lens_from_power(power=power, bending=bending, th=self.gap.thi,
                               med=self.gap.medium)
        cv1_new, cv2_new, _ = lens

        self.s1.profile_cv = cv1_new
        self.s2.profile_cv = cv2_new

    def do_flip(self):
        r, t = self.tfrm
        thi = self.gap.thi
        if self.is_flipped:
            r_new = np.matmul(rot_around_y, r).T
            t_new = t - r_new.dot(np.array([0, 0, thi]))
        else:
            t_new = t + r.dot(np.array([0, 0, thi]))
            r_new = np.matmul(r, rot_around_y)
        self.tfrm = r_new, t_new

    def update_size(self):
        extents = np.union1d(self.s1.get_y_aperture_extent(),
                             self.s2.get_y_aperture_extent())
        self.edge_extent = (extents[0], extents[-1])
        self.sd = max(self.s1.surface_od(), self.s2.surface_od())
        return self.sd

    def extent(self):
        if hasattr(self, 'edge_extent'):
            return self.edge_extent
        else:
            return (-self.sd, self.sd)

    def render_shape(self):
        s1 = self.s1
        s2 = self.s2
        flat1_pkg = self.do_flat1, self.flat1
        flat2_pkg = self.do_flat2, self.flat2
        poly, self.profile_polys = render_lens_shape(
            s1, s1.profile, s2, s2.profile, self.gap.thi,
            self.extent(), self.sd, self.is_flipped, 
            hole_sd=self.hole_sd, flat1_pkg=flat1_pkg, flat2_pkg=flat2_pkg)

        return poly

    def render_handles(self, opt_model):
        self.handles = {}
        thi = self.gap.thi

        shape = self.render_shape()
        color = calc_render_color_for_material(self.gap.medium)

        for i, poly in enumerate(shape):
            self.handles['shape'+str(i+1)] = GraphicsHandle(
                poly, self.tfrm, 'polygon', color
                )

        extent = self.extent()

        for i, poly_segs in enumerate(self.profile_polys, start=1):
            for poly_seg in poly_segs:
                gh = GraphicsHandle(poly_seg, self.tfrm, 'polyline')
                self.handles[f's{i}_profile'] = gh

        poly_s1 = self.profile_polys[0]
        poly_s2 = self.profile_polys[1]

        if self.hole_sd is None:
            poly_s1 = poly_s1[0]
            poly_s2 = poly_s2[0]
            poly_sd_upr = []
            poly_sd_upr.append([poly_s1[-1][0], extent[1]])
            poly_sd_upr.append([poly_s2[0][0], extent[1]])
            self.handles['sd_upr'] = GraphicsHandle(poly_sd_upr, self.tfrm,
                                                    'polyline')
            poly_sd_lwr = []
            poly_sd_lwr.append([poly_s2[-1][0], extent[0]])
            poly_sd_lwr.append([poly_s1[0][0], extent[0]])
            self.handles['sd_lwr'] = GraphicsHandle(poly_sd_lwr, self.tfrm,
                                                    'polyline')
        else:
            poly_s1_lwr, poly_s1_upr = poly_s1
            poly_s2_lwr, poly_s2_upr = poly_s2
            poly_sd_upr = []
            poly_sd_upr.append([poly_s1_upr[-1][0], extent[1]])
            poly_sd_upr.append([poly_s2_upr[0][0], extent[1]])
            self.handles['sd_upr'] = GraphicsHandle(poly_sd_upr, self.tfrm,
                                                    'polyline')
            poly_sd_lwr = []
            poly_sd_lwr.append([poly_s2_lwr[-1][0], extent[0]])
            poly_sd_lwr.append([poly_s1_lwr[0][0], extent[0]])
            self.handles['sd_lwr'] = GraphicsHandle(poly_sd_lwr, self.tfrm,
                                                    'polyline')
            
        poly_ct = []
        poly_ct.append([0., 0.])
        poly_ct.append([thi, 0.])
        self.handles['ct'] = GraphicsHandle(poly_ct, self.tfrm, 'polyline')

        return self.handles

    def handle_actions(self):
        self.actions = {}

        shape_actions = {}
        shape_actions['pt'] = BendAction(self)
        shape_actions['y'] = AttrAction(self, 'sd')
        shape_actions['glass'] = ReplaceGlassAction(self.gap)
        self.actions['shape'] = shape_actions

        s1_prof_actions = {}
        s1_prof_actions['pt'] = SagAction(self.s1)
        self.actions['s1_profile'] = s1_prof_actions

        s2_prof_actions = {}
        s2_prof_actions['pt'] = SagAction(self.s2)
        self.actions['s2_profile'] = s2_prof_actions

        sd_upr_action = {}
        sd_upr_action['y'] = AttrAction(self, 'sd')
        self.actions['sd_upr'] = sd_upr_action

        sd_lwr_action = {}
        sd_lwr_action['y'] = AttrAction(self, 'sd')
        self.actions['sd_lwr'] = sd_lwr_action

        ct_action = {}
        ct_action['x'] = AttrAction(self.gap, 'thi')
        self.actions['ct'] = ct_action

        return self.actions

class SurfaceInterface(Part):

    label_format = 'S{}'
    serial_number = 0
    default_ele_token = 'surface'

    def __init__(self, ifc=None, ele_def_pkg=None, tfrm=None, idx=0, sd=1., 
                 z_dir=1.0, label=None):
        if label is None:
            SurfaceInterface.serial_number += 1
            self.label = SurfaceInterface.label_format.format(
                SurfaceInterface.serial_number)
        else:
            self.label = label

        if tfrm is not None:
            self.tfrm = tfrm
        else:
            self.tfrm = (np.identity(3), np.array([0., 0., 0.]))
        
        if ifc is not None:
            self.s = ifc
            self.s_indx = idx
            self.profile = ifc.profile
            self.ele_token = SurfaceInterface.default_ele_token
        elif ele_def_pkg is not None:
            seq_model, ele_def = ele_def_pkg
            self.sync_to_ele_def(seq_model, ele_def)
            
        self.z_dir = z_dir
        self.sd = sd
        self.hole_sd = None
        self.flat = None
        self.do_flat = 'if concave'
        self.medium_name = 'Surface'
        self.handles = {}
        self.actions = {}

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['parent']
        del attrs['s']
        del attrs['handles']
        del attrs['actions']
        encode_obj_reference(self, 'profile', attrs)
        if hasattr(self, 'profile_polys'):
            del attrs['profile_polys']
        return attrs

    def __str__(self):
        return f"Surface: {self.profile!r}, sd={self.sd:.4f}"

    def listobj_str(self):
        o_str = f"part: {type(self).__name__}, "
        o_str += self.profile.listobj_str()
        o_str += f"sd={self.sd:.4f}\n"
        return o_str

    def sync_to_restore(self, ele_model, surfs, gaps, tfrms, 
                        profile_dict, parts_dict):
        self.parent = ele_model
        if not hasattr(self, 'tfrm'):
            self.tfrm = tfrms[self.s_indx]
        self.s = surfs[self.s_indx]
        sync_obj_reference(self, 'profile', profile_dict, self.s.profile)
        if not hasattr(self, 'medium_name'):
            self.medium_name = 'Surface'
        if not hasattr(self, 'ele_token'):
            self.ele_token = SurfaceInterface.default_ele_token
        if not hasattr(self, 'hole_sd'):
            self.hole_sd = None
        self.handles = {}
        self.actions = {}

    def sync_to_seq(self, seq_model):
        self.s_indx = seq_model.ifcs.index(self.s)
        self.profile = self.s.profile

    def sync_to_ele_def(self, seq_model, ele_def):
        ele_type, idx_list, gap_list = ele_def
        ele_token, ele_module, ele_class = ele_type
        self.ele_token = ele_token
        self.s_indx = idx_list[0]
        self.s = seq_model.ifcs[idx_list[0]]
        self.profile = self.s.profile

    def tree(self, **kwargs):
        default_label_prefix = kwargs.get('default_label_prefix', 'S')
        default_tag = kwargs.get('default_tag', '#element#surface')
        tag = default_tag + kwargs.get('tag', '')
        # Interface branch
        m = Node(default_label_prefix, id=self, tag=tag)
        p = Node('p', id=self.profile, tag='#profile', parent=m)
        Node(f'i{self.s_indx}', id=self.s, tag='#ifc', parent=p)

        # Gap branch = None

        return m

    def reference_interface(self):
        return self.s

    def reference_idx(self):
        return self.s_indx

    def idx_list(self):
        seq_model = self.parent.opt_model['seq_model']
        self.s_indx = seq_model.ifcs.index(self.s)
        return [self.s_indx]

    def profile_list(self):
        return [self.profile]

    def gap_list(self):
        return []

    def do_flip(self):
        r, t = self.tfrm
        if self.is_flipped:
            r_new = np.matmul(rot_around_y, r).T
        else:
            r_new = np.matmul(r, rot_around_y)
        self.tfrm = r_new, t

    def update_size(self):
        self.edge_extent = self.s.get_y_aperture_extent()
        self.sd = self.s.surface_od()
        return self.sd

    def extent(self):
        if hasattr(self, 'edge_extent'):
            return self.edge_extent
        else:
            self.edge_extent = self.s.get_y_aperture_extent()
            return self.edge_extent

    def render_shape(self):
        is_concave_s = self.s.profile_cv < 0.0
        
        self.profile_polys = []

        computed_flat = compute_flat(self.s, self.sd)
        if use_flat(self.do_flat, is_concave_s):
            if self.flat is None:
                flat = computed_flat
            else:
                flat = self.flat
        else:
            flat = None
        poly, = full_profile(self.profile, self.is_flipped, self.extent(), 
                             flat, hole_id=self.hole_sd)
        self.profile_polys.append(poly)

        return poly

    def render_handles(self, opt_model):
        self.handles = {}

        render_color = (158, 158, 158, 64)
        shape = self.render_shape()

        self.handles['shape'] = GraphicsHandle(
            shape, self.tfrm, 'polyline', render_color
            )

        poly_s1 = self.profile_polys[0]

        if self.hole_sd is None:
            pt_sd_upr = deepcopy(poly_s1[-1])
            self.handles['sd_upr'] = GraphicsHandle(pt_sd_upr, self.tfrm,
                                                    'vertex')
            pt_sd_lwr = deepcopy(poly_s1[0])
            self.handles['sd_lwr'] = GraphicsHandle(pt_sd_lwr, self.tfrm,
                                                    'vertex')
        else:
            poly_s1_lwr, poly_s1_upr = poly_s1
            pt_sd_upr = deepcopy(poly_s1_upr[-1])
            self.handles['sd_upr'] = GraphicsHandle(pt_sd_upr, self.tfrm,
                                                    'vertex')
            pt_sd_lwr = deepcopy(poly_s1_lwr[0])
            self.handles['sd_lwr'] = GraphicsHandle(pt_sd_lwr, self.tfrm,
                                                    'vertex')
        return self.handles

    def handle_actions(self):
        self.actions = {}

        shape_actions = {}
        shape_actions['pt'] = SagAction(self.s)
        self.actions['shape'] = shape_actions

        s_prof_actions = {}
        s_prof_actions['pt'] = SagAction(self.s)
        self.actions['s_profile'] = s_prof_actions

        sd_upr_action = {}
        sd_upr_action['y'] = AttrAction(self, 'edge_extent[1]')
        self.actions['sd_upr'] = sd_upr_action

        sd_lwr_action = {}
        sd_lwr_action['y'] = AttrAction(self, 'edge_extent[0]')
        self.actions['sd_lwr'] = sd_lwr_action

        return self.actions


class Mirror(SurfaceInterface):

    label_format = 'M{}'
    serial_number = 0
    default_ele_token = 'mirror'

    def __init__(self, ifc=None, ele_def_pkg=None, 
                 thi=None, label=None, **kwargs):
        if label is None:
            Mirror.serial_number += 1
            label = Mirror.label_format.format(Mirror.serial_number)
        
        if ifc is not None:
            super().__init__(ifc=ifc, label=label, **kwargs)
            self.ele_token = Mirror.default_ele_token
        elif ele_def_pkg is not None:
            super().__init__(ele_def_pkg=ele_def_pkg, label=label, **kwargs)

        self.render_color = (158, 158, 158, 64)
        self.thi = thi
        self.medium_name = 'Mirror'
    
    def __str__(self):
        thi = self.get_thi()
        fmt = f"Mirror: {self.profile!r}, t={thi:.4f}, sd={self.sd:.4f}"
        return fmt

    def listobj_str(self):
        o_str = f"part: {type(self).__name__}, "
        o_str += self.profile.listobj_str()
        o_str += f"t={self.get_thi():.4f}, sd={self.sd:.4f}\n"
        return o_str
    
    def sync_to_restore(self, ele_model, surfs, gaps, tfrms, 
                        profile_dict, parts_dict):
        if not hasattr(self, 'medium_name'):
            self.medium_name = 'Mirror'
        if not hasattr(self, 'ele_token'):
            self.ele_token = Mirror.default_ele_token
        super().sync_to_restore(ele_model, surfs, gaps, tfrms, 
                                profile_dict, parts_dict)

    def get_thi(self):
        thi = self.thi
        if self.thi is None:
            thi = 0.05*self.sd
        return thi
    
    def tree(self, **kwargs):
        kwargs['default_label_prefix'] = 'M'
        kwargs['default_tag'] = '#element#mirror'
        return super().tree(**kwargs)

    def substrate_offset(self):
        thi = self.get_thi()
        # We want to extend the mirror substrate along the same direction
        # of the incoming ray. The mirror's z_dir is following reflection so
        # flip the sign to get the preceding direction.
        offset = -self.z_dir*thi
        return offset

    def render_shape(self):
        s = self.s
        thi = self.substrate_offset()
        flat_pkg = self.do_flat, self.flat

        poly, self.profile_polys = render_lens_shape(
            s, s.profile, s, s.profile, thi, 
            self.extent(), self.sd, self.is_flipped, apply_tfrm=False,
            hole_sd=self.hole_sd, flat1_pkg=flat_pkg, flat2_pkg=flat_pkg)

        return poly

    def render_handles(self, opt_model):
        self.handles = {}

        shape = self.render_shape()
        for i, poly in enumerate(shape):
            self.handles['shape'+str(i+1)] = GraphicsHandle(
                poly, self.tfrm, 'polygon', self.render_color
                )

        extent = self.extent()

        poly_s1 = self.profile_polys[0]
        for poly_seg in poly_s1:
            gh1 = GraphicsHandle(poly_seg, self.tfrm, 'polyline')
            self.handles['s_profile'] = gh1

        poly_s2 = self.profile_polys[1]

        if self.hole_sd is None:
            poly_s1 = poly_s1[0]
            poly_s2 = poly_s2[0]
            poly_sd_upr = []
            poly_sd_upr.append([poly_s1[-1][0], extent[1]])
            poly_sd_upr.append([poly_s2[0][0], extent[1]])
            self.handles['sd_upr'] = GraphicsHandle(poly_sd_upr, self.tfrm,
                                                    'polyline')
            poly_sd_lwr = []
            poly_sd_lwr.append([poly_s2[-1][0], extent[0]])
            poly_sd_lwr.append([poly_s1[0][0], extent[0]])
            self.handles['sd_lwr'] = GraphicsHandle(poly_sd_lwr, self.tfrm,
                                                    'polyline')
        else:
            poly_s1_lwr, poly_s1_upr = poly_s1
            poly_s2_lwr, poly_s2_upr = poly_s2
            poly_sd_upr = []
            poly_sd_upr.append([poly_s1_upr[-1][0], extent[1]])
            poly_sd_upr.append([poly_s2_upr[0][0], extent[1]])
            self.handles['sd_upr'] = GraphicsHandle(poly_sd_upr, self.tfrm,
                                                    'polyline')
            poly_sd_lwr = []
            poly_sd_lwr.append([poly_s2_lwr[-1][0], extent[0]])
            poly_sd_lwr.append([poly_s1_lwr[0][0], extent[0]])
            self.handles['sd_lwr'] = GraphicsHandle(poly_sd_lwr, self.tfrm,
                                                    'polyline')
        return self.handles


class CementedElement(Part):
    """Cemented element domain model. Manage rendering and selection/editing.

    A CementedElement consists of 3 or more Surfaces, 2 or more Gaps, and
    edge_extent information.

    Attributes:
        idxs: list of seq_model interface indices (depends on is_flipped)
        ifcs: list of :class:`~rayoptics.seq.interface.Interface`
        gaps: list of thickness and material :class:`~rayoptics.seq.gap.Gap`
        tfrm: global transform to element origin, (Rot3, trans3)
        medium_name: the material filling the gap
        flats: semi-diameter of flat if ifc is concave, or None
        handles: dict of graphical entities
        actions: dict of actions associated with the graphical handles
    """
    clut = rgbt.RGBTable(filename='red_blue64.csv',
                         data_range=[10.0, 100.])

    label_format = 'CE{}'
    serial_number = 0
    default_ele_token = 'cemented'

    def __init__(self, ifc_list=None, ele_def_pkg=None, label=None):
        if label is None:
            CementedElement.serial_number += 1
            self.label = CementedElement.label_format.format(
                CementedElement.serial_number)
        else:
            self.label = label

        if ifc_list is not None:
            g_tfrm = ifc_list[0][4]
            if g_tfrm is not None:
                self.tfrm = g_tfrm
            else:
                self.tfrm = (np.identity(3), np.array([0., 0., 0.]))
            self.idxs = []
            self.ifcs = []
            self.profiles = []
            self.gaps = []
            for interface in ifc_list:
                i, ifc, g, z_dir, g_tfrm = interface
                self.idxs.append(i)
                self.ifcs.append(ifc)
                self.profiles.append(ifc.profile)
                if g is not None:
                    self.gaps.append(g)
            if len(self.gaps) == len(self.ifcs):
                self.gaps.pop()
            self.medium_name = self._construct_medium_name()
            self.ele_token = CementedElement.default_ele_token
        elif ele_def_pkg is not None:
            seq_model, ele_def = ele_def_pkg
            self.sync_to_ele_def(seq_model, ele_def)
            self.tfrm = (np.identity(3), np.array([0., 0., 0.]))

        self._sd = self.update_size()
        self.hole_sd = None
        self.flats = [None]*len(self.ifcs)
        self.do_flat_0 = 'if concave'  # alternatives are 'never', 'always',
        self.do_flat_k = 'if concave'  # or 'if convex'

        self.do_render_shape = True # if true, render elements, otherwise surfaces

        self.handles = {}
        self.actions = {}

    @property
    def sd(self):
        """Semi-diameter """
        return self._sd

    @sd.setter
    def sd(self, semidiam):
        self._sd = semidiam
        self.edge_extent = (-semidiam, semidiam)

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['parent']
        del attrs['ifcs']
        del attrs['gaps']
        del attrs['flats']
        del attrs['handles']
        del attrs['actions']
        attrs['profile_ids'] = [str(id(p)) for p in self.profiles]
        del attrs['profiles']
        if hasattr(self, 'profile_polys'):
            del attrs['profile_polys']

        return attrs

    def __str__(self):
        fmt = 'CementedElement: {}'
        return fmt.format(self.idxs)

    def _construct_medium_name(self):
        medium_name = ''
        for g in self.gaps:
            if medium_name != '':
                medium_name += ', '
            medium_name += g.medium.name()
        return medium_name

    def listobj_str(self):
        ele_type = self.ele_token, type(self).__module__, type(self).__name__
        idx_list = tuple(i for i in self.idx_list())
        gap_list = tuple(idx_list[i] for i, g in enumerate(self.gap_list()))

        o_str = f"{ele_type[0]}: {ele_type[2]}\n"
        o_str += f"idx={idx_list},   gaps={gap_list}\n"
        return o_str
        fmt = f"Element: {self.s1.profile!r}, {self.s2.profile!r}, t={self.gap.thi:.4f}, sd={self.sd:.4f}, glass: {self.gap.medium.name()}"

    def sync_to_restore(self, ele_model, surfs, gaps, tfrms, 
                        profile_dict, parts_dict):
        # when restoring, we want to use the stored indices to look up the
        # new object instances
        self.parent = ele_model
        if not hasattr(self, 'tfrm'):
            self.tfrm = tfrms[self.idxs[0]]
            
        self.ifcs = [surfs[i] for i in self.idxs]
        if hasattr(self, 'profile_ids'):
            self.profiles = []
            for p_id in self.profile_ids:
                self.profiles.append(profile_dict[p_id])
            delattr(self, 'profile_ids')
        else:
            self.profiles = [ifc.profile for ifc in self.ifcs]
        if self.is_flipped:
            self.gaps = [gaps[i] for i in self.idxs[1:]]
        else:
            self.gaps = [gaps[i] for i in self.idxs[:-1]]
        self.flats = [None]*len(self.ifcs)
        if not hasattr(self, 'do_flat_0'):
            self.do_flat_0 = 'if concave'
        if not hasattr(self, 'do_flat_k'):
            self.do_flat_k = 'if concave'
        if not hasattr(self, 'hole_sd'):
            self.hole_sd = None
        if not hasattr(self, 'medium_name'):
            self.medium_name = self._construct_medium_name()
        if not hasattr(self, 'ele_token'):
            self.ele_token = CementedElement.default_ele_token
        if not hasattr(self, 'do_render_shape'):
            self.do_render_shape = True
        self.handles = {}
        self.actions = {}

    def sync_to_seq(self, seq_model):
        # when updating, we want to use the stored object instances to get the
        # current indices into the interface list (e.g. to handle insertion and
        # deletion of interfaces)
        self.idxs = [seq_model.ifcs.index(ifc) for ifc in self.ifcs]
        self.profiles = [ifc.profile for ifc in self.ifcs]

    def sync_to_ele_def(self, seq_model, ele_def):
        ele_type, idx_list, gap_list = ele_def
        ele_token, ele_module, ele_class = ele_type
        self.ele_token = ele_token

        if self.ele_token == 'cemented':
            self.idxs = [idx for idx in idx_list]
            self.ifcs = [seq_model.ifcs[idx] for idx in idx_list]
            self.profiles = [ifc.profile for ifc in self.ifcs]

            self.gaps = [seq_model.gaps[i] for i in gap_list]
            self.medium_name = self._construct_medium_name()

        elif self.ele_token == 'mangin':
            num_gaps = (idx_list[-1] - idx_list[0])>>1
            self.idxs = [idx for idx in idx_list]
            self.ifcs = [seq_model.ifcs[idx] for idx in idx_list]
            self.profiles = [ifc.profile for ifc in self.ifcs[:num_gaps+1]]

            self.gaps = [seq_model.gaps[i] for i in gap_list]
            self.medium_name = self._construct_medium_name()

    def tree(self, **kwargs):
        default_tag = '#element#cemented'
        tag = default_tag + kwargs.get('tag', '')
        zdir = kwargs.get('z_dir', 1)
        ce = Node(self.label, id=self, tag=tag)
        if self.ele_token == 'cemented':
            for i, sg in enumerate(zip_longest(self.ifcs, self.gaps)):
                i1 = i + 1
                ifc, gap = sg
                pid = f'p{i1}'
                p = Node(pid, id=self.profiles[i], tag='#profile', parent=ce)
                Node(f'i{self.idxs[i]}', id=ifc, tag='#ifc', parent=p)
                # Gap branch
                if gap is not None:
                    t = Node(f't{i1}', id=gap, tag='#thic', parent=ce)
                    Node(f'g{self.idxs[i]}', id=(gap, zdir),
                        tag='#gap', parent=t)
        elif self.ele_token == 'mangin':
            idx1 = self.idxs[0]
            idxk = self.idxs[-1]
            len_ifcs = len(self.ifcs)
            len_gaps = len(self.gaps)
            num_gaps = (idxk-idx1)>>1
            for i in range(num_gaps):
                i1 = i + 1
                ifc = self.ifcs[i]
                gap = self.gaps[i]
                pid = f'p{i1}'
                p = Node(pid, id=self.profiles[i], tag='#profile', parent=ce)
                Node(f'i{self.idxs[i]}', id=self.ifcs[i], tag='#ifc', parent=p)
                i2 = len_ifcs - i1
                Node(f'i{self.idxs[i2]}', id=self.ifcs[i2], 
                     tag='#ifc', parent=p)
                # Gap branch
                t = Node(f't{i1}', id=self.gaps[i], tag='#thic', parent=ce)
                Node(f'g{self.idxs[i]}', id=(self.gaps[i], zdir),
                    tag='#gap', parent=t)
                i2 = len_gaps - i1
                Node(f'g{self.idxs[i2]}', id=(self.gaps[i2], -zdir),
                    tag='#gap', parent=t)
            i += 1
            i1 += 1
            p = Node(f'p{i1}', id=self.profiles[i], tag='#profile', parent=ce)
            Node(f'i{self.idxs[i]}', id=self.ifcs[i], tag='#ifc', parent=p)

        return ce

    def idx_list(self):
        seq_model = self.parent.opt_model['seq_model']
        self.idxs = [seq_model.ifcs.index(ifc) for ifc in self.ifcs]
        return self.idxs

    def reference_idx(self):
        return self.idxs[0]

    def reference_interface(self):
        return self.ifcs[0]

    def profile_list(self):
       return self.profiles

    def gap_list(self):
        return self.gaps

    def do_flip(self):
        r, t = self.tfrm
        thi = 0
        for g in self.gaps:
            thi += g.thi
        if self.is_flipped:
            r_new = np.matmul(rot_around_y, r).T
            t_new = t - r_new.dot(np.array([0, 0, thi]))
        else:
            t_new = t + r.dot(np.array([0, 0, thi]))
            r_new = np.matmul(r, rot_around_y)
        self.tfrm = r_new, t_new

    def update_size(self):
        extents = np.union1d(self.ifcs[0].get_y_aperture_extent(),
                             self.ifcs[-1].get_y_aperture_extent())
        self.edge_extent = (extents[0], extents[-1])
        self.sd = max([ifc.surface_od() for ifc in self.ifcs])
        return self.sd

    def compute_inner_flat(self, idx, sd, k):
        ''' compute flats, if needed, for the inner cemented surfaces. 

        Args:
            idx: index of inner surface in profile list
            sd: the semi-diameter of the cemented element
            k: final, k-th, profile index
        
        This function is needed to handle the cases where one of the outer
        surfaces has a flat and the inner surface would intersect this flat.
        All inner cemented surfaces are assumed to be spherical.
        See model US007277232_Example04P.roa
        See also cv_fisheye.roa
        '''
        def sphere_sag_to_zone(sag, c):
            if c == 0.:
                return sd
            else:
                R = abs(1/c)
                try:
                    zone = sqrt(2*sag*R - sag**2)
                except ValueError:
                    zone = R
            return zone

        p = self.profiles[idx]
        R = 0 if p.cv == 0 else abs(1/p.cv)
        ifc = self.ifcs[idx]
        ca = ifc.surface_od()

        flat_0 = self.flats[0]
        sag0 = self.profiles[0].sag(0., flat_0) if flat_0 else 0.0
        flat_k = self.flats[k]
        sagk = self.profiles[k].sag(0., flat_k) if flat_k else 0.0
        
        thi_b4 = thi_aftr = 0.
        for i in range(idx):
            thi_b4 += self.gaps[i].thi
        for i in range(idx, len(self.gaps)):
            thi_aftr += self.gaps[i].thi

        flat_i = None
        if p.cv < 0.0:
            # if there's a first flat, check for intersection
            if self.flats[0] is not None:
                flat_i = sphere_sag_to_zone(sag0 + thi_b4, p.cv)
            # check if radius is smaller than semi-diameter, add flat if needed
            elif R < sd:
                flat_i = ca if ca < R else R
        elif p.cv > 0.0:
            # if there's a last flat, check for intersection
            if self.flats[k] is not None:
                flat_i = sphere_sag_to_zone(sagk + thi_aftr, p.cv)

        # check if radius is smaller than semi-diameter, add flat if needed
        if R < sd:
            flat_i = ca if ca < R else R

        if flat_i is not None and flat_i > sd:
            flat_i = sd

        return flat_i

    def extent(self):
        if hasattr(self, 'edge_extent'):
            return self.edge_extent
        else:
            return (-self.sd, self.sd)

    def render_shape(self):
        '''return a tuple of polylines of the lenses of the cemented element. 
        
        Flats on the outer surfaces of the cemented assembly are checked for 
        intersections with the internal surfaces. The first outer surface
        is the zeroth interface; the other outer surface is interface `k`.
        In the case of the mangin assembly, the k-th surface is in the middle
        of the interface list.
        '''
        ifcs = self.ifcs
        profiles = self.profiles
        gaps = self.gaps
        polygon_list = []

        len_gaps = len(gaps)
        if self.ele_token == 'mangin':
            num_gaps = len_gaps>>1
            k = num_gaps
        else:
            k = -1
        
        # examine all profiles for possible (or required) flats.
        # only consider first profile for a flat if it is concave
        is_concave_0 = self.profiles[0].cv < 0.0
        if use_flat(self.do_flat_0, is_concave_0):
            self.flats[0] = compute_flat(self.ifcs[0], self.sd)
        else:
            self.flats[0] = None
        # for the mangin case, the first and last surfaces are the
        # same. Sync the flat definitions.
        if self.ele_token == 'mangin':
            self.flats[-1] = self.flats[0]

        # only consider last profile for a flat if it is concave
        # note that "last profile" for a mangin mirror is 
        # the middle profile in the list.
        is_concave_k = self.profiles[k].cv > 0.0
        if use_flat(self.do_flat_k, is_concave_k):
            self.flats[k] = compute_flat(self.ifcs[k], self.sd)
        else:
            self.flats[k] = None

        flat1_pkg = 'always', self.flats[0]
        for i, gap in enumerate(gaps):
            # compute flats for inner profiles that intersect outer flats
            if i+1<len(gaps) and i+1 != k:
                self.flats[i+1] = self.compute_inner_flat(i+1, self.sd, k)

            if self.flats[i+1] is not None:
                flat2_pkg = 'always', self.flats[i+1]
            else:
                flat2_pkg = None

            poly_list, profile_polys = render_lens_shape(
                ifcs[i], profiles[i], ifcs[i+1], profiles[i+1], gap.thi,
                self.extent(), self.sd, self.is_flipped, hole_sd=self.hole_sd, 
                flat1_pkg=flat1_pkg, flat2_pkg=flat2_pkg
                )

            if i>0:
                r, t = trns.forward_transform(ifcs[i-1], zdist, ifcs[i])

                t_new = np.matmul(r_prev, t) + t_prev
                r_new = np.matmul(r_prev, r)
            else:
                r_new, t_new = np.identity(3), np.array([0., 0., 0.])

            poly_tfrmd_list = []
            for polyline in poly_list:
                poly = np.array(polyline)
                poly_tfrmd = transform_poly((r_new, t_new), poly)
                poly_tfrmd_list.append(poly_tfrmd.tolist())

            polygon_list.append(tuple(poly_tfrmd_list))

            zdist = gap.thi
            flat1_pkg = flat2_pkg
            r_prev, t_prev = r_new, t_new

        return tuple(polygon_list)

    def render_as_surfs(self):
        '''return a tuple of polylines of the surfaces of the cemented element. '''
        ifcs = self.ifcs
        profiles = self.profiles
        gaps = self.gaps
        len_gaps = len(gaps)
        if self.ele_token == 'mangin':
            num_gaps = len_gaps>>1
            k = num_gaps
        else:
            k = -1
        polygon_list = []

        # examine all profiles for possible (or required) flats.
        # only consider first profile for a flat if it is concave
        is_concave_0 = self.profiles[0].cv < 0.0
        if use_flat(self.do_flat_0, is_concave_0):
            self.flats[0] = compute_flat(self.ifcs[0], self.sd)
        else:
            self.flats[0] = None

        # only consider last profile for a flat if it is concave
        is_concave_k = self.profiles[-1].cv > 0.0
        if use_flat(self.do_flat_k, is_concave_k):
            self.flats[-1] = compute_flat(self.ifcs[-1], self.sd)
        else:
            self.flats[-1] = None

        for i, ifc in enumerate(ifcs):
            # compute flats for inner profiles that intersect outer flats
            self.flats[i] = self.compute_inner_flat(i, self.sd, k)

            if self.flats[i] is not None:
                flat_pkg = 'always', self.flats[i]
            else:
                flat_pkg = None

            poly_list = render_surf_shape(ifc, profiles[i], self.extent(),
                                          self.sd, self.is_flipped, 
                                          hole_sd=self.hole_sd, 
                                          flat_pkg=flat_pkg)

            if i>0:
                r, t = trns.forward_transform(b4_ifc, zdist, ifc)

                t_new = np.matmul(r_prev, t) + t_prev
                r_new = np.matmul(r_prev, r)
            else:
                r_new, t_new = np.identity(3), np.array([0., 0., 0.])

            poly_tfrmd_list = []
            for polyline in poly_list:
                poly = np.array(polyline)
                poly_tfrmd = transform_poly((r_new, t_new), poly)
                poly_tfrmd_list.append(poly_tfrmd.tolist())

            polygon_list.append(tuple(poly_tfrmd_list))

            if i<len(gaps):
                zdist = gaps[i].thi
            b4_ifc = ifc
            r_prev, t_prev = r_new, t_new

        return tuple(polygon_list)
    
    def render_handles(self, opt_model):
        self.handles = {}

        if self.do_render_shape:
            shape = self.render_shape()
        else:
            shape = self.render_as_surfs()

        self.handles['shape'] = GraphicsHandle(shape, self.tfrm, 'polyline')

        if self.do_render_shape:
            for i, gap in enumerate(self.gaps):
                polygon_list = shape[i]
                color = calc_render_color_for_material(gap.medium)
                for j, poly in enumerate(polygon_list):
                    self.handles['shape'+f"{i+1}"+f"{j+1}"] = GraphicsHandle(
                        poly, self.tfrm, 'polygon', color
                        )

        return self.handles
    
    def handle_actions(self):
        self.actions = {}

        return self.actions


class ThinElement(Part):

    label_format = 'TL{}'
    serial_number = 0
    default_ele_token = 'thin_lens'

    def __init__(self, ifc=None, ele_def_pkg=None, tfrm=None, idx=0, sd=None,
                 label=None):
        if label is None:
            ThinElement.serial_number += 1
            self.label = ThinElement.label_format.format(
                ThinElement.serial_number)
        else:
            self.label = label

        if ifc is not None:
            self.intrfc = ifc
            self.intrfc_indx = idx
            self.medium_name = 'Thin Element'
            self.ele_token = ThinElement.default_ele_token
        elif ele_def_pkg is not None:
            seq_model, ele_def = ele_def_pkg
            self.sync_to_ele_def(seq_model, ele_def)

        self.render_color = (192, 192, 192)
        if tfrm is not None:
            self.tfrm = tfrm
        else:
            self.tfrm = (np.identity(3), np.array([0., 0., 0.]))
        if sd is not None:
            self.sd = sd
        else:
            self.sd = self.intrfc.max_aperture
        self.handles = {}
        self.actions = {}

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['parent']
        del attrs['intrfc']
        del attrs['handles']
        del attrs['actions']
        return attrs

    def __str__(self):
        return str(self.intrfc)

    def tree(self, **kwargs):
        default_tag = '#element#thinlens'
        tag = default_tag + kwargs.get('tag', '')
        tle = Node('TL', id=self, tag=tag)
        Node('tl', id=self.intrfc, tag='#ifc', parent=tle)
        return tle

    def sync_to_restore(self, ele_model, surfs, gaps, tfrms, 
                        profile_dict, parts_dict):
        self.parent = ele_model
        if not hasattr(self, 'tfrm'):
            self.tfrm = tfrms[self.intrfc_indx]
        self.intrfc = surfs[self.intrfc_indx]
        if not hasattr(self, 'medium_name'):
            self.medium_name = 'Thin Element'
        if not hasattr(self, 'ele_token'):
            self.ele_token = ThinElement.default_ele_token
        self.handles = {}
        self.actions = {}

        ro_version = ele_model.opt_model.ro_version
        if version.parse(ro_version) < version.parse("0.7.0a"):
            ThinElement.serial_number += 1
            self.label = ThinElement.label_format.format(ThinElement.serial_number)

    def sync_to_seq(self, seq_model):
        self.intrfc_indx = seq_model.ifcs.index(self.intrfc)

    def sync_to_ele_def(self, seq_model, ele_def):
        ele_type, idx_list, gap_list = ele_def
        ele_token, ele_module, ele_class = ele_type
        self.ele_token = ele_token
        self.intrfc_indx = idx_list[0]
        self.intrfc = seq_model.ifcs[idx_list[0]]

    def reference_interface(self):
        return self.intrfc

    def reference_idx(self):
        return self.intrfc_indx

    def profile_list(self):
        return []

    def idx_list(self):
        seq_model = self.parent.opt_model['seq_model']
        self.intrfc_indx = seq_model.ifcs.index(self.intrfc)
        return [self.intrfc_indx]

    def gap_list(self):
        return []

    def do_flip(self):
        r, t = self.tfrm
        if self.is_flipped:
            r_new = np.matmul(rot_around_y, r).T
        else:
            r_new = np.matmul(r, rot_around_y)
        self.tfrm = r_new, t

    def update_size(self):
        self.sd = self.intrfc.surface_od()
        return self.sd

    def render_shape(self):
        poly = self.intrfc.full_profile((-self.sd, self.sd))
        return poly

    def render_handles(self, opt_model):
        self.handles = {}
        shape = self.render_shape()
        self.handles['shape'] = GraphicsHandle(shape, self.tfrm, 'polyline', 
                                               self.render_color)
        return self.handles

    def handle_actions(self):
        self.actions = {}
        return self.actions

class DummyInterface(Part):

    label_format = 'D{}'
    serial_number = 0
    default_ele_token = 'dummy'

    def __init__(self, ifc=None, ele_def_pkg=None, idx=0, sd=None, tfrm=None,
                 label=None):
        if label is None:
            DummyInterface.serial_number += 1
            self.label = DummyInterface.label_format.format(
                DummyInterface.serial_number)
        else:
            self.label = label

        self.render_color = (192, 192, 192)
        if tfrm is not None:
            self.tfrm = tfrm
        else:
            self.tfrm = (np.identity(3), np.array([0., 0., 0.]))

        self.ele_token = DummyInterface.default_ele_token
        if ifc is not None:
            self.ref_ifc = ifc
            self.idx = idx
            self.profile = ifc.profile
        elif ele_def_pkg is not None:
            seq_model, ele_def = ele_def_pkg
            self.sync_to_ele_def(seq_model, ele_def)

        self.medium_name = 'Interface'
        if sd is not None:
            self.sd = sd
        else:
            self.sd = self.ref_ifc.max_aperture
        self.handles = {}
        self.actions = {}

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['parent']
        del attrs['ref_ifc']
        del attrs['handles']
        del attrs['actions']
        encode_obj_reference(self, 'profile', attrs)
        return attrs

    def __str__(self):
        return str(self.ref_ifc)

    def sync_to_restore(self, ele_model, surfs, gaps, tfrms, 
                        profile_dict, parts_dict):
        self.parent = ele_model
        if not hasattr(self, 'tfrm'):
            self.tfrm = tfrms[self.idx]
        self.ref_ifc = surfs[self.idx]
        sync_obj_reference(self, 'profile', profile_dict, self.ref_ifc.profile)
        if not hasattr(self, 'medium_name'):
            self.medium_name = 'Interface'
        if not hasattr(self, 'ele_token'):
            if self.label == 'Object':
                self.ele_token = 'object'
            elif self.label == 'Image':
                self.ele_token = 'image'
            else:
                self.ele_token = DummyInterface.default_ele_token
        self.handles = {}
        self.actions = {}

    def sync_to_seq(self, seq_model):
        self.idx = seq_model.ifcs.index(self.ref_ifc)
        self.profile = self.ref_ifc.profile

    def sync_to_ele_def(self, seq_model, ele_def):
        ele_type, idx_list, gap_list = ele_def
        ele_token, ele_module, ele_class = ele_type
        self.ele_token = ele_token
        if ele_token == 'object':
            self.label = 'Object'
        elif ele_token == 'image':
            self.label = 'Image'
        self.idx = idx_list[0]
        self.ref_ifc = seq_model.ifcs[self.idx]
        self.profile = self.ref_ifc.profile

    def tree(self, **kwargs):
        default_tag = '#dummyifc'
        if self.ele_token == 'object':
            default_tag += '#object'
        elif self.ele_token == 'image':
            default_tag += '#image'
        tag = default_tag + kwargs.get('tag', '')
        di = Node('DI', id=self, tag=tag)
        p = Node('p', id=self.profile, tag='#profile', parent=di)
        Node(f'i{self.idx}', id=self.ref_ifc, tag='#ifc', parent=p)
        return di

    def reference_interface(self):
        return self.ref_ifc

    def reference_idx(self):
        return self.idx

    def interface_list(self):
        return [self.ref_ifc]

    def profile_list(self):
        return [self.profile]

    def idx_list(self):
        seq_model = self.parent.opt_model['seq_model']
        self.idx = seq_model.ifcs.index(self.ref_ifc)
        return [self.idx]

    def gap_list(self):
        return []

    def do_flip(self):
        r, t = self.tfrm
        if self.is_flipped:
            r_new = np.matmul(rot_around_y, r).T
        else:
            r_new = np.matmul(r, rot_around_y)
        self.tfrm = r_new, t

    def update_size(self):
        self.sd = self.ref_ifc.surface_od()
        return self.sd

    def render_shape(self):
        poly, = full_profile(self.profile, self.is_flipped, 
                             (-self.sd, self.sd))
        return poly

    def render_handles(self, opt_model):
        self.handles = {}

        shape = self.render_shape()
        self.handles['shape'] = GraphicsHandle(shape, self.tfrm, 'polyline')

        return self.handles

    def handle_actions(self):
        self.actions = {}

        def get_adj_spaces():
            seq_model = self.parent.opt_model.seq_model
            if self.idx > 0:
                before = seq_model.gaps[self.idx-1].thi
            else:
                before = None
            if self.idx < seq_model.get_num_surfaces() - 1:
                after = seq_model.gaps[self.idx].thi
            else:
                after = None
            return (before, after)

        def set_adj_spaces(cur_value, change):
            seq_model = self.parent.opt_model.seq_model
            if cur_value[0] is not None:
                seq_model.gaps[self.idx-1].thi = cur_value[0] + change
            if cur_value[1] is not None:
                seq_model.gaps[self.idx].thi = cur_value[1] - change

        slide_action = {}
        slide_action['x'] = Action(get_adj_spaces, set_adj_spaces)
        self.actions['shape'] = slide_action

        return self.actions


class Space(Part):

    label_format = 'SP{}'
    serial_number = 0
    default_ele_token = 'space'

    def __init__(self, label=None, g=None, ele_def_pkg=None, idx=0, tfrm=None, 
                 z_dir=1, **kwargs):
        if label is None:
            Space.serial_number += 1
            self.label = Space.label_format.format(Space.serial_number)
        else:
            self.label = label

        if tfrm is not None:
            self.tfrm = tfrm
        else:
            self.tfrm = (np.identity(3), np.array([0., 0., 0.]))

        if g is not None:
            self.gap = g
            self.idx = idx
            self.s1 = None
            self.s2 = None
            self.medium_name = self.gap.medium.name()
            self.ele_token = Space.default_ele_token
        elif ele_def_pkg is not None:
            seq_model, ele_def = ele_def_pkg
            self.sync_to_ele_def(seq_model, ele_def)
            
        # self.render_color = (237, 243, 254, 64)  # light blue
        self.render_color = (0, 243, 0, 64)  # light blue
        self.z_dir = z_dir

        self.handles = {}
        self.actions = {}

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['parent']
        del attrs['s1']
        del attrs['s2']
        del attrs['gap']
        del attrs['handles']
        del attrs['actions']
        return attrs

    def __str__(self):
        return str(self.gap)

    def sync_to_restore(self, ele_model, surfs, gaps, tfrms, 
                        profile_dict, parts_dict):
        self.parent = ele_model
        self.gap = gaps[self.idx]
        self.s1 = surfs[self.idx]
        self.s2 = surfs[self.idx+1]
        if not hasattr(self, 'tfrm'):
            self.tfrm = tfrms[self.idx]
        if not hasattr(self, 'render_color'):
            self.render_color = (237, 243, 254, 64)  # light blue
        if not hasattr(self, 'medium_name'):
            self.medium_name = self.gap.medium.name()
        if not hasattr(self, 'ele_token'):
            self.ele_token = Space.default_ele_token
        self.handles = {}
        self.actions = {}

    def sync_to_seq(self, seq_model):
        self.idx = idx = seq_model.gaps.index(self.gap)
        self.z_dir = seq_model.z_dir[self.idx]
        self.s1 = seq_model.ifcs[idx]
        self.s2 = seq_model.ifcs[idx+1]

    def sync_to_ele_def(self, seq_model, ele_def):
        ele_type, idx_list, gap_list = ele_def
        ele_token, ele_module, ele_class = ele_type
        self.ele_token = ele_token
        self.idx = idx = gap_list[0]
        self.gap = seq_model.gaps[self.idx]
        self.medium_name = self.gap.medium.name()
        self.s1 = seq_model.ifcs[idx]
        self.s2 = seq_model.ifcs[idx+1]

    def tree(self, **kwargs):
        default_label_prefix = kwargs.get('default_label_prefix', 'SP')
        default_tag = kwargs.get('default_tag', '#space')
        if hasattr(self, 'parent'):
            seq_model = self.parent.opt_model['seq_model']
            if self.idx == 0:
                default_tag += '#object'
            elif self.idx == len(seq_model.gaps)-1:
                default_tag += '#image'
        tag = default_tag + kwargs.get('tag', '')
        sp = Node(default_label_prefix, id=self, tag=tag)
        t = Node('t', id=self.gap, tag='#thic', parent=sp)
        zdir = kwargs.get('z_dir', self.z_dir)
        Node(f'g{self.idx}', id=(self.gap, zdir), tag='#gap', parent=t)
        return sp

    def reference_interface(self):
        return None

    def reference_idx(self):
        return self.idx

    def profile_list(self):
        return []

    def idx_list(self):
        return []

    def gap_list(self):
        return [self.gap]

    def do_flip(self):
        r, t = self.tfrm
        thi = self.gap.thi
        if self.is_flipped:
            r_new = np.matmul(rot_around_y, r).T
            t_new = t - r_new.dot(np.array([0, 0, thi]))
        else:
            t_new = t + r.dot(np.array([0, 0, thi]))
            r_new = np.matmul(r, rot_around_y)
        self.tfrm = r_new, t_new

    def update_size(self):
        if self.s1 is not None and self.s2 is not None:
            extents = np.union1d(self.s1.get_y_aperture_extent(),
                                self.s2.get_y_aperture_extent())
            self.edge_extent = (extents[0], extents[-1])
            self.sd = max(self.s1.surface_od(), self.s2.surface_od())
            return self.sd
        else:
            return None

    def extent(self):
        if hasattr(self, 'edge_extent'):
            return self.edge_extent
        else:
            return (-self.sd, self.sd)

    def render_shape(self):
        s1 = self.s1
        s2 = self.s2
        poly_pkg, profile_polys = render_lens_shape(
            s1, s1.profile, s2, s2.profile, self.gap.thi, 
            self.extent(), self.sd, self.is_flipped
            )
        poly, = poly_pkg
        return poly

    def render_handles(self, opt_model):
        self.handles = {}

        shape = self.render_shape()
        color = calc_render_color_for_material(self.gap.medium)
        self.handles['shape'] = GraphicsHandle(shape, self.tfrm, 'polygon',
                                               color)

        poly_ct = []
        poly_ct.append([0., 0.])
        poly_ct.append([self.gap.thi, 0.])

        # Modify the tfrm to account for any decenters following
        #  the reference ifc.
        tfrm = self.tfrm
        decenter = opt_model.seq_model.ifcs[self.idx].decenter
        if decenter is not None:
            r_global, t_global = tfrm
            r_after_ifc, t_after_ifc = decenter.tform_after_surf()
            t = r_global.dot(t_after_ifc) + t_global
            r = r_global if r_after_ifc is None else r_global.dot(r_after_ifc)
            tfrm = r, t

        self.handles['ct'] = GraphicsHandle(poly_ct, tfrm, 'polyline')

        return self.handles

    def handle_actions(self):
        self.actions = {}

        ct_action = {}
        ct_action['x'] = AttrAction(self.gap, 'thi')
        self.actions['ct'] = ct_action

        return self.actions


class AirGap(Space):
    label_format = 'AG{}'
    serial_number = 0
    default_ele_token = 'air'

    def __init__(self, g=None, ele_def_pkg=None, label=None, **kwargs):
        if label is None:
            AirGap.serial_number += 1
            label = AirGap.label_format.format(AirGap.serial_number)
        
        self.ele_token = AirGap.default_ele_token
        if g is not None:
            super().__init__(g=g, label=label, **kwargs)
        elif ele_def_pkg is not None:
            super().__init__(ele_def_pkg=ele_def_pkg, label=label, **kwargs)
    
    def sync_to_restore(self, ele_model, surfs, gaps, tfrms, 
                        profile_dict, parts_dict):
        if not hasattr(self, 'ele_token'):
            self.ele_token = AirGap.default_ele_token
        super().sync_to_restore(ele_model, surfs, gaps, tfrms, 
                                profile_dict, parts_dict)
    def tree(self, **kwargs):
        kwargs['default_label_prefix'] = 'AG'
        kwargs['default_tag'] = '#space#airgap'
        return super().tree(**kwargs)

    def render_handles(self, opt_model):
        self.handles = {}

        poly_ct = []
        poly_ct.append([0., 0.])
        poly_ct.append([self.gap.thi, 0.])

        # Modify the tfrm to account for any decenters following
        #  the reference ifc.
        tfrm = self.tfrm
        decenter = opt_model.seq_model.ifcs[self.idx].decenter
        if decenter is not None:
            r_global, t_global = tfrm
            r_after_ifc, t_after_ifc = decenter.tform_after_surf()
            t = r_global.dot(t_after_ifc) + t_global
            r = r_global if r_after_ifc is None else r_global.dot(r_after_ifc)
            tfrm = r, t

        self.handles['ct'] = GraphicsHandle(poly_ct, tfrm, 'polyline')

        return self.handles

class Assembly(Part):

    label_format = 'ASM{}'
    serial_number = 0
    default_ele_token = 'asm'

    def __init__(self, part_list, idx=0, tfrm=None, label=None):
        if label is None:
            Assembly.serial_number += 1
            self.label = Assembly.label_format.format(Assembly.serial_number)
        else:
            self.label = label

        if tfrm is not None:
            self.tfrm = tfrm
        else:
            self.tfrm = (np.identity(3), np.array([0., 0., 0.]))

        self.parts = part_list
        self.idx = idx
        self.ele_token = Assembly.default_ele_token
        self.handles = {}
        self.actions = {}

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['parent']
        del attrs['parts']
        del attrs['handles']
        del attrs['actions']
        part_ids = [str(id(p)) for p in self.parts]
        attrs['part_ids'] = part_ids
        return attrs

    def __str__(self):
        part_labels = [p.label for p in self.parts]
        return f"{self.label}: {part_labels}"

    def sync_to_restore(self, ele_model, surfs, gaps, tfrms, 
                        profile_dict, parts_dict):
        self.parent = ele_model
        self.parts = [parts_dict[pid] for pid in self.part_ids]
        delattr(self, 'part_ids')
        if not hasattr(self, 'ele_token'):
            self.ele_token = Assembly.default_ele_token
        self.handles = {}
        self.actions = {}

    def sync_to_seq(self, seq_model):
        self.tfrm = seq_model.gbl_tfrms[self.reference_idx()]

    def sync_to_ele_def(self, seq_model, ele_def):
        ele_type, idx_list, gap_list = ele_def
        ele_token, ele_module, ele_class = ele_type
        self.ele_token = ele_token
        part_tree = self.parent.opt_model['part_tree']

        ref_idx = idx_list[0] if len(idx_list)>0 else gap_list[0]
        for p in self.parts:
            if ref_idx == p.reference_idx():
                p_node = part_tree.node(p)
        self.medium_name = self._construct_medium_name()

    def tree(self, **kwargs):
        if 'part_tree' in kwargs:
            part_tree = kwargs.get('part_tree')
        else:
            part_tree = self.parent.opt_model['part_tree']
        default_tag = '#group#assembly'
        tag = default_tag + kwargs.get('tag', '')
        asm = Node('ASM', id=self, tag=tag)
        child_nodes = [part_tree.node(p) for p in self.parts]
        asm.children = child_nodes
        return asm

    def idx_list(self):
        idxs = []
        for p in self.parts:
            idxs += p.idx_list()
        return idxs

    def reference_idx(self):
        idxs = self.idx_list()
        self.idx = idxs[0]
        return self.idx

    def reference_interface(self):
        seq_model = self.parent.opt_model['seq_model']
        ref_idx = self.reference_idx()
        return seq_model.ifcs[ref_idx]

    def profile_list(self):
        profiles = []
        for p in self.parts:
            profiles += p.profile_list()
        return profiles

    def gap_list(self):
        gaps = []
        for p in self.parts:
            gaps += p.gap_list()
        return gaps

    def do_flip(self):
        sm = self.parent.opt_model['seq_model']
        idxs = self.idx_list()
        idx1, idx2 = idxs[0], idxs[-1]
        flip_pt = 0.5*(sm.gbl_tfrms[idx2][1] + sm.gbl_tfrms[idx1][1])
        flip_pt_tfrm = sm.gbl_tfrms[idx1][0], flip_pt
        do_flip_with_part_list(self.parts, flip_pt_tfrm)

    def update_size(self):
        pass

    def render_shape(self):
        pass

    def render_handles(self, opt_model):
        self.handles = {}
        return self.handles

    def handle_actions(self):
        self.actions = {}
        return self.actions


# --- Element model
class ElementModel:
    """Maintain the element based representation of the optical model

    Attributes:
        opt_model: the :class:`~optical.opticalmodel.OpticalModel`
        elements: list of element type things

    """

    def __init__(self, opt_model, **kwargs):
        self.opt_model = opt_model
        self.elements: List[Part] = []

    def reset(self):
        self.__init__(self.opt_model)

    def __json_encode__(self):
        attrs = dict(vars(self))
        attrs['serial_numbers'] = self.save_serial_numbers()
        del attrs['opt_model']
        del attrs['elements']
        return attrs

    def sync_to_restore(self, opt_model):
        self.opt_model = opt_model

        profile_dict = {}
        if hasattr(opt_model, 'profile_dict'):
            profile_dict = opt_model.profile_dict

        parts_dict = {}
        if not hasattr(self, 'elements'):
            if hasattr(opt_model, 'parts_dict'):
                self.elements = []
                for e in opt_model.parts_dict.values():
                    self.add_element(e)
                parts_dict = opt_model.parts_dict

        seq_model = opt_model.seq_model
        surfs = seq_model.ifcs
        gaps = seq_model.gaps
        tfrms = seq_model.compute_global_coords(1)

        if hasattr(self, 'serial_numbers'):
            self.restore_serial_numbers(self.serial_numbers)
            delattr(self, 'serial_numbers')
        else:
            self.reset_serial_numbers()

        for i, e in enumerate(self.elements, start=1):
            e.sync_to_restore(self, surfs, gaps, tfrms, 
                              profile_dict, parts_dict)
            if not hasattr(e, 'label'):
                e.label = e.label_format.format(i)
        self.sequence_elements()

    def reset_serial_numbers(self):
        Element.serial_number = 0
        Mirror.serial_number = 0
        CementedElement.serial_number = 0
        ThinElement.serial_number = 0
        SurfaceInterface.serial_number = 0
        DummyInterface.serial_number = 0
        Space.serial_number = 0
        AirGap.serial_number = 0
        Assembly.serial_number = 0

    def save_serial_numbers(self)->dict:
        serial_numbers = {
            'Element': Element.serial_number,
            'Mirror': Mirror.serial_number,
            'CementedElement': CementedElement.serial_number,
            'ThinElement': ThinElement.serial_number,
            'SurfaceInterface': SurfaceInterface.serial_number,
            'DummyInterface': DummyInterface.serial_number,
            'Space': Space.serial_number,
            'AirGap': AirGap.serial_number,
            'Assembly': Assembly.serial_number,
        }
        return serial_numbers

    def restore_serial_numbers(self, serial_numbers):
        Element.serial_number = serial_numbers.get('Element', 0)
        Mirror.serial_number = serial_numbers.get('Mirror', 0)
        CementedElement.serial_number = \
            serial_numbers.get('CementedElement', 0)
        ThinElement.serial_number = serial_numbers.get('ThinElement', 0)
        SurfaceInterface.serial_number = \
            serial_numbers.get('SurfaceInterface', 0)
        DummyInterface.serial_number = serial_numbers.get('DummyInterface', 0)
        Space.serial_number = serial_numbers.get('Space', 0)
        AirGap.serial_number = serial_numbers.get('AirGap', 0)
        Assembly.serial_number = serial_numbers.get('Assembly', 0)

    def update_model(self, **kwargs):
        """ dynamically build element list from part_tree. """
        opm = self.opt_model
    
        info = parttree.sequence_to_elements(opm['sm'], opm['em'], opm['pt'])
        
        # use the seq_model to sort the element_model list
        self.sequence_elements()

        # unless updated by the ele_model, sync it to the seq_model.
        src_model = kwargs.get('src_model', None)
        if src_model is not self:
            self.sync_to_seq(opm['sm'])

    def apply_scale_factor(self, scale_factor):
        """ Apply scale factor by resyncing with the sequential model. """
        seq_model = self.opt_model['seq_model']
        self.sync_to_seq(seq_model)

    def sync_to_seq(self, seq_model):
        """ Update element positions and ref_idx using the sequential model. """
        tfrms = seq_model.compute_global_coords(1)

        # update the elements
        for e in self.elements:
            e.update_size()
            e.sync_to_seq(seq_model)
            r, t = tfrms[e.reference_idx()]
            r_new = np.matmul(rot_around_y, r).T if e.is_flipped else r
            e.tfrm = r_new, t

    def sequence_elements(self):
        """ Sort elements in order of reference interfaces in seq_model """
        seq_model = self.opt_model.seq_model

        # sort by element reference interface sequential index
        self.elements.sort(key=lambda e: e.reference_idx()+0.5 
                           if isinstance(e, Space) else e.reference_idx())

        # Make sure z_dir matches the sequential model. Used to get
        # the correct substrate offset.
        if hasattr(seq_model, 'z_dir'):
            for e in self.elements:
                if hasattr(e, 'z_dir'):
                    e.z_dir = seq_model.z_dir[e.reference_idx()]

    def add_element(self, e: Part):
        e.parent = self
        self.elements.append(e)

    def remove_element(self, e: Part):
        e.parent = None
        self.elements.remove(e)

    def remove_node(self, e_node):
        part_tree = self.opt_model.part_tree
        nodes = part_tree.nodes_with_tag(tag='#element#airgap#dummyifc',
                                         root=e_node)
        eles = [n.id for n in nodes]
        for e in eles:
            self.remove_element(e)

    def get_num_elements(self):
        return len(self.elements)

    def list_model(self, tag: str = '#element#assembly#dummyifc'):
        nodes = self.opt_model['part_tree'].nodes_with_tag(tag=tag)
        for i, node in enumerate(nodes):
            ele = node.id
            print("%d: %s (%s): %s" %
                  (i, ele.label, type(ele).__name__, ele))

    def list_elements(self):
        for i, ele in enumerate(self.elements):
            print("%d: %s (%s): %s" %
                  (i, ele.label, type(ele).__name__, ele))

    def element_type(self, i):
        return type(self.elements[i]).__name__

    def build_ele_sg_lists(self):
        part_tag = '#element#space#airgap#dummyifc'
        nodes = self.opt_model['part_tree'].nodes_with_tag(tag=part_tag)
        eles = [n.id for n in nodes]
        ele_list = []
        ele_dict = {}
        seq_model = self.opt_model['seq_model']
        for e in eles:
            ele_type, idx_list, gap_list = build_ele_def(e, seq_model)
            ele_list.append((ele_type, idx_list, gap_list))
            ele_dict[(ele_type, idx_list, gap_list)] = e
        return ele_list, ele_dict

    def list_ele_sg(self, part_tree, seq_model):
        ele_list, ele_dict = self.build_ele_sg_lists()
        for elem in ele_list:
            ele_type, idx_list, gap_list = elem
            e = ele_dict[elem]
            print(f"{e.label}: {ele_type[0]} {idx_list} {gap_list}")


def build_ele_def(e, seq_model):
    """Package defining element info, including effect of flipping. """
    def guarded_gap_idx(g):
        try:
            return seq_model.gaps.index(g)
        except ValueError:
            return str(id(g))
    ele_type = e.ele_token, type(e).__module__, type(e).__name__
    idx_list = tuple(i for i in e.idx_list())
    gap_list = tuple(guarded_gap_idx(g) for g in e.gap_list())
    if e.is_flipped:
        # reverse the list contents to match sequential order
        idx_list = tuple(idx for idx in idx_list[::-1])
        gap_list = tuple(g for g in gap_list[::-1])
    return ele_type, idx_list, gap_list
