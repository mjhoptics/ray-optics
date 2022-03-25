#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Module for element modeling

.. Created on Sun Jan 28 16:27:01 2018

.. codeauthor: Michael J. Hayford
"""

from collections import namedtuple
import itertools
from packaging import version

from abc import abstractmethod
from typing import Protocol, ClassVar, List, Dict, Any

import numpy as np

from anytree import Node  # type: ignore

import rayoptics

import rayoptics.util.rgbtable as rgbt
import rayoptics.oprops.thinlens as thinlens
from rayoptics.elem import parttree
from rayoptics.elem.profiles import SurfaceProfile, Spherical, Conic
from rayoptics.elem.surface import Surface
from rayoptics.seq.gap import Gap
from rayoptics.seq.medium import Glass, glass_decode

# from rayoptics.optical.opticalmodel import OpticalModel
from rayoptics.seq.sequential import SequentialModel
from rayoptics.seq.interface import Interface

import rayoptics.gui.appcmds as cmds
from rayoptics.gui.actions import (Action, AttrAction, SagAction, BendAction,
                                   ReplaceGlassAction)

import opticalglass.glassfactory as gfact  # type: ignore
import opticalglass.glasspolygons as gp  # type: ignore

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
    tle = ThinElement(tl)
    tree = tle.tree()
    return [[tl, None, None, 1, +1]], [tle], tree


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
    me = Mirror(m, sd=sd, **ele_kwargs)

    tree = me.tree()

    return [[m, None, None, 1, -1]], [me], tree


def lens_from_power(power=0., bending=0., th=None, sd=1.,
                    med=None, nom_wvl='d'):
    if med is None:
        med = Glass()
    rndx = med.rindex(nom_wvl)

    if th is None:
        th = sd/5
    
    if bending == -1:
        cv2 = -power/(rndx - 1)
        cv1 = 0
    else:
        B = (bending - 1)/(bending + 1)
        a = (rndx - 1)*(th/rndx)*B
        b = 1 - B
        c = -power/(rndx - 1)
        cv1 = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
        cv2 = cv1*B

    return cv1, cv2, th, rndx, sd
   

def create_lens(power=0., bending=0., th=None, sd=1., med=None, 
                lens=None, **kwargs):
    if med is None:
        med = Glass()
    if lens is None:
        lens = lens_from_power(power=power, bending=bending, th=th, sd=sd,
                               med=med)
    else:
        cv1, cv2, th, glass, sd = lens
        med = gfact.create_glass(glass)
        rndx = med.calc_rindex('d')
        lens = cv1, cv2, th, rndx, sd

    s1 = Surface(profile=Spherical(c=cv1), max_ap=sd, delta_n=(rndx - 1))
    s2 = Surface(profile=Spherical(c=cv2), max_ap=sd, delta_n=(1 - rndx))
    g = Gap(t=th, med=med)
    le = Element(s1, s2, g, sd=sd)
    tree = le.tree()

    return [[s1, g, None, rndx, 1], [s2, None, None, 1, 1]], [le], tree


def achromat(power, Va, Vb):
    """Compute lens powers for a thin doublet achromat, given their V-numbers."""
    power_a = (Va/(Va - Vb))*power
    power_b = (Vb/(Vb - Va))*power
    return power_a, power_b


def create_cemented_doublet(power=0., bending=0., th=None, sd=1.,
                            glasses=('N-BK7,Schott', 'N-F2,Schott'),
                            **kwargs):
    from opticalglass.spectral_lines import get_wavelength  # type: ignore
    from opticalglass import glass
    wvls = np.array([get_wavelength(w) for w in ['d', 'F', 'C']])
    gla_a = gfact.create_glass(glasses[0])
    rndx_a = gla_a.calc_rindex(wvls)
    Va, PcDa = glass.calc_glass_constants(*rndx_a)
    gla_b = gfact.create_glass(glasses[1])
    rndx_b = gla_b.calc_rindex(wvls)
    Vb, PcDb = glass.calc_glass_constants(*rndx_b)

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
    ce = CementedElement(ifc_list)
    tree = ce.tree()

    return [[s1, g1, None, rndx_a, 1],
            [s2, g2, None, rndx_b, 1],
            [s3, None, None, 1, 1]], [ce], tree


def create_dummy_plane(sd=1., **kwargs):
    s = Surface(**kwargs)
    se = DummyInterface(s, sd=sd)
    tree = se.tree()
    return [[s, None, None, 1, +1]], [se], tree


def create_air_gap(t=0., **kwargs):
    g = Gap(t=t)
    ag = AirGap(g, **kwargs)
    kwargs.pop('label', None)
    tree = ag.tree(**kwargs)
    return g, ag, tree


def create_from_file(filename, **kwargs):
    opm = cmds.open_model(filename, post_process_imports=False)
    sm = opm['seq_model']
    osp = opm['optical_spec']
    em = opm['ele_model']
    pt = opm['part_tree']
    ar = opm['analysis_results']
    if len(pt.nodes_with_tag(tag='#element')) == 0:
        parttree.elements_from_sequence(em, sm, pt)

    if 'power' in kwargs:
        desired_power = kwargs['power']
        cur_power = ar['parax_data'].fod.power
        # scale_factor is linear, power is 1/linear
        #  so use reciprocal of power to compute scale_factor
        scale_factor = cur_power/desired_power
        sm.apply_scale_factor(scale_factor)

    # extract the system definition, minus object and image
    seq = [list(node) for node in sm.path(start=1, stop=-1)]
    seq[-1][1] = None
    # get the top level nodes of the input system, minus object and image
    part_nodes = pt.nodes_with_tag(tag='#element#airgap#assembly',
                                   not_tag='#object#image',
                                   node_list=pt.root_node.children)
    parts = [part_node.id for part_node in part_nodes]

    if (len(part_nodes) == 1 and '#assembly' in part_nodes[0].tag):
        asm_node = part_nodes[0]
        print("found root assembly node")
    else:
        # create an Assembly from the top level part list
        label = kwargs.get('label', None)
        tfrm = kwargs.get('tfrm', opm['seq_model'].gbl_tfrms[1])
        asm = Assembly(parts, idx=1, label=label, tfrm=tfrm)
        asm_node = asm.tree(part_tree=opm['part_tree'], tag='#file')
    asm_node.parent = None

    return seq, parts, part_nodes


def create_assembly_from_seq(opt_model, idx1, idx2, **kwargs):
    part_list, node_list = parttree.part_list_from_seq(opt_model, idx1, idx2)
    label = kwargs.get('label', None)
    tfrm = kwargs.get('tfrm', opt_model['seq_model'].gbl_tfrms[idx1])
    asm = Assembly(part_list, idx=idx1, label=label, tfrm=tfrm)
    asm_node = asm.tree(part_tree=opt_model['part_tree'])

    return asm, asm_node


def calc_render_color_for_material(matl):
    """ get element color based on V-number of glass"""
    try:
        gc = float(matl.glass_code())
    except AttributeError:
        return (255, 255, 255, 64)  # white
    else:
        # set element color based on V-number
        indx, vnbr = glass_decode(gc)
        dsg, rgb = gp.find_glass_designation(indx, vnbr)
        if rgb is None:
            return [228, 237, 243, 64]  # ED designation
#            rgb = Element.clut.get_color(vnbr)
        return rgb
    

def full_profile(profile, is_flipped, edge_extent,
                 flat_id=None, hole_id=None, dir=1, steps=6):
    do_orig = False
    if do_orig:
        return full_profile_orig(profile, is_flipped, edge_extent,
                                 flat_id, dir, steps)
    else:
        return full_profile_new(profile, is_flipped, edge_extent,
                                flat_id, hole_id, dir, steps)


def full_profile_orig(profile, is_flipped, edge_extent,
                      flat_id=None, dir=1, steps=6):
    """Produce a 2d segmented approximation to the *profile*. 

    profile: optical profile to be sampled
    edge_extent: tuple with symmetric or asymetric bounds
    flat_id: if not None, inside diameter of flat zone
    dir: sampling direction, +1 for up, -1 for down
    steps: number of profile curve samples
    """
    from rayoptics.raytr.traceerror import TraceError
    if flat_id is None:
        return profile.profile(edge_extent, dir, steps)
    else:
        if len(edge_extent) == 1:
            sd_upr = edge_extent[0]
            sd_lwr = -edge_extent[0]
        else:
            sd_upr = edge_extent[1]
            sd_lwr = edge_extent[0]
        if dir < 0:
            sd_lwr, sd_upr = sd_upr, sd_lwr

        prf = []
        try:
            sag = profile.sag(0, flat_id)
        except TraceError:
            sag = None
        else:
            prf.append([sag, sd_lwr])
        prf += profile.profile((flat_id,), dir, steps)
        if sag is not None:
            prf.append([sag, sd_upr])
        return prf


def full_profile_new(profile, is_flipped, edge_extent,
                     flat_id=None, hole_id=None, dir=1, steps=6):
    """Produce a 2d segmented approximation to the *profile*. 

    profile: optical profile to be sampled
    edge_extent: tuple with symmetric or asymetric bounds
    flat_id: if not None, inside diameter of flat zone
    hole_id: if not None, inside diameter of centered surface hole
    dir: sampling direction, +1 for up, -1 for down
    steps: number of profile curve samples
    """
    from rayoptics.raytr.traceerror import TraceError
    def process_edges(edge_extent, dir):
        if len(edge_extent) == 1:
            sd_upr = edge_extent[0]
            sd_lwr = -edge_extent[0]
        else:
            sd_upr = edge_extent[1]
            sd_lwr = edge_extent[0]
        if dir < 0:
            sd_lwr, sd_upr = sd_upr, sd_lwr
        return sd_lwr, sd_upr

    if flat_id is None:
        if hole_id is None:
            prf = profile.profile(edge_extent, dir, steps)
        else:
            sd_lwr, sd_upr = process_edges(edge_extent, dir)
            prf = profile.profile((sd_upr, hole_id), dir, steps)
            prf += profile.profile((-hole_id, sd_lwr), dir, steps)
            
    else:
        sd_lwr, sd_upr = process_edges(edge_extent, dir)
        prf = []

        # compute top part of flat
        try:
            sag = profile.sag(0, flat_id)
        except TraceError:
            sag = None
        else:
            prf.append([sag, sd_lwr])

        if hole_id is None:
            prf += profile.profile((flat_id,), dir, steps)
        else:
            prf += profile.profile((flat_id, hole_id), dir, steps)
            prf += profile.profile((-hole_id, -flat_id), dir, steps)

        if sag is not None:
            prf.append([sag, sd_upr])

    if is_flipped:
        prf = [[-pt[0], pt[1]] for pt in prf]
    return prf


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
class Part(Protocol):
    """Abstract base class for all types of elements. """
    label_format: ClassVar[str]
    label: str
    parent: Any
    is_flipped: bool = False

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
        flat1, flat2: semi-diameter of flat or None. Setting to None will result in 
                      re-evaluation of flat ID
        do_flat1, do_flat2: 'if concave', 'always', 'never', 'if convex'
        handles: dict of graphical entities
        actions: dict of actions associated with the graphical handles
    """
    clut = rgbt.RGBTable(filename='red_blue64.csv',
                         data_range=[10.0, 100.])

    label_format = 'E{}'
    serial_number = 0

    def __init__(self, s1, s2, g, tfrm=None, idx=0, idx2=1, sd=1.,
                 label=None):
        if label is None:
            Element.serial_number += 1
            self.label = Element.label_format.format(Element.serial_number)
        else:
            self.label = label
        if tfrm is not None:
            self.tfrm = tfrm
        else:
            self.tfrm = (np.identity(3), np.array([0., 0., 0.]))
        self.s1 = s1
        self.profile1 = s1.profile
        self.s1_indx = idx
        self.s2 = s2
        self.s2_indx = idx2
        self.profile2 = s2.profile
        self.gap = g
        self.medium_name = self.gap.medium.name()
        self._sd = sd
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
        return attrs
        
    def __str__(self):
        fmt = 'Element: {!r}, {!r}, t={:.4f}, sd={:.4f}, glass: {}'
        return fmt.format(self.s1.profile, self.s2.profile, self.gap.thi,
                          self.sd, self.gap.medium.name())

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
        if not hasattr(self, 'do_flat1'):
            self.do_flat1 = 'if concave'
        if not hasattr(self, 'do_flat2'):
            self.do_flat2 = 'if concave'
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
        seq_model = self.parent.opt_model['seq_model']
        self.s1_indx = seq_model.ifcs.index(self.s1)
        self.s2_indx = seq_model.ifcs.index(self.s2)
        return [self.s1_indx, self.s2_indx]

    def reference_idx(self):
        return self.s1_indx

    def reference_interface(self):
        return self.s1

    def profile_list(self):
        return [self.profile1, self.profile2]

    def gap_list(self):
        return [self.gap]

    def get_bending(self):
        cv1 = self.s1.profile_cv
        cv2 = self.s2.profile_cv
        delta_cv = cv1 - cv2
        bending = 0.
        if delta_cv != 0.0:
            bending = (cv1 + cv2)/delta_cv
        return bending

    def set_bending(self, bending):
        cv1 = self.s1.profile_cv
        cv2 = self.s2.profile_cv
        delta_cv = cv1 - cv2
        cv2_new = 0.5*(bending - 1.)*delta_cv
        cv1_new = bending*delta_cv - cv2_new
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

    def compute_flat(self, s):
        ca = s.surface_od()
        if (1.0 - ca/self.sd) >= 0.05:
            flat = ca
        else:
            flat = None
        return flat

    def extent(self):
        if hasattr(self, 'edge_extent'):
            return self.edge_extent
        else:
            return (-self.sd, self.sd)

    def render_shape(self):
        def use_flat(do_flat, is_concave):
            if do_flat == 'always':
                return True
            elif do_flat == 'is concave' and is_concave:
                return True
            elif do_flat == 'is convex' and not is_concave:
                return True
            return False
        is_concave_s1 = self.s1.profile_cv < 0.0
        is_concave_s2 = self.s2.profile_cv > 0.0

        if use_flat(self.do_flat1, is_concave_s1):
            if self.flat1 is None:
                flat1 = self.flat1 = self.compute_flat(self.s1)
            else:
                flat1 = self.flat1
        else:
            flat1 = None
        poly = full_profile(self.profile1, self.is_flipped, self.extent(), 
                            flat1)

        if use_flat(self.do_flat2, is_concave_s2):
            if self.flat2 is None:
                flat2 = self.flat2 = self.compute_flat(self.s2)
            else:
                flat2 = self.flat2
        else:
            flat2 = None
        poly2 = full_profile(self.profile2, self.is_flipped, self.extent(),
                             flat2, dir=-1)

        for p in poly2:
            p[0] += self.gap.thi
        poly += poly2
        poly.append(poly[0])
        return poly

    def render_handles(self, opt_model):
        self.handles = {}
        thi = self.gap.thi

        shape = self.render_shape()
        color = calc_render_color_for_material(self.gap.medium)
        self.handles['shape'] = GraphicsHandle(shape, self.tfrm, 'polygon',
                                               color)

        extent = self.extent()
        if self.flat1 is not None:
            extent_s1 = self.flat1,
        else:
            extent_s1 = extent
        poly_s1 = full_profile(self.profile1, self.is_flipped, extent_s1)
        gh1 = GraphicsHandle(poly_s1, self.tfrm, 'polyline')
        self.handles['s1_profile'] = gh1

        if self.flat2 is not None:
            extent_s2 = self.flat2,
        else:
            extent_s2 = extent
        poly_s2 = full_profile(self.profile2, self.is_flipped, extent_s2, 
                               dir=-1)
        r, t = self.tfrm
        t_new = t + np.matmul(r, np.array([0, 0, thi]))
        gh2 = GraphicsHandle(poly_s2, (r, t_new), 'polyline')
        self.handles['s2_profile'] = gh2

        poly_sd_upr = []
        poly_sd_upr.append([poly_s1[-1][0], extent[1]])
        poly_sd_upr.append([poly_s2[0][0]+thi, extent[1]])
        self.handles['sd_upr'] = GraphicsHandle(poly_sd_upr, self.tfrm,
                                                'polyline')

        poly_sd_lwr = []
        poly_sd_lwr.append([poly_s2[-1][0]+thi, extent[0]])
        poly_sd_lwr.append([poly_s1[0][0], extent[0]])
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


class Mirror(Part):

    label_format = 'M{}'
    serial_number = 0

    def __init__(self, ifc, tfrm=None, idx=0, sd=1., thi=None, z_dir=1.0,
                 label=None):
        if label is None:
            Mirror.serial_number += 1
            self.label = Mirror.label_format.format(Mirror.serial_number)
        else:
            self.label = label

        self.render_color = (158, 158, 158, 64)
        if tfrm is not None:
            self.tfrm = tfrm
        else:
            self.tfrm = (np.identity(3), np.array([0., 0., 0.]))
        self.s = ifc
        self.s_indx = idx
        self.profile = ifc.profile
        self.z_dir = z_dir
        self.sd = sd
        self.flat = None
        self.thi = thi
        self.medium_name = 'Mirror'
        self.handles = {}
        self.actions = {}

    def get_thi(self):
        thi = self.thi
        if self.thi is None:
            thi = 0.05*self.sd
        return thi

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['parent']
        del attrs['s']
        del attrs['handles']
        del attrs['actions']
        encode_obj_reference(self, 'profile', attrs)
        return attrs

    def __str__(self):
        thi = self.get_thi()
        fmt = 'Mirror: {!r}, t={:.4f}, sd={:.4f}'
        return fmt.format(self.profile, thi, self.sd)

    def listobj_str(self):
        o_str = f"part: {type(self).__name__}, "
        o_str += self.profile.listobj_str()
        o_str += f"t={self.get_thi():.4f}, sd={self.sd:.4f}\n"
        return o_str

    def sync_to_restore(self, ele_model, surfs, gaps, tfrms, 
                        profile_dict, parts_dict):
        self.parent = ele_model
        if not hasattr(self, 'tfrm'):
            self.tfrm = tfrms[self.s_indx]
        self.s = surfs[self.s_indx]
        sync_obj_reference(self, 'profile', profile_dict, self.s.profile)
        if not hasattr(self, 'medium_name'):
            self.medium_name = 'Mirror'
        self.handles = {}
        self.actions = {}

    def sync_to_seq(self, seq_model):
        self.s_indx = seq_model.ifcs.index(self.s)
        self.profile = self.s.profile

    def tree(self, **kwargs):
        default_tag = '#element#mirror'
        tag = default_tag + kwargs.get('tag', '')
        # Interface branch
        m = Node('M', id=self, tag=tag)
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

    def substrate_offset(self):
        thi = self.get_thi()
        # We want to extend the mirror substrate along the same direction
        # of the incoming ray. The mirror's z_dir is following reflection so
        # flip the sign to get the preceding direction.
        offset = -self.z_dir*thi
        return offset

    def render_shape(self):
        poly = full_profile(self.profile, self.is_flipped, self.extent(),
                            self.flat)
        poly2 = full_profile(self.profile, self.is_flipped, self.extent(),
                             self.flat, dir=-1)

        offset = self.substrate_offset()

        for p in poly2:
            p[0] += offset
        poly += poly2
        poly.append(poly[0])
        return poly

    def render_handles(self, opt_model):
        self.handles = {}
        # ifcs_gbl_tfrms = opt_model.seq_model.gbl_tfrms

        self.handles['shape'] = GraphicsHandle(self.render_shape(), self.tfrm,
                                               'polygon', self.render_color)

        poly = full_profile(self.profile, self.is_flipped, self.extent())
        self.handles['s_profile'] = GraphicsHandle(poly, self.tfrm,
                                                   # ifcs_gbl_tfrms[self.s_indx],
                                                   'polyline')

        offset = self.substrate_offset()

        poly_sd_upr = []
        poly_sd_upr.append(poly[-1])
        poly_sd_upr.append([poly[-1][0]+offset, poly[-1][1]])
        self.handles['sd_upr'] = GraphicsHandle(poly_sd_upr, self.tfrm,
                                                'polyline')

        poly_sd_lwr = []
        poly_sd_lwr.append(poly[0])
        poly_sd_lwr.append([poly[0][0]+offset, poly[0][1]])
        self.handles['sd_lwr'] = GraphicsHandle(poly_sd_lwr, self.tfrm,
                                                'polyline')

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

    def __init__(self, ifc_list, label=None):
        if label is None:
            CementedElement.serial_number += 1
            self.label = CementedElement.label_format.format(
                CementedElement.serial_number)
        else:
            self.label = label

        g_tfrm = ifc_list[0][4]
        if g_tfrm is not None:
            self.tfrm = g_tfrm
        else:
            self.tfrm = (np.identity(3), np.array([0., 0., 0.]))
        self.idxs = []
        self.ifcs = []
        self.profiles = []
        self.gaps = []
        self.medium_name = ''
        for interface in ifc_list:
            i, ifc, g, z_dir, g_tfrm = interface
            self.idxs.append(i)
            self.ifcs.append(ifc)
            self.profiles.append(ifc.profile)
            if g is not None:
                self.gaps.append(g)
                if self.medium_name != '':
                    self.medium_name += ', '
                self.medium_name += g.medium.name()

        if len(self.gaps) == len(self.ifcs):
            self.gaps.pop()
            self.medium_name = self.medium_name.rpartition(',')[0]

        self._sd = self.update_size()
        self.flats = [None]*len(self.ifcs)

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
        del attrs['profile_polys']
        del attrs['handles']
        del attrs['actions']
        attrs['profile_ids'] = [str(id(p)) for p in self.profiles]
        del attrs['profiles']

        return attrs

    def __str__(self):
        fmt = 'CementedElement: {}'
        return fmt.format(self.idxs)

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
        if not hasattr(self, 'medium_name'):
            self.medium_name = self.gap.medium.name()
        self.handles = {}
        self.actions = {}

    def sync_to_seq(self, seq_model):
        # when updating, we want to use the stored object instances to get the
        # current indices into the interface list (e.g. to handle insertion and
        # deletion of interfaces)
        self.idxs = [seq_model.ifcs.index(ifc) for ifc in self.ifcs]
        self.profiles = [ifc.profile for ifc in self.ifcs]

    def tree(self, **kwargs):
        default_tag = '#element#cemented'
        tag = default_tag + kwargs.get('tag', '')
        zdir = kwargs.get('z_dir', 1)
        ce = Node(self.label, id=self, tag=tag)
        for i, sg in enumerate(itertools.zip_longest(self.ifcs, self.gaps)):
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

    def compute_flat(self, s):
        ca = s.surface_od()
        if (1.0 - ca/self.sd) >= 0.05:
            flat = ca
        else:
            flat = None
        return flat

    def extent(self):
        if hasattr(self, 'edge_extent'):
            return self.edge_extent
        else:
            return (-self.sd, self.sd)

    def render_shape(self):
        if self.profiles[0].cv < 0.0:
            self.flats[0] = self.compute_flat(self.ifcs[0])
        else:
            self.flats[0] = None
        if self.profiles[-1].cv > 0.0:
            self.flats[-1] = self.compute_flat(self.ifcs[-1])
        else:
            self.flats[-1] = None

        # generate the profile polylines
        self.profile_polys = []
        sense = 1
        for profile, flat in zip(self.profiles, self.flats):
            poly = full_profile(profile, self.is_flipped, self.extent(), 
                                flat, dir=sense)
            self.profile_polys.append(poly)
            sense = -sense

        # offset the profiles wrt the element origin
        thi = 0
        for i, poly_profile in enumerate(self.profile_polys[1:]):
            thi += self.gaps[i].thi
            for p in poly_profile:
                p[0] += thi

        # just return outline
        poly_shape = []
        poly_shape += self.profile_polys[0]
        poly_shape += self.profile_polys[-1]
        poly_shape.append(poly_shape[0])

        return poly_shape

    def render_handles(self, opt_model):
        self.handles = {}

        shape = self.render_shape()
        # self.handles['shape'] = GraphicsHandle(shape, self.tfrm, 'polygon')

        for i, gap in enumerate(self.gaps):
            poly = []
            poly += self.profile_polys[i]
            poly += self.profile_polys[i+1]
            poly.append(self.profile_polys[i][0])
            color = calc_render_color_for_material(gap.medium)
            self.handles['shape'+str(i+1)] = GraphicsHandle(poly, self.tfrm,
                                                            'polygon', color)
        return self.handles

    def handle_actions(self):
        self.actions = {}

        return self.actions


class ThinElement(Part):

    label_format = 'TL{}'
    serial_number = 0

    def __init__(self, ifc, tfrm=None, idx=0, sd=None, label=None):
        if label is None:
            ThinElement.serial_number += 1
            self.label = ThinElement.label_format.format(
                ThinElement.serial_number)
        else:
            self.label = label

        self.render_color = (192, 192, 192)
        if tfrm is not None:
            self.tfrm = tfrm
        else:
            self.tfrm = (np.identity(3), np.array([0., 0., 0.]))
        self.intrfc = ifc
        self.intrfc_indx = idx
        self.medium_name = 'Thin Element'
        if sd is not None:
            self.sd = sd
        else:
            self.sd = ifc.max_aperture
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
        self.handles = {}
        self.actions = {}

        ro_version = ele_model.opt_model.ro_version
        if version.parse(ro_version) < version.parse("0.7.0a"):
            ThinElement.serial_number += 1
            self.label = ThinElement.label_format.format(ThinElement.serial_number)

    def sync_to_seq(self, seq_model):
        self.intrfc_indx = seq_model.ifcs.index(self.intrfc)

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
        self.handles['shape'] = GraphicsHandle(shape, self.tfrm, 'polygon',
                                               self.render_color)
        return self.handles

    def handle_actions(self):
        self.actions = {}
        return self.actions

class DummyInterface(Part):

    label_format = 'D{}'
    serial_number = 0

    def __init__(self, ifc, idx=0, sd=None, tfrm=None, label=None):
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
        self.ref_ifc = ifc
        self.idx = idx
        self.profile = ifc.profile
        self.medium_name = 'Interface'
        if sd is not None:
            self.sd = sd
        else:
            self.sd = ifc.max_aperture
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
        self.handles = {}
        self.actions = {}

    def sync_to_seq(self, seq_model):
        self.idx = seq_model.ifcs.index(self.ref_ifc)
        self.profile = self.ref_ifc.profile

    def tree(self, **kwargs):
        default_tag = '#dummyifc'
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
        poly = full_profile(self.profile, self.is_flipped, (-self.sd, self.sd))
        return poly

    def render_handles(self, opt_model):
        self.handles = {}

        self.handles['shape'] = GraphicsHandle(self.render_shape(), self.tfrm,
                                               'polyline')

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


class AirGap(Part):

    label_format = 'AG{}'
    serial_number = 0

    def __init__(self, g, idx=0, tfrm=None, label=None, z_dir=1, **kwargs):
        if label is None:
            AirGap.serial_number += 1
            self.label = AirGap.label_format.format(AirGap.serial_number)
        else:
            self.label = label

        if tfrm is not None:
            self.tfrm = tfrm
        else:
            self.tfrm = (np.identity(3), np.array([0., 0., 0.]))

        self.render_color = (237, 243, 254, 64)  # light blue
        self.gap = g
        self.z_dir = z_dir
        self.medium_name = self.gap.medium.name()
        self.idx = idx
        self.handles = {}
        self.actions = {}

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['parent']
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
        if not hasattr(self, 'tfrm'):
            self.tfrm = tfrms[self.idx]
        if not hasattr(self, 'render_color'):
            self.render_color = (237, 243, 254, 64)  # light blue
        if not hasattr(self, 'medium_name'):
            self.medium_name = self.gap.medium.name()
        self.handles = {}
        self.actions = {}

        ro_version = ele_model.opt_model.ro_version
        if version.parse(ro_version) < version.parse("0.7.0a"):
            AirGap.serial_number += 1
            self.label = AirGap.label_format.format(AirGap.serial_number)

    def sync_to_seq(self, seq_model):
        self.idx = seq_model.gaps.index(self.gap)
        self.z_dir = seq_model.z_dir[self.idx]

    def tree(self, **kwargs):
        default_tag = '#airgap'
        tag = default_tag + kwargs.get('tag', '')
        ag = Node(self.label, id=self, tag=tag)
        t = Node('t', id=self.gap, tag='#thic', parent=ag)
        zdir = kwargs.get('z_dir', self.z_dir)
        Node(f'g{self.idx}', id=(self.gap, zdir), tag='#gap', parent=t)
        return ag

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
        pass

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

    def handle_actions(self):
        self.actions = {}

        ct_action = {}
        ct_action['x'] = AttrAction(self.gap, 'thi')
        self.actions['ct'] = ct_action

        return self.actions


class Assembly(Part):

    label_format = 'ASM{}'
    serial_number = 0

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
        self.handles = {}
        self.actions = {}

    def sync_to_seq(self, seq_model):
        self.tfrm = seq_model.gbl_tfrms[self.reference_idx()]

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
        opt_model: the :class:`~rayoptics.optical.opticalmodel.OpticalModel`
        elements: list of element type things

    """

    def __init__(self, opt_model, **kwargs):
        self.opt_model = opt_model
        self.elements: List[Part] = []

    def reset(self):
        self.__init__(self.opt_model)

    def __json_encode__(self):
        attrs = dict(vars(self))
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
                self.elements = [e for e in opt_model.parts_dict.values()]
                parts_dict = opt_model.parts_dict

        seq_model = opt_model.seq_model
        surfs = seq_model.ifcs
        gaps = seq_model.gaps
        tfrms = seq_model.compute_global_coords(1)

        self.reset_serial_numbers()
        for i, e in enumerate(self.elements, start=1):
            e.sync_to_restore(self, surfs, gaps, tfrms, 
                              profile_dict, parts_dict)
            if not hasattr(e, 'label'):
                e.label = e.label_format.format(i)
        self.sequence_elements()
        # self.relabel_airgaps()

    def reset_serial_numbers(self):
        Element.serial_number = 0
        Mirror.serial_number = 0
        CementedElement.serial_number = 0
        ThinElement.serial_number = 0
        DummyInterface.serial_number = 0
        AirGap.serial_number = 0
        Assembly.serial_number = 0

    def airgaps_from_sequence(self, seq_model, tfrms):
        """ add airgaps and dummy interfaces to an older version model """
        for e in self.elements:
            if isinstance(e, AirGap):
                return  # found an AirGap, model probably OK

        num_elements = 0
        seq_model = self.opt_model.seq_model
        for i, g in enumerate(seq_model.gaps):
            if g.medium.name().lower() == 'air':
                if i > 0:
                    s = seq_model.ifcs[i]
                    tfrm = tfrms[i]
                    num_elements = self.process_airgap(
                        seq_model, i, g, s, tfrm,
                        num_elements, add_ele=False)

    def add_dummy_interface_at_image(self, seq_model, tfrms):
        if len(self.elements) and self.elements[-1].label == 'Image':
            return

        s = seq_model.ifcs[-1]
        idx = seq_model.get_num_surfaces() - 1
        di = DummyInterface(s, sd=s.surface_od(), tfrm=tfrms[-1], idx=idx,
                            label='Image')
        self.opt_model.part_tree.add_element_to_tree(di, tag='#image')
        self.add_element(di)

    def update_model(self, **kwargs):
        # dynamically build element list from part_tree
        part_tree = self.opt_model['part_tree']
        part_tag = '#element#airgap#dummyifc#assembly'
        nodes = part_tree.nodes_with_tag(tag=part_tag)
        elements = [n.id for n in nodes]

        # hook or unhook elements from ele_model
        cur_set = set(self.elements)
        new_set = set(elements)
        added_ele = list(new_set.difference(cur_set))
        for e in added_ele:
            e.parent = self
        removed_ele = list(cur_set.difference(new_set))
        for e in removed_ele:
            e.parent = None

        self.elements = elements
        self.sequence_elements()

        src_model = kwargs.get('src_model', None)
        if src_model is not self:
            self.sync_to_seq(self.opt_model['seq_model'])

    def sync_to_seq(self, seq_model):
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
        self.elements.sort(key=lambda e: e.reference_idx())

        # Make sure z_dir matches the sequential model. Used to get
        # the correct substrate offset.
        if hasattr(seq_model, 'z_dir'):
            for e in self.elements:
                if hasattr(e, 'z_dir'):
                    e.z_dir = seq_model.z_dir[e.reference_idx()]

    def relabel_airgaps(self):
        for i, e in enumerate(self.elements):
            if isinstance(e, AirGap):
                eb = self.elements[i-1].label
                ea = self.elements[i+1].label
                e.label = AirGap.label_format.format(eb + '-' + ea)

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
