#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2017 Michael J. Hayford
""" Module for optical surface related classes

    Surface
        Container of profile, extent, position and orientation information of
        the surface

    DecenterData
        Maintains data and actions to support 4 types of position and
        orientation changes.

        - DEC: pos and orientation applied prior to surface
        - REV: pos and orientation applied following surface in reverse
        - DAR: pos and orientation applied prior to surface and then returned to initial frame
        - BEN: used for fold mirrors, orientation applied before and after surface

    Aperture
        - Circular
        - Rectangular
        - Elliptical

.. Created on Sat Sep 16 09:22:05 2017

.. codeauthor: Michael J. Hayford
"""

from enum import Enum, auto
from math import sqrt
import numpy as np

from rayoptics.seq import interface
from . import profiles
from rayoptics.optical.model_enums import get_decenter_for_type
from rayoptics.raytr.traceerror import TraceError
from rayoptics.util import misc_math


class InteractionMode(Enum):
    """ enum for different interact_mode specifications

    Retained to restore old files

    .. deprecated:: 0.4.5
    """
    Transmit = auto()  #: propagate in transmission at this interface
    Reflect = auto()   #: propagate in reflection at this interface


class Surface(interface.Interface):
    """ Container of profile, extent, position and orientation. 

    Attributes:
        label: optional label
        profile: :class:`~.elem.profiles.SurfaceProfile`
        clear_apertures: list of :class:`Aperture`
        edge_apertures: list of :class:`Aperture`
        """

    def __init__(self, lbl='', profile=None,
                 clear_apertures=None, edge_apertures=None,
                 **kwargs):
        super().__init__(**kwargs)
        self.label = lbl
        if profile:
            self.profile = profile
        else:
            self.profile = profiles.Spherical()
        self.clear_apertures = clear_apertures if clear_apertures else []
        self.edge_apertures = edge_apertures if edge_apertures else []

    def __repr__(self):
        if len(self.label) > 0:
            return "{!s}(lbl={!r}, profile={!r}, interact_mode='{!s}')" \
                   .format(type(self).__name__,
                           self.label, self.profile, self.interact_mode)
        else:
            return "{!s}(profile={!r}, interact_mode='{!s}')" \
                   .format(type(self).__name__,
                           self.profile, self.interact_mode)

    def interface_type(self):
        return type(self.profile).__name__

    def listobj_str(self):
        o_str = f"{self.label}: " if self.label != "" else ""
        o_str += f"{self.interact_mode}\n"
        if self.profile is not None:
            o_str += self.profile.listobj_str()
        if hasattr(self, 'phase_element') and self.phase_element is not None:
            o_str += self.phase_element.listobj_str()
        if hasattr(self, 'decenter') and self.decenter is not None:
            o_str += self.decenter.listobj_str()

        o_str += f"surface_od={self.surface_od()}\n"
        for ca in self.clear_apertures:
            o_str += ca.listobj_str()
        for ea in self.edge_apertures:
            o_str += ea.listobj_str()

        return o_str

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['profile']
        attrs['profile_id'] = str(id(self.profile))
        return attrs

    def update(self):
        super().update()
        self.profile.update()

    def sync_to_restore(self, opt_model):
        super().sync_to_restore(opt_model)
        if not hasattr(self, 'profile'):
            self.profile = opt_model.profile_dict[self.profile_id]
            delattr(self, 'profile_id')
        for ca in self.clear_apertures:
            ca.sync_to_restore(opt_model)
        for ea in self.edge_apertures:
            ea.sync_to_restore(opt_model)

    @property
    def profile_cv(self):
        return self.profile.cv

    @profile_cv.setter
    def profile_cv(self, cv):
        self.profile.cv = cv

    @property
    def optical_power(self):
        return self.delta_n * self.profile.cv

    @optical_power.setter
    def optical_power(self, pwr):
        self.profile.cv = pwr/self.delta_n if self.delta_n != 0.0 else 0.0

    def set_optical_power(self, pwr, n_before, n_after):
        self.delta_n = n_after - n_before
        self.optical_power = pwr

    def apply_scale_factor(self, scale_factor):
        super().apply_scale_factor(scale_factor)
        self.profile.apply_scale_factor(scale_factor)
        abs_scale_factor = abs(scale_factor)
        for e in self.edge_apertures:
            e.apply_scale_factor(abs_scale_factor)
        for ca in self.clear_apertures:
            ca.apply_scale_factor(abs_scale_factor)

    def from_first_order(self, nu_before, nu_after, y):
        pass

    def update_following_reflection(self):
        self.apply_scale_factor(-1)

    def z_sag(self, pt):
        return self.profile.sag(0., pt[1])

    def set_z_sag(self, pt):
        self.profile.cv = self.calc_cv_from_zsag(pt)

    def calc_cv_from_zsag(self, pt):
        x, y = pt
        cv = 2*x / (x**2 + y**2)
        return cv

    def flip(self):
        self.profile.flip()

    def set_max_aperture(self, max_ap):
        super().set_max_aperture(max_ap)
        for ap in self.clear_apertures:
            # print(f"max ap: {max_ap}")
            ap.set_dimension(max_ap, max_ap)

    def surface_od(self):
        od = 0
        if len(self.edge_apertures) > 0:
            for e in self.edge_apertures:
                edg = e.max_dimension()
                if edg > od:
                    od = edg
        elif len(self.clear_apertures) > 0:
            for ca in self.clear_apertures:
                ap = ca.max_dimension()
                if ap > od:
                    od = ap
        else:
            od = self.max_aperture

        return od

    def point_inside(self, x: float, y: float, fuzz: float = 1e-5) -> bool:
        is_inside = True
        if len(self.clear_apertures) > 0:
            for ca in self.clear_apertures:
                is_inside = is_inside and ca.point_inside(x, y, fuzz)
                if not is_inside:
                    return is_inside
        else:
            return super().point_inside(x, y, fuzz)

        return is_inside

    def edge_pt_target(self, rel_dir):
        """ Get a target for ray aiming to aperture boundaries.
        
        """
        if len(self.clear_apertures) > 0:
            # hardwire to 1 aperture until the need for more arises
            return self.clear_apertures[0].edge_pt_target(rel_dir)
        else:
            return super().edge_pt_target(rel_dir)

    def get_y_aperture_extent(self):
        """ returns [y_min, y_max] for the union of apertures """
        od = [1.0e10, -1.0e10]
        if len(self.edge_apertures) > 0:
            for e in self.edge_apertures:
                edg = e.bounding_box()
                if edg[0][1] < od[0]:
                    od[0] = edg[0][1]
                if edg[1][1] > od[1]:
                    od[1] = edg[1][1]
        elif len(self.clear_apertures) > 0:
            for ca in self.clear_apertures:
                ap = ca.bounding_box()
                if ap[0][1] < od[0]:
                    od[0] = ap[0][1]
                if ap[1][1] > od[1]:
                    od[1] = ap[1][1]
        else:
            od = [-self.max_aperture, self.max_aperture]

        return od

    def full_profile(self, edge_extent, flat_id=None, dir=1, steps=6):
        if flat_id is None:
            return self.profile.profile(edge_extent, dir, steps)
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
                sag = self.profile.sag(0, flat_id)
            except TraceError:
                sag = None
            else:
                prf.append([sag, sd_lwr])
            prf += self.profile.profile((flat_id,), dir, steps)
            if sag is not None:
                prf.append([sag, sd_upr])
            return prf

    def intersect(self, p0, d, eps=1.0e-12, z_dir=1.0):
        return self.profile.intersect(p0, d, eps, z_dir)

    def normal(self, p):
        return self.profile.normal(p)


class DecenterData():
    """ Maintains data and actions for position and orientation changes.

        - 'decenter': pos and orientation applied prior to surface
        - 'reverse': pos and orientation applied following surface in reverse
        - 'dec and return': pos and orientation applied prior to surface and then returned to initial frame
        - 'bend':  used for fold mirrors, orientation applied before and after surface

    """

    def __init__(self, dtype, x=0., y=0., alpha=0., beta=0., gamma=0.):
        self.dtype = dtype
        # x, y, z vertex decenter
        self.dec = np.array([x, y, 0.])
        # alpha, beta, gamma euler angles
        self.euler = np.array([alpha, beta, gamma])
        # x, y, z rotation point offset
        self.rot_pt = np.array([0., 0., 0.])
        self.rot_mat = None

    def __json_decode__(self, **attrs):
        for a_key, a_val in attrs.items():
            if a_key == 'dtype':
                self._dtype = (a_val if isinstance(a_val, str)
                               else get_decenter_for_type(a_val))
            else:
                setattr(self, a_key, a_val)

    def __repr__(self):
        return "%r: Decenter: %r, Tilt: %r" % (self.dtype, self.dec,
                                               self.euler)

    def listobj_str(self):
        o_str = f"decenter type: {self.dtype}\n"
        o_str += f"decenter: {self.dec}\n"
        o_str += f"euler angles: {self.euler}\n"
        return o_str

    @property
    def dtype(self):
        return self._dtype

    @dtype.setter
    def dtype(self, value):
        self._dtype = (value if isinstance(value, str)
                       else get_decenter_for_type(value))

    def update(self):
        if self.euler.any():
            self.rot_mat = misc_math.euler2rot3d(self.euler)
        else:
            self.rot_mat = None

    def apply_scale_factor(self, scale_factor):
        self.dec *= scale_factor
        self.rot_pt *= scale_factor

    def tform_before_surf(self):
        if self.dtype != 'reverse':
            return self.rot_mat, self.dec
        else:
            return None, np.array([0., 0., 0.])

    def tform_after_surf(self):
        if self.dtype == 'reverse' or self.dtype == 'dec and return':
            rt = self.rot_mat
            if self.rot_mat is not None:
                rt = self.rot_mat.transpose()
            return rt, -self.dec
        elif self.dtype == 'bend':
            return self.rot_mat, np.array([0., 0., 0.])
        else:
            return None, np.array([0., 0., 0.])


class Aperture():
    def __init__(self, x_offset=0.0, y_offset=0.0, rotation=0.0):
        self.x_offset = x_offset
        self.y_offset = y_offset
        self.rotation = rotation

    def sync_to_restore(self, opt_model):
        if not hasattr(self, 'x_offset'):
            self.x_offset = 0.0
        if not hasattr(self, 'y_offset'):
            self.y_offset = 0.0
        if not hasattr(self, 'rotation'):
            self.rotation = 0.0

    def listobj_str(self):
        o_str = ""
        if self.x_offset != 0. or self.x_offset != 0. or self.rotation != 0.:
            o_str = (f"x_offset={self.x_offset}   y_offset={self.y_offset}"
                     f"   rotation={self.rotation}\n")
        return o_str

    def dimension(self):
        pass

    def set_dimension(self, x, y):
        pass

    def max_dimension(self):
        x, y = self.dimension()
        return sqrt(x*x + y*y)

    def point_inside(self, x: float, y: float, fuzz: float = 1e-5) -> bool:
        pass

    def bounding_box(self):
        center = np.array([self.x_offset, self.y_offset])
        extent = np.array(self.dimension())
        return center-extent, center+extent

    def apply_scale_factor(self, scale_factor):
        self.x_offset *= scale_factor
        self.y_offset *= scale_factor

    def tform(self, x, y):
        x -= self.x_offset
        y -= self.y_offset
        return x, y


class Circular(Aperture):
    def __init__(self, radius=1.0, **kwargs):
        super().__init__(**kwargs)
        self.radius = radius

    def listobj_str(self):
        o_str = f"ca: radius={self.radius}\n"
        o_str += super().listobj_str()
        return o_str

    def dimension(self):
        return (self.radius, self.radius)

    def set_dimension(self, x, y):
        self.radius = x

    def max_dimension(self):
        return self.radius

    def point_inside(self, x: float, y: float, fuzz: float = 1e-5) -> bool:
        x, y = self.tform(x, y)
        return sqrt(x*x + y*y) <= self.radius + fuzz

    def edge_pt_target(self, rel_dir):
        """ Get a target for ray aiming to aperture boundaries.
        
        """
        edge_pt = np.array([self.radius*rel_dir[0], 
                            self.radius*rel_dir[1]])
        return edge_pt

    def apply_scale_factor(self, scale_factor):
        super().apply_scale_factor(scale_factor)
        self.radius *= scale_factor


class Rectangular(Aperture):
    def __init__(self, x_half_width=1.0, y_half_width=1.0, **kwargs):
        super().__init__(**kwargs)
        self.x_half_width = x_half_width
        self.y_half_width = y_half_width

    def listobj_str(self):
        o_str = (f"ca: {type(self).__name__}: x_half_width={self.x_half_width}"
                 f"   y_half_width={self.y_half_width}\n")
        o_str += super().listobj_str()
        return o_str

    def dimension(self):
        return (self.x_half_width, self.y_half_width)

    def set_dimension(self, x, y):
        self.x_half_width = abs(x)
        self.y_half_width = abs(y)

    def point_inside(self, x: float, y: float, fuzz: float = 1e-5) -> bool:
        x, y = self.tform(x, y)
        return (abs(x) <= self.x_half_width + fuzz 
                and abs(y) <= self.y_half_width + fuzz)

    def edge_pt_target(self, rel_dir):
        """ Get a target for ray aiming to aperture boundaries. """
        edge_pt = np.array([self.x_half_width*rel_dir[0], 
                            self.y_half_width*rel_dir[1]])
        return edge_pt

    def apply_scale_factor(self, scale_factor):
        super().apply_scale_factor(scale_factor)
        self.x_half_width *= scale_factor
        self.y_half_width *= scale_factor


class Elliptical(Aperture):
    def __init__(self, x_half_width=1.0, y_half_width=1.0, **kwargs):
        super().__init__(**kwargs)
        self.x_half_width = x_half_width
        self.y_half_width = y_half_width

    def listobj_str(self):
        o_str = (f"ca: {type(self).__name__}: x_half_width={self.x_half_width}"
                 f"   y_half_width={self.y_half_width}\n")
        o_str += super().listobj_str()
        return o_str

    def dimension(self):
        return (self.x_half_width, self.y_half_width)

    def set_dimension(self, x, y):
        self.x_half_width = abs(x)
        self.y_half_width = abs(y)

    def apply_scale_factor(self, scale_factor):
        super().apply_scale_factor(scale_factor)
        self.x_half_width *= scale_factor
        self.y_half_width *= scale_factor
