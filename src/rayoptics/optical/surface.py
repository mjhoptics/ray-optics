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

from . import interface
from . import profiles
import transforms3d as t3d
from .model_enums import DecenterType as dec
from rayoptics.optical.traceerror import TraceError


class InteractionMode(Enum):
    """ enum for different interact_mode specifications

    Retained to restore old files

    .. deprecated:: 0.4.5
    """
    Transmit = auto()  #: propagate in transmission at this interface
    Reflect = auto()   #: propagate in reflection at this interface


class Surface(interface.Interface):
    """ Container of profile, extent, position and orientation. """

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
            return "{!s}(lbl={!r}, profile={!r}, interact_mode={!s})" \
                   .format(type(self).__name__,
                           self.label, self.profile, self.interact_mode)
        else:
            return "{!s}(profile={!r}, interact_mode={!s})" \
                   .format(type(self).__name__,
                           self.profile, self.interact_mode)

    def interface_type(self):
        return type(self.profile).__name__

    def update(self):
        super().update()
        self.profile.update()

    def sync_to_restore(self, opt_model):
        super().sync_to_restore(opt_model)
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
        self.max_aperture *= scale_factor
        self.profile.apply_scale_factor(scale_factor)
        for e in self.edge_apertures:
            e.apply_scale_factor(scale_factor)
        for ca in self.clear_apertures:
            ca.apply_scale_factor(scale_factor)

    def from_first_order(self, nu_before, nu_after, y):
        pass

    def z_sag(self, pt):
        return self.profile.sag(0., pt[1])

    def set_z_sag(self, pt):
        self.profile.cv = self.calc_cv_from_zsag(pt)

    def calc_cv_from_zsag(self, pt):
        x, y = pt
        cv = 2*x / (x**2 + y**2)
        return cv

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

        - LOCAL: pos and orientation applied prior to surface
        - REV:   pos and orientation applied following surface in reverse
        - DAR:   pos and orientation applied prior to surface and then returned to initial frame
        - BEND:  used for fold mirrors, orientation applied before and after surface

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

    def __repr__(self):
        return "%r: Decenter: %r, Tilt: %r" % (self.dtype.name, self.dec,
                                               self.euler)

    def update(self):
        def convertl2r(self):
            return np.array([-self.euler[0], -self.euler[1], self.euler[2]])
        if self.euler.any():
            self.rot_mat = t3d.euler.euler2mat(*np.deg2rad(convertl2r(self)))
        else:
            self.rot_mat = None

    def apply_scale_factor(self, scale_factor):
        self.dec *= scale_factor
        self.rot_pt *= scale_factor

    def tform_before_surf(self):
        if self.dtype is not dec.REV:
            return self.rot_mat, self.dec
        else:
            return None, np.array([0., 0., 0.])

    def tform_after_surf(self):
        if self.dtype is dec.REV or self.dtype is dec.DAR:
            rt = self.rot_mat
            if self.rot_mat is not None:
                rt = self.rot_mat.transpose()
            return rt, -self.dec
        elif self.dtype is dec.BEND:
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

    def dimension(self):
        pass

    def set_dimension(self, x, y):
        pass

    def max_dimension(self):
        x, y = self.dimension()
        return sqrt(x*x + y*y)

    def bounding_box(self):
        center = np.array([self.x_offset, self.y_offset])
        extent = np.array(self.dimension())
        return center-extent, center+extent

    def apply_scale_factor(self, scale_factor):
        self.x_offset *= scale_factor
        self.y_offset *= scale_factor


class Circular(Aperture):
    def __init__(self, radius=1.0, **kwargs):
        super().__init__(**kwargs)
        self.radius = radius

    def dimension(self):
        return (self.radius, self.radius)

    def set_dimension(self, x, y):
        self.radius = sqrt(x*x + y*y)

    def max_dimension(self):
        return self.radius

    def apply_scale_factor(self, scale_factor):
        super().apply_scale_factor(scale_factor)
        self.radius *= scale_factor


class Rectangular(Aperture):
    def __init__(self, x_half_width=1.0, y_half_width=1.0, **kwargs):
        super().__init__(**kwargs)
        self.x_half_width = x_half_width
        self.y_half_width = y_half_width

    def dimension(self):
        return (self.x_half_width, self.y_half_width)

    def set_dimension(self, x, y):
        self.x_half_width = abs(x)
        self.y_half_width = abs(y)

    def apply_scale_factor(self, scale_factor):
        super().apply_scale_factor(scale_factor)
        self.x_half_width *= scale_factor
        self.y_half_width *= scale_factor


class Elliptical(Aperture):
    def __init__(self, x_half_width=1.0, y_half_width=1.0, **kwargs):
        super().__init__(**kwargs)
        self.x_half_width = x_half_width
        self.y_half_width = y_half_width

    def dimension(self):
        return (self.x_half_width, self.y_half_width)

    def set_dimension(self, x, y):
        self.x_half_width = abs(x)
        self.y_half_width = abs(y)

    def apply_scale_factor(self, scale_factor):
        super().apply_scale_factor(scale_factor)
        self.x_half_width *= scale_factor
        self.y_half_width *= scale_factor
