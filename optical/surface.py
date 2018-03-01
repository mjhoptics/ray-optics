#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 09:22:05 2017

@author: Mike
"""


from . import profiles
from math import sqrt
import numpy as np
import transforms3d as t3d


class Surface:
    def __init__(self, lbl=''):
        self.label = lbl
        self.refract_mode = ''
        self.profile = profiles.Spherical()
        self.decenter = None
        self.clear_apertures = []
        self.edge_apertures = []

    def __repr__(self):
        if len(self.label) > 0:
            return "Surface(%r: %r)" % (self.label, self.profile)
        else:
            return "Surface(%r)" % (self.profile)

    def full_profile(self, sd, flat_id=None, dir=1, steps=6):
        if flat_id is None:
            return self.profile.profile(sd, dir, steps)
        else:
            prf = []
            sag = self.profile.sag(0, flat_id)
            prf.append([sag, -dir*sd])
            prf += self.profile.profile(flat_id, dir, steps)
            prf.append([sag, dir*sd])
            return prf

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
        return od


class DecenterData():
    def __init__(self):
        self.type = 'DEC'
        # x, y, z vertex decenter
        self.dec = np.array([0, 0, 0])
        # alpha, beta, gamma euler angles
        self.euler = np.array([0, 0, 0])
        # x, y, z rotation point offset
        self.rot_pt = np.array([0, 0, 0])
        self.rot_mat = None

    def __repr__(self):
        return "%r: Decenter: %r, Tilt: %r" % (self.type, self.dec,
                                               self.euler)

    def update(self):
        def convertl2r(self):
            return np.array([-self.euler[0], -self.euler[1], self.euler[2]])
        if self.euler.any():
            self.rot_mat = t3d.euler.euler2mat(*np.deg2rad(convertl2r(self)))
        else:
            self.rot_mat = None

    def tform_before_surf(self):
        if self.type is not 'REV':
            return self.rot_mat, self.dec
        else:
            return None, np.array([0., 0., 0.])

    def tform_after_surf(self):
        if self.type is 'REV' or self.type is 'DAR':
            return self.rot_mat.transpose(), -self.dec
        elif self.type is 'BEN':
            return self.rot_mat, np.array([0., 0., 0.])
        else:
            return None, np.array([0., 0., 0.])


class Aperture():
    def __init__(self):
        self.type = ''
        self.x_offset = 0.0
        self.y_offset = 0.0
        self.rotation = 0.0

    def dimension(self):
        pass

    def set_dimension(self, x, y):
        pass

    def max_dimension(self):
        x, y = self.dimension()
        return sqrt(x*x + y*y)


class Circular(Aperture):
    def __init__(self, r_=1.0):
        self.type = 'Circular'
        self.radius = r_

    def dimension(self):
        return (self.radius, self.radius)

    def set_dimension(self, x, y):
        self.radius = sqrt(x*x + y*y)

    def max_dimension(self):
        return self.radius


class Rectangular(Aperture):
    def __init__(self, x_=1.0, y_=1.0):
        self.type = 'Rectangular'
        self.x_half_width = x_
        self.y_half_width = y_

    def dimension(self):
        return (self.x_half_width, self.y_half_width)

    def set_dimension(self, x, y):
        self.x_half_width = abs(x)
        self.y_half_width = abs(y)


class Elliptical(Aperture):
    def __init__(self, x_=1.0, y_=1.0):
        self.type = 'Elliptical'
        self.x_half_width = x_
        self.y_half_width = y_

    def dimension(self):
        return (self.x_half_width, self.y_half_width)

    def set_dimension(self, x, y):
        self.x_half_width = abs(x)
        self.y_half_width = abs(y)
