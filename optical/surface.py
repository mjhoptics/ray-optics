#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 09:22:05 2017

@author: Mike
"""


import profiles


class Surface:
    def __init__(self, lbl=''):
        self.label = lbl
        self.refract_mode = ''
        self.profile = profiles.Spherical()
        self.decenter = DecenterData()
        self.clear_apertures = {}
        self.edge_apertures = {}

    def __repr__(self):
        if len(self.label) > 0:
            return "Surface(%r: %r)" % (self.label, self.profile)
        else:
            return "Surface(%r)" % (self.profile)


class DecenterData():
    def __init__(self):
        self.type = 'DEC'
        self.x_decenter = 0.0
        self.y_decenter = 0.0
        self.z_decenter = 0.0
        self.alpha = 0.0
        self.beta = 0.0
        self.gamma = 0.0
        self.x_offset = 0.0
        self.y_offset = 0.0
        self.z_offset = 0.0


class Aperture():
    def __init__(self):
        self.type = ''
        self.x_offset = 0.0
        self.y_offset = 0.0
        self.rotation = 0.0


class Circular(Aperture):
    def __init__(self):
        self.type = 'Circular'
        self.radius = 1.0


class Rectangular(Aperture):
    def __init__(self):
        self.type = 'Rectangular'
        self.x_half_width = 1.0
        self.y_half_width = 1.0


class Elliptical(Aperture):
    def __init__(self):
        self.type = 'Elliptical'
        self.x_half_width = 1.0
        self.y_half_width = 1.0
