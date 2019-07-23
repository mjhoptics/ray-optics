#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" optical model constants

.. Created on Wed May 23 16:00:55 2018

.. codeauthor: Michael J. Hayford
"""

# sequential ray trace data structures
Intfc, Gap, Tfrm, Indx, Zdir = range(5)

Surf, Gap = range(2)

# paraxial optics data structures
# 4 parts of paraxial model: axial ray, principal ray,
#                            lens data, optical invariant
ax, pr, lns, inv = range(4)
# paraxial ray data at an interface: height, n*slope, n*angle of incidence
ht, slp, aoi = range(3)
# lens data: power, reduced distance, refractive index (n), refract mode
pwr, tau, indx, rmd = range(4)

# ray trace segment data
# p: intersection point with interface
# d: direction cosine exiting the interface
# dst: distance from intersection point to next interface
# nrml: surface normal at intersection point
# phase: optical phase introduced at the interface
p, d, dst, nrml, phase = range(5)

# ray package
# ray: list of ray segment data
# op:  optical path wrt equally inclined chords to the optical axis
# wvl: wavelength (in nm) that the ray was traced in
ray, op, wvl = range(3)
