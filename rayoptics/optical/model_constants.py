#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" optical model constants
Created on Wed May 23 16:00:55 2018

@author: Michael J. Hayford
"""

# sequential model data structures
Surf, Gap, Indx = range(3)


# paraxial optics data structures
# 4 parts of paraxial model: axial ray, principal ray,
#                            lens data, optical invariant
ax, pr, lns, inv = range(4)
# paraxial ray data at an interface: height, n*slope, n*angle of incidence
ht, slp, aoi = range(3)
# lens data: power, reduced distance, refractive index (n), refract mode
pwr, tau, indx, rmd = range(4)
