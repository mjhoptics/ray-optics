#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Top level model classes

Created on Wed Mar 14 11:08:28 2018

@author: Michael J. Hayford
"""

from . import sequential as seq
from . import elements as ele


class SystemSpec:
    dims = ('M', 'CM', 'MM', 'IN', 'FT')

    def __init__(self):
        self.title = ''
        self.initials = ''
        self.dimensions = 'MM'
        self.aperture_override = ''
        self.temperature = 20.0
        self.pressure = 760.0


class OpticalModel:
    """ Top level container for optical model. """
    def __init__(self):
        self.radius_mode = False
        self.system_spec = SystemSpec()
        self.seq_model = seq.SequentialModel(self)
        self.ele_model = ele.ElementModel(self)

    def reset(self):
        rdm = self.radius_mode
        self.__init__()
        self.radius_mode = rdm

    def update_model(self):
        self.seq_model.update_model()
