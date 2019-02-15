#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" optical model enums

.. Created on Tue Dec  4 11:32:57 2018

.. codeauthor: Michael J. Hayford
"""

from enum import Enum


class PupilType(Enum):
    """ enum for different aperture specifications """
    FNO = 0  #: image space f/#
    EPD = 1  #: entrance pupil diameter
    NA = 2   #: image space numerical aperture
    NAO = 3  #: object space numerical aperture


class FieldType(Enum):
    """ enum for different field specifications """
    OBJ_ANG = 0  #: object space angle in degrees
    OBJ_HT = 1   #: object height
    IMG_HT = 2   #: image height


class DimensionType(Enum):
    """ enum for different linear dimensions """
    MM = 0  #: millimeters
    CM = 1  #: centimeters
    M = 2   #: meters
    IN = 3  #: inches
    FT = 4  #: feet
