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
    EPD = 0  #: entrance pupil diameter
    NAO = 1  #: object space numerical aperture
    FNO = 2  #: image space f/#
    NA = 3   #: image space numerical aperture


class FieldType(Enum):
    """ enum for different field specifications """
    OBJ_ANG = 0  #: object space angle in degrees
    OBJ_HT = 1   #: object height
    IMG_HT = 2   #: image height
    IMG_ANG = 3   #: image space angle in degrees


class DimensionType(Enum):
    """ enum for different linear dimensions """
    MM = 0  #: millimeters
    CM = 1  #: centimeters
    M = 2   #: meters
    IN = 3  #: inches
    FT = 4  #: feet


class DecenterType(Enum):
    """ enum for different tilt and decenter types """
    LOCAL = 0  #: pos and orientation applied prior to surface
    REV = 1    #: pos and orientation applied following surface in reverse
    DAR = 2    #: pos and orientation applied prior to surface and then returned to initial frame
    BEND = 3   #: used for fold mirrors, orientation applied before and after surface
