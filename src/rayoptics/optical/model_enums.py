#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" **DEPRECATED**: optical model enums

    The enums in this module are deprecated in favor of strings conveying the 
    same information. The functions in this module are used to convert enums 
    into the corresponding strings.

.. Created on Tue Dec  4 11:32:57 2018

.. codeauthor: Michael J. Hayford
"""

from enum import Enum


class DimensionType(Enum):
    """ **DEPRECATED**: enum for different linear dimensions """
    MM = 0  #: millimeters
    CM = 1  #: centimeters
    M = 2   #: meters
    IN = 3  #: inches
    FT = 4  #: feet


def get_dimension_for_type(dimension_type):
    if dimension_type == DimensionType.MM:
        dimension_key = 'mm'
    elif dimension_type == DimensionType.CM:
        dimension_key = 'cm'
    elif dimension_type == DimensionType.M:
        dimension_key = 'meters'
    elif dimension_type == DimensionType.IN:
        dimension_key = 'inches'
    elif dimension_type == DimensionType.FT:
        dimension_key = 'feet'
    return dimension_key


class DecenterType(Enum):
    """ **DEPRECATED**: enum for different tilt and decenter types """
    LOCAL = 0  #: pos and orientation applied prior to surface
    REV = 1    #: pos and orientation applied following surface in reverse
    DAR = 2    #: pos and orientation applied prior to surface and then returned to initial frame
    BEND = 3   #: used for fold mirrors, orientation applied before and after surface


def get_decenter_for_type(decenter_type):
    if decenter_type == DecenterType.LOCAL:
        decenter_key = 'decenter'
    elif decenter_type == DecenterType.REV:
        decenter_key = 'reverse'
    elif decenter_type == DecenterType.DAR:
        decenter_key = 'dec and return'
    elif decenter_type == DecenterType.BEND:
        decenter_key = 'bend'
    return decenter_key
