#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" optical model enums
Created on Tue Dec  4 11:32:57 2018

@author: Michael J. Hayford
"""

from enum import Enum


class PupilType(Enum):
    FNO = 0
    EPD = 1
    NA = 2
    NAO = 3


class FieldType(Enum):
    OBJ_ANG = 0
    OBJ_HT = 1
    IMG_HT = 2


class DimensionType(Enum):
    MM = 0
    CM = 1
    M = 2
    IN = 3
    FT = 4
