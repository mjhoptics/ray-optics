#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" type hints for rayoptics

Created on Tue Dec  7 10:48:54 2021

.. codeauthor: Michael J. Hayford
"""
from typing import Iterator
from typing import Literal
from typing import Optional
from rayoptics.coord_geometry_types import Vec3d, Dir3d, Tfm3d

DCenterTypes = Literal['decenter', 'reverse', 'dec and return', 'bend']
InteractMode = Literal['transmit', 'reflect', 'dummy', 'phantom']

Z_DIR = Literal[-1, 1]

PathSeg = tuple['Interface', 'Gap', Tfm3d, float, Z_DIR]
SeqPath = list[PathSeg]
Path = Iterator[PathSeg]

RaySeg = tuple[Vec3d, Dir3d, float, Dir3d]
RayPkg = tuple[list[RaySeg], float, float]
#OptionalTraceError = 'TraceError' | None
RayResult = tuple[RayPkg | None, Optional['TraceError']]
