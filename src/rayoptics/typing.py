#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" type hints for rayoptics

Created on Tue Dec  7 10:48:54 2021

.. codeauthor: Michael J. Hayford
"""
from typing import Iterator
from typing import Literal
import numpy.typing as npt

Vec2d = npt.NDArray
Dir2d = npt.NDArray
Mat2d = npt.NDArray
Tfm2d = tuple[Mat2d|None, Vec2d]

Vec3d = npt.NDArray
Dir3d = npt.NDArray
Mat3d = npt.NDArray
Tfm3d = tuple[Mat3d|None, Vec3d]

Vec4d = npt.NDArray
Dir4d = npt.NDArray
Mat4d = npt.NDArray
Tfm2d = tuple[Mat4d|None, Vec4d]

V2d = npt.ArrayLike
M2d = npt.ArrayLike
T2d = tuple[M2d|None, V2d]

V3d = npt.ArrayLike
D3d = npt.ArrayLike
M3d = npt.ArrayLike
T3d = tuple[M3d|None, V3d]

V4d = npt.ArrayLike
M4d = npt.ArrayLike
T4d = tuple[M4d|None, V4d]

DCenterTypes = Literal['decenter', 'reverse', 'dec and return', 'bend']
InteractMode = Literal['transmit', 'reflect', 'dummy']

Z_DIR = Literal[-1, 1]

PathSeg = tuple["Interface", "Gap", Tfm3d, float, Z_DIR]
SeqPath = list[PathSeg]
Path = Iterator[PathSeg]

RaySeg = tuple[Vec3d, Dir3d, float, Dir3d]
RayPkg = tuple[list[RaySeg], float, float]
RayResult = tuple[RayPkg|None, "TraceError"|None]
