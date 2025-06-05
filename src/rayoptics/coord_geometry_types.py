#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" type hints for for vectors and matrices

These type hints are provided to ensure a consist convention of distinguishing between numpy arrays and array-like arrays.

2D, 3D and 4D cases are provided.

VecNd is used for coordinates
DirNd is used for vector directions, unit length
MatNd is a N x N matrix
TfmNd is used to package together a rotation matrix and translation vector

Created on Mon Jun  2 10:48:54 2025

.. codeauthor: Michael J. Hayford
"""
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
