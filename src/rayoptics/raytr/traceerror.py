#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Support for ray trace exception handling

.. Created on Wed Oct 24 15:22:40 2018

.. codeauthor: Michael J. Hayford
"""


class TraceError(Exception):
    """ Exception raised when ray tracing a model """


class TraceMissedSurfaceError(TraceError):
    """ Exception raised when ray misses an interface """
    def __init__(self, ifc=None, prev_seg=None):
        self.ifc = ifc
        self.prev_seg = prev_seg


class TraceTIRError(TraceError):
    """ Exception raised when ray TIRs at an interface """
    def __init__(self, inc_dir, normal, prev_indx, follow_indx):
        self.ifc = None
        self.int_pt = None
        self.inc_dir = inc_dir
        self.normal = normal
        self.prev_indx = prev_indx
        self.follow_indx = follow_indx


class TraceEvanescentRayError(TraceError):
    """ Exception raised when ray diffracts evanescently at an interface """
    def __init__(self, ifc, int_pt, inc_dir, normal, prev_indx, follow_indx):
        self.ifc = ifc
        self.int_pt = int_pt
        self.inc_dir = inc_dir
        self.normal = normal
        self.prev_indx = prev_indx
        self.follow_indx = follow_indx


class TraceRayBlockedError(TraceError):
    """ Exception raised when ray is blocked by an aperture on an interface """
    def __init__(self, ifc, int_pt):
        self.ifc = ifc
        self.int_pt = int_pt
