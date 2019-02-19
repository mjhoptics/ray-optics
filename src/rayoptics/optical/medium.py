#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2017 Michael J. Hayford
""" Module for simple optical media definitions

.. Created on Fri Sep 15 17:06:17 2017

.. codeauthor: Michael J. Hayford
"""


def glass_encode(n, v):
    return str(1000*round((n - 1), 3) + round(v/100, 3))


def glass_decode(gc):
    return round(1.0 + (int(gc)/1000), 3), round(100.0*(gc - int(gc)), 3)


class Medium:
    """ Constant refractive index medium. """
    def __init__(self, nd, lbl):
        self.label = lbl
        self.n = nd

    def __repr__(self):
        return 'Medium ' + self.label + ': ' + str(self.n)

    def name(self):
        return self.label

    def rindex(self, wv_nm):
        return self.n


class Air(Medium):
    """ Optical definition for air (low fidelity definition) """
    def __init__(self):
        self.label = 'air'
        self.n = 1.0

    def __repr__(self):
        return 'Air'


class Glass(Medium):
    """ Optical medium defined by a glass code, i.e. index - V number pair """
    def __init__(self, nd=1.5168, vd=64.17, mat='N-BK7'):
        self.label = mat
        if mat == 'N-BK7':
            self.n = 1.5168
            self.v = 64.17
        else:
            self.n = nd
            self.v = vd

    def __repr__(self):
        return 'Glass ' + self.label + ': ' + glass_encode(self.n, self.v)

    def glass_code(self):
        return str(1000*round((self.n - 1), 3) + round(self.v/100, 3))

    def name(self):
        if self.label == '':
            return glass_encode(self.n, self.v)
        else:
            return self.label

    def rindex(self, wv_nm):
        return self.n
