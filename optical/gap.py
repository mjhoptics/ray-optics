#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 17:06:17 2017

@author: Mike
"""


import medium as m


class Gap:
    def __init__(self, t=0.0, med=m.Air()):
        self.thi = t
        self.medium = med

    def __repr__(self):
        return "Gap(t=%r, medium=%r)" % (self.thi, self.medium)
