#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 16:27:01 2018

@author: Mike
"""


class Element():
    def __init__(self, tfrm, s1, g, s2, sd):
        self.tfrm = tfrm
        self.s1 = s1
        self.s2 = s2
        self.thi = g.thi
        self.medium = g.medium
        self.sd = sd
        self.flat1 = None
        self.flat2 = None

    def shape(self):
        poly = self.s1.full_profile(self.sd, self.flat1)
        poly2 = self.s2.full_profile(self.sd, self.flat2, -1)
        for p in poly2:
            p[0] += self.thi
        poly += poly2
        poly.append(poly[0])
        return poly


class ElementModel:

    def __init__(self):
        self.elements = []

    def reset(self):
        self.__init__()

    def elements_from_sequence(self, seq_model):
        tfrms = seq_model.compute_global_coords(1)
        for i, g in enumerate(seq_model.gaps):
            if g.medium.name().lower() == 'air':
                # close off element
                s2 = seq_model.surfs[i+1]
                if s2.refract_mode is 'REFL':
                    tfrm = tfrms[i+1]
                    sd = s2.surface_od()
                    self.elements.append(Element(tfrm, None, g, s2, sd))
            else:
                tfrm = tfrms[i]
                s1 = seq_model.surfs[i]
                s2 = seq_model.surfs[i+1]
                sd = max(s1.surface_od(), s2.surface_od())
                self.elements.append(Element(tfrm, s1, g, s2, sd))

    def get_num_elements(self):
        return len(self.elements)

    def list_elements(self):
        for i, ele in enumerate(self.elements):
            if ele.s1 is not None:
                print(ele.s1.profile,
                      ele.s2.profile,
                      ele.thi, ele.medium.name())
            else:
                print('REFL',
                      ele.s2.profile)
