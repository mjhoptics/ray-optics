#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Module for element modeling

Created on Sun Jan 28 16:27:01 2018

@author: Michael J. Hayford
"""

import itertools
import util.rgbtable as rgbt
from optical.model_constants import Surf, Gap
import optical.thinlens


class Element():
    clut = rgbt.RGBTable(filename='util/red_blue64.csv',
                         data_range=[10.0, 100.])

    def __init__(self, tfrm, s1, g, s2, sd):
        self.tfrm = tfrm
        self.s1 = s1
        self.s2 = s2
        self.thi = g.thi
        self.medium = g.medium
        self.sd = sd
        self.flat1 = None
        self.flat2 = None
        try:
            gc = float(self.medium.glass_code())
        except AttributeError:
            self.render_color = (255, 255, 255)  # white
        else:
            # set element color based on V-number
            vnbr = round(100.0*(gc - int(gc)), 3)
            self.render_color = Element.clut.get_color(vnbr)

    def reference_interface(self):
        return self.s1

    def update_size(self):
        self.sd = max(self.s1.surface_od(), self.s2.surface_od())
        return self.sd

    def shape(self):
        poly = self.s1.full_profile(self.sd, self.flat1)
        poly2 = self.s2.full_profile(self.sd, self.flat2, -1)
        for p in poly2:
            p[0] += self.thi
        poly += poly2
        poly.append(poly[0])
        return poly


class Mirror():
    def __init__(self, tfrm, s, sd, thi=None):
        self.render_color = (192, 192, 192)
        self.tfrm = tfrm
        self.s = s
        self.sd = sd
        self.flat = None
        if thi is None:
            self.thi = 0.05*self.sd
        else:
            self.thi = thi

    def reference_interface(self):
        return self.s

    def update_size(self):
        self.sd = self.s.surface_od()
        return self.sd

    def shape(self):
        poly = self.s.full_profile(self.sd, self.flat)
        poly2 = self.s.full_profile(self.sd, self.flat, -1)
        for p in poly2:
            p[0] += self.thi
        poly += poly2
        poly.append(poly[0])
        return poly


class ThinElement():
    def __init__(self, tfrm, intrfc):
        self.render_color = (192, 192, 192)
        self.tfrm = tfrm
        self.intrfc = intrfc
        self.sd = intrfc.od

    def reference_interface(self):
        return self.intrfc

    def update_size(self):
        self.sd = self.intrfc.surface_od()
        return self.sd

    def shape(self):
        poly = self.intrfc.full_profile(self.sd)
        return poly


class ElementModel:

    def __init__(self, opt_model):
        self.parent = opt_model
        self.elements = []

    def reset(self):
        self.__init__()

    def elements_from_sequence(self, seq_model):
        tfrms = seq_model.compute_global_coords(1)
        for i, g in enumerate(seq_model.gaps):
            if isinstance(seq_model.surfs[i], optical.thinlens.ThinLens):
                self.elements.append(ThinElement(tfrms[i], seq_model.surfs[i]))
                return

            if g.medium.name().lower() == 'air':
                # close off element
                s2 = seq_model.surfs[i+1]
                if s2.refract_mode is 'REFL':
                    tfrm = tfrms[i+1]
                    sd = s2.surface_od()
                    self.elements.append(Mirror(tfrm, s2, sd))
            else:
                tfrm = tfrms[i]
                s1 = seq_model.surfs[i]
                s2 = seq_model.surfs[i+1]
                sd = max(s1.surface_od(), s2.surface_od())
                self.elements.append(Element(tfrm, s1, g, s2, sd))

    def update_model(self):
        seq_model = self.parent.seq_model
        tfrms = seq_model.compute_global_coords(1)
        for e in self.elements:
            e.update_size()
            intrfc = e.reference_interface()
            try:
                i = seq_model.surfs.index(intrfc)
            except ValueError:
                print("Interface {} not found".format(intrfc.lbl))
            else:
                e.tfrm = tfrms[i]

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
