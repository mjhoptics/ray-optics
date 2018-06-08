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

    def __init__(self, tfrm, s1, s1_indx, g, s2, s2_indx, sd):
        self.tfrm = tfrm
        self.s1 = s1
        self.s1_indx = s1_indx
        self.s2 = s2
        self.s2_indx = s2_indx
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

    def __json_encode__(self):
        attrs = dict(self.__dict__)
        del attrs['tfrm']
        del attrs['s1']
        del attrs['s2']
        return attrs

    def sync_to_restore(self, surfs, tfrms):
        self.tfrm = tfrms[self.s1_indx]
        self.s1 = surfs[self.s1_indx]
        self.s2 = surfs[self.s2_indx]

    def reference_interface(self):
        return self.s1

    def sync_to_update(self, surfs):
        self.s1_indx = surfs.index(self.s1)
        self.s2_indx = surfs.index(self.s2)

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
    def __init__(self, tfrm, s, s_indx, sd, thi=None):
        self.render_color = (192, 192, 192)
        self.tfrm = tfrm
        self.s = s
        self.s_indx = s_indx
        self.sd = sd
        self.flat = None
        if thi is None:
            self.thi = 0.05*self.sd
        else:
            self.thi = thi

    def __json_encode__(self):
        attrs = dict(self.__dict__)
        del attrs['tfrm']
        del attrs['s']
        return attrs

    def sync_to_restore(self, surfs, tfrms):
        self.tfrm = tfrms[self.s_indx]
        self.s = surfs[self.s_indx]

    def reference_interface(self):
        return self.s

    def sync_to_update(self, surfs):
        self.s_indx = surfs.index(self.s)

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
    def __init__(self, tfrm, intrfc, intrfc_indx):
        self.render_color = (192, 192, 192)
        self.tfrm = tfrm
        self.intrfc = intrfc
        self.intrfc_indx = intrfc_indx
        self.sd = intrfc.od

    def __json_encode__(self):
        attrs = dict(self.__dict__)
        del attrs['tfrm']
        del attrs['intrfc']
        return attrs

    def sync_to_restore(self, surfs, tfrms):
        self.tfrm = tfrms[self.intrfc_indx]
        self.intrfc = surfs[self.intrfc_indx]

    def reference_interface(self):
        return self.intrfc

    def sync_to_update(self, surfs):
        self.intrfc_indx = surfs.index(self.intrfc)

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

    def __json_encode__(self):
        attrs = dict(self.__dict__)
        del attrs['parent']
        return attrs

    def elements_from_sequence(self, seq_model):
        tfrms = seq_model.compute_global_coords(1)
        for i, g in enumerate(seq_model.gaps):
            if isinstance(seq_model.surfs[i], optical.thinlens.ThinLens):
                te = ThinElement(tfrms[i], seq_model.surfs[i], i)
                self.elements.append(te)
                return

            if g.medium.name().lower() == 'air':
                # close off element
                s2 = seq_model.surfs[i+1]
                if s2.refract_mode is 'REFL':
                    tfrm = tfrms[i+1]
                    sd = s2.surface_od()
                    self.elements.append(Mirror(tfrm, s2, i+1, sd))
            else:
                tfrm = tfrms[i]
                s1 = seq_model.surfs[i]
                s2 = seq_model.surfs[i+1]
                sd = max(s1.surface_od(), s2.surface_od())
                self.elements.append(Element(tfrm, s1, i, g, s2, i+1, sd))

    def sync_to_restore(self, opt_model):
        self.parent = opt_model
        seq_model = opt_model.seq_model
        surfs = seq_model.surfs
        tfrms = seq_model.compute_global_coords(1)
        for e in self.elements:
            e.sync_to_restore(surfs, tfrms)

    def update_model(self):
        seq_model = self.parent.seq_model
        tfrms = seq_model.compute_global_coords(1)
        for e in self.elements:
            e.update_size()
            e.sync_to_update(seq_model.surfs)
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
