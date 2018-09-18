#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
"""
Created on Sat Jun 30 10:15:37 2018

@author: Michael J. Hayford
"""

import warnings
import matplotlib.cbook

from matplotlib.figure import Figure
from matplotlib.patches import Polygon

import numpy as np

Fit_All, User_Scale = range(2)


warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)


def rgb2mpl(rgb):
    if len(rgb) == 3:
        return [rgb[0]/255., rgb[1]/255., rgb[2]/255.]
    elif len(rgb) == 4:
        return [rgb[0]/255., rgb[1]/255., rgb[2]/255., rgb[3]/255.]


backgrnd_color = rgb2mpl([237, 243, 254])  # light blue


def bbox_from_poly(poly):
    minx, miny = np.min(poly, axis=0)
    maxx, maxy = np.max(poly, axis=0)
    return np.array([[minx, miny], [maxx, maxy]])


class LensLayoutFigure(Figure):

    def __init__(self, opt_model,
                 do_draw_frame=False,
                 scale_type=Fit_All,
                 user_scale_value=1.0,
                 oversize_factor=0.05,
                 **kwargs):
        self.opt_model = opt_model
        self.user_scale_value = user_scale_value
        self.scale_type = scale_type
        self.linewidth = 0.5
        self.do_draw_frame = do_draw_frame
        self.oversize_factor = oversize_factor

        Figure.__init__(self, **kwargs)

        self.set_facecolor(backgrnd_color)

        self.patches = []

        ele_model = self.opt_model.ele_model
        if len(ele_model.elements) == 0:
            ele_model.elements_from_sequence(self.opt_model.seq_model)
        self.update_data()

    def refresh(self):
        self.update_data()
        self.plot()

    def update_data(self):
        self.patches = []
        self.ele_bbox = self.update_element_model(self.opt_model.ele_model)
        self.ray_bbox = self.update_ray_model()
        self.sys_bbox = bbox_from_poly(np.concatenate((self.ele_bbox,
                                                       self.ray_bbox)))
        return self

    def update_element_model(self, ele_model):
        bbox_list = []
        for e in ele_model.elements:
            p, bbox = self.update_element_shape(e)
            self.patches.append(p)
            if len(bbox_list) == 0:
                bbox_list = bbox
            else:
                bbox_list = np.vstack((bbox_list, bbox))
        ele_bbox = bbox_from_poly(bbox_list)
        return ele_bbox

    def update_element_shape(self, e):
        t = np.array([e.tfrm[1][2], -e.tfrm[1][1]])
        poly = np.array(e.shape())
        poly += t
        bbox = bbox_from_poly(poly)
        p = Polygon(poly, closed=True, fc=rgb2mpl(e.render_color), ec='black')
        p.set_linewidth(self.linewidth)
        return p, bbox

    def system_length(self):
        seq_model = self.opt_model.seq_model
        img_dist = abs(seq_model.optical_spec.parax_data[2].img_dist)
        ele_length = self.ele_bbox[1][0] - self.ele_bbox[0][0]
        return ele_length+img_dist

    def update_ray_model(self, start_surf=1):
        start_offset = 0.05*self.system_length()

        bbox_list = []
        fov = self.opt_model.seq_model.optical_spec.field_of_view
        for fi, f in enumerate(fov.fields):
            p, bbox = self.update_ray_fan_shape(fi, start_offset)
            self.patches.append(p)
            if len(bbox_list) == 0:
                bbox_list = bbox
            else:
                bbox_list = np.vstack((bbox_list, bbox))
        ray_bbox = bbox_from_poly(bbox_list)
        return ray_bbox

    def update_ray_fan_shape(self, field_num, start_offset):
        offset = start_offset
        seq_model = self.opt_model.seq_model
        optical_spec = seq_model.optical_spec
        tfrms = seq_model.transforms
        fld, wvl, foc = optical_spec.lookup_fld_wvl_focus(field_num)
        rayset = optical_spec.trace_boundary_rays_at_field(seq_model, fld, wvl)

        # If the object distance (tfrms[0][1][2]) is greater than the
        #  start_offset, then modify rayset start to match start_offset.
        # Remember object transformation for resetting at the end.
        tfrtm0 = tfrms[0]

        if abs(tfrms[0][1][2]) > start_offset:
            r, t = seq_model.setup_shift_of_ray_bundle(offset)
            tfrms[0] = (r, t)
            seq_model.shift_start_of_ray_bundle(rayset, offset, r, t)

        poly1 = []
        for i, r in enumerate(rayset[3][0][0:]):
            rot, trns = tfrms[i]
            p = rot.dot(r[0]) + trns
#            print(i, r[0], rot, trns, p)
#            print("r3", i, p[2], p[1])
            poly1.append([p[2], p[1]])

        poly2 = []
        for i, r in enumerate(rayset[4][0][0:]):
            rot, trns = tfrms[i]
            p = rot.dot(r[0]) + trns
#            print(i, r[0], rot, trns, p)
#            print("r4", i, p[2], p[1])
            poly2.append([p[2], p[1]])

        poly2.reverse()
        poly1.extend(poly2)
        bbox = bbox_from_poly(poly1)

        rndr_clr = rgb2mpl([254, 197, 254, 64])  # magenta, 25%

        p = Polygon(poly1, fc=rndr_clr, ec='black')
        p.set_linewidth(self.linewidth)

        tfrms[0] = tfrtm0

        return p, bbox

    def scale_bounds(self, oversize_factor):
        inc_x = oversize_factor*(self.sys_bbox[1][0] - self.sys_bbox[0][0])
        inc_y = oversize_factor*(self.sys_bbox[1][1] - self.sys_bbox[0][1])
        incr = max(inc_x, inc_y)
        bb = self.sys_bbox
        self.ax.set_xlim(bb[0][0]-incr, bb[1][0]+incr)
        self.ax.set_ylim(bb[0][1]-incr, bb[1][1]+incr)

    def draw_frame(self, do_draw_frame):
        if do_draw_frame:
            self.ax.set_axis_on()
            self.tight_layout(pad=0.)
        else:
            self.ax.set_axis_off()
            self.tight_layout(pad=0.)

    def plot(self):
        try:
            self.ax.cla()
        except AttributeError:
            self.ax = self.add_subplot(1, 1, 1, aspect=1.0)

        for p in self.patches:
            self.ax.add_patch(p)

        self.scale_bounds(self.oversize_factor)
        self.draw_frame(self.do_draw_frame)
        self.ax.set_facecolor(backgrnd_color)

        self.canvas.draw()

        return self
