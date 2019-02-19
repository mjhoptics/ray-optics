#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
"""

.. Created on Sat Jun 30 10:15:37 2018

.. codeauthor: Michael J. Hayford
"""

import warnings
import matplotlib.cbook

from matplotlib.figure import Figure
from matplotlib.patches import Polygon

import numpy as np

import rayoptics.gui.layout as layout

Fit_All, User_Scale = range(2)


warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)


def rgb2mpl(rgb):
    if len(rgb) == 3:
        return [rgb[0]/255., rgb[1]/255., rgb[2]/255.]
    elif len(rgb) == 4:
        return [rgb[0]/255., rgb[1]/255., rgb[2]/255., rgb[3]/255.]


backgrnd_color = rgb2mpl([237, 243, 254])  # light blue


class LensLayoutFigure(Figure):

    def __init__(self, opt_model,
                 do_draw_frame=False,
                 scale_type=Fit_All,
                 user_scale_value=1.0,
                 oversize_factor=0.05,
                 offset_factor=0.05,
                 do_draw_rays=True,
                 **kwargs):
        self.opt_model = opt_model
        self.user_scale_value = user_scale_value
        self.scale_type = scale_type
        self.linewidth = 0.5
        self.do_draw_frame = do_draw_frame
        self.do_draw_rays = do_draw_rays
        self.oversize_factor = oversize_factor
        self.offset_factor = offset_factor

        Figure.__init__(self, **kwargs)

        self.set_facecolor(backgrnd_color)

        ele_model = self.opt_model.ele_model
        if len(ele_model.elements) == 0:
            ele_model.elements_from_sequence(self.opt_model.seq_model)

        self.update_data()

    def refresh(self):
        self.update_data()
        self.plot()

    def update_data(self):
        self.patches = []

        self.ele_shapes = self.create_element_model(self.opt_model.ele_model)
        self.ele_bbox = self.update_patches(self.ele_shapes)

        if self.do_draw_rays:
            offset = self.offset_factor
            self.ray_shapes = self.create_ray_model(offset_factor=offset)
            self.ray_bbox = self.update_patches(self.ray_shapes)

            concat_bbox = np.concatenate((self.ele_bbox, self.ray_bbox))
            self.sys_bbox = layout.bbox_from_poly(concat_bbox)
        else:
            self.sys_bbox = self.ele_bbox

        return self

    def system_length(self):
        osp = self.opt_model.optical_spec
        img_dist = abs(osp.parax_data[2].img_dist)
        ele_length = self.ele_bbox[1][0] - self.ele_bbox[0][0]
        return ele_length+img_dist

    def update_patches(self, shapes):
        bbox_list = []
        for shape in shapes:
            poly, bbox = shape[0](shape[1])
            self.patches.append(poly)
            if len(bbox_list) == 0:
                bbox_list = bbox
            else:
                bbox_list = np.vstack((bbox_list, bbox))
        bbox = layout.bbox_from_poly(bbox_list)
        return bbox

    def create_element_model(self, ele_model):
        elements = []
        for e in ele_model.elements:
            oe = layout.OpticalElement(e)
            elements.append((self.update_element_shape, oe))
        return elements

    def update_element_shape(self, oe):
        poly, bbox = oe.update_shape()
        p = Polygon(poly, closed=True, fc=rgb2mpl(oe.render_color()),
                    ec='black')
        p.set_linewidth(self.linewidth)
        return p, bbox

    def create_ray_model(self, start_surf=1, offset_factor=0.05):
        start_offset = offset_factor*self.system_length()

        ray_bundles = []
        fov = self.opt_model.optical_spec.field_of_view
        wvl = self.opt_model.seq_model.central_wavelength()
        for fld in fov.fields:
            rb = layout.RayBundle(self.opt_model,
                                  fld, wvl, start_offset)
            ray_bundles.append((self.update_ray_fan_shape, rb))
        return ray_bundles

    def update_ray_fan_shape(self, rb):
        rndr_clr = rgb2mpl([254, 197, 254, 64])  # magenta, 25%

        poly, bbox = rb.update_shape()
        p = Polygon(poly, fc=rndr_clr, ec='black')
        p.set_linewidth(self.linewidth)

        return p, bbox

    def scale_bounds(self, oversize_factor):
        bbox = layout.scale_bounds(self.sys_bbox, oversize_factor)
        self.ax.set_xlim(bbox[0][0], bbox[1][0])
        self.ax.set_ylim(bbox[0][1], bbox[1][1])

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
