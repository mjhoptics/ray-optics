#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Static optical system layout

.. Created on Sat Jun 30 10:15:37 2018

.. codeauthor: Michael J. Hayford
"""

from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon

import numpy as np

from rayoptics.gui import layout
from rayoptics.gui import util
from rayoptics.util.rgb2mpl import rgb2mpl, backgrnd_color


class LensLayoutFigure(Figure):
    """  Static optical system layout

    Attributes:
        opt_model: parent optical model
        do_draw_frame: if True, draw frame around system layout
        oversize_factor: what fraction to oversize the system bounding box
        offset_factor: how much to draw rays before first surface
        do_draw_rays: if True, draw edge rays
    """
    def __init__(self, opt_model,
                 do_draw_frame=False,
                 oversize_factor=0.05,
                 offset_factor=0.05,
                 do_draw_rays=True,
                 **kwargs):
        self.opt_model = opt_model
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
        self.artists = []

        self.ele_shapes = self.create_element_model(self.opt_model.ele_model)
        self.ele_bbox = self.update_patches(self.ele_shapes)

        if self.do_draw_rays:
            offset = self.offset_factor
            self.ray_shapes = self.create_ray_model(offset_factor=offset)
            self.ray_bbox = self.update_patches(self.ray_shapes)

            concat_bbox = np.concatenate((self.ele_bbox, self.ray_bbox))
            self.sys_bbox = util.bbox_from_poly(concat_bbox)
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
            self.artists.append(poly)
            if len(bbox_list) == 0:
                bbox_list = bbox
            else:
                bbox_list = np.vstack((bbox_list, bbox))
        bbox = util.bbox_from_poly(bbox_list)
        return bbox

    def create_element_model(self, ele_model):
        elements = []
        for e in ele_model.elements:
            oe = layout.create_optical_element(self.opt_model, e)
            handles = oe.update_shape(self)
            if 'shape' in handles:
                elements.append((self.update_element_shape, oe))
        return elements

    def update_element_shape(self, oe):
        handles = oe.update_shape(self)
        return handles['shape']

    def create_ray_model(self, start_surf=1, offset_factor=0.05):
        start_offset = offset_factor*self.system_length()

        ray_bundles = []
        fov = self.opt_model.optical_spec.field_of_view
        wvl = self.opt_model.seq_model.central_wavelength()
        for i, fld in enumerate(fov.fields):
            fld_label = fov.index_labels[i]
            rb = layout.RayBundle(self.opt_model, fld, fld_label,
                                  wvl, start_offset)
            ray_bundles.append((self.update_ray_fan_shape, rb))
        return ray_bundles

    def create_polygon(self, poly, **kwargs):
        rgb_color = kwargs.pop('fill_color')
        p = Polygon(poly, closed=True, fc=rgb2mpl(rgb_color),
                    ec='black', **kwargs)
        p.set_linewidth(self.linewidth)
        return p

    def create_polyline(self, poly, **kwargs):
        x = poly.T[0]
        y = poly.T[1]
        p = Line2D(x, y, linewidth=2)
        return p

    def update_ray_fan_shape(self, rb):
        handles = rb.update_shape(self)
        return handles['shape']

    def scale_bounds(self, oversize_factor):
        bbox = util.scale_bounds(self.sys_bbox, oversize_factor)
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

        for p in self.artists:
            if isinstance(p, Line2D):
                self.ax.add_line(p)
            elif isinstance(p, Patch):
                self.ax.add_patch(p)
            else:
                self.ax.add_artist(p)

        self.scale_bounds(self.oversize_factor)
        self.draw_frame(self.do_draw_frame)
        self.ax.set_facecolor(backgrnd_color)

        self.canvas.draw()

        return self
