#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
""" Interactive paraxial layout figure

.. Created on Thu Mar 14 10:20:33 2019

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


class InteractiveLayout(Figure):

    def __init__(self, opt_model,
                 do_draw_frame=False,
                 scale_type=Fit_All,
                 user_scale_value=1.0,
                 oversize_factor=0.05,
                 offset_factor=0.05,
                 do_draw_rays=True,
                 **kwargs):
        self.layout = layout.LensLayout(opt_model)
        self.user_scale_value = user_scale_value
        self.scale_type = scale_type
        self.linewidth = 0.5
        self.do_draw_frame = do_draw_frame
        self.do_draw_rays = do_draw_rays
        self.oversize_factor = oversize_factor
        self.offset_factor = offset_factor

        Figure.__init__(self, **kwargs)

        self.set_facecolor(backgrnd_color)

        self.update_data()

    def connect(self):
        'connect to all the events we need'
        self.cidpick = self.canvas.mpl_connect('pick_event', self.on_pick)
#        self.cidpress = self.canvas.mpl_connect('button_press_event',
#                                                self.on_press)
        self.cidrelease = self.canvas.mpl_connect('button_release_event',
                                                  self.on_release)
#        self.cidmotion = self.canvas.mpl_connect('motion_notify_event',
#                                                 self.on_motion)

    def disconnect(self):
        'disconnect all the stored connection ids'
        self.canvas.mpl_disconnect(self.cidpick)
        self.canvas.mpl_disconnect(self.cidpress)
        self.canvas.mpl_disconnect(self.cidrelease)
        self.canvas.mpl_disconnect(self.cidmotion)

    def refresh(self):
        self.update_data()
        self.plot()

    def update_data(self):
        self.patches = []

        self.ele_shapes = self.layout.create_element_model(self)
        self.ele_bbox = self.update_patches(self.ele_shapes)

        if self.do_draw_rays:
            system_length = self.layout.system_length(self.ele_bbox)
            start_offset = self.offset_factor*system_length
            self.ray_shapes = self.layout.create_ray_model(self, start_offset)
            self.ray_bbox = self.update_patches(self.ray_shapes)

            concat_bbox = np.concatenate((self.ele_bbox, self.ray_bbox))
            self.sys_bbox = layout.bbox_from_poly(concat_bbox)
        else:
            self.sys_bbox = self.ele_bbox

        return self

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

    def create_polygon(self, poly, rgb_color, **kwargs):
        p = Polygon(poly, closed=True, fc=rgb2mpl(rgb_color),
                    ec='black', **kwargs)
        p.set_linewidth(self.linewidth)
        return p

    def create_polyline(self, poly, **kwargs):
        p = Polygon(poly, closed=False, ec='black', **kwargs)
        p.set_linewidth(self.linewidth)
        return p

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
            p.set_picker(True)
            self.ax.add_patch(p)

        self.scale_bounds(self.oversize_factor)
        self.draw_frame(self.do_draw_frame)
        self.ax.set_facecolor(backgrnd_color)

        self.connect()
        self.canvas.draw()

        return self

    def on_press(self, event):
#        hit, props = self.line.contains(event)
        print("on_press")
#        if hit:
#            if self.eline.press is None:
#                self.eline.on_press(event)
#            if self.eline.press:
#                v, e = self.eline.press
#                if v:
#                    self.actions['press'](self, v)
#                    hit_vertex, x_hit, y_hit, x_data, y_data, mkr = v
#                    self.vertex = hit_vertex + self.data_slice.start
#                    print("vertex selected", hit_vertex)
#                elif e:
#                    hit_edge, x_hit, y_hit, x_data, y_data, mkr = e
#                    self.actions['press'](self, e)
#                    print("edge selected", e[0])

#            print('on_press', event.button, event.x, event.y,
#                  event.xdata, event.ydata, event.key,
#                  len(hit_list), hit_list)

    def on_motion(self, event):
        pass
#        if self.vertex:
#            self.actions['drag'](self, event)

    def on_release(self, event):
        'on release we reset the press data'
        print("on_release")
#        if self.vertex:
#            self.actions['release'](self, event)

    def on_pick(self, event):
        obj = event.artist
        me = event.mouseevent
        id = obj.get_gid()
        print("on_pick", type(obj).__name__, id, me.name, me.x, me.y,
              me.button, me.key, me.xdata, me.ydata, me.dblclick)
