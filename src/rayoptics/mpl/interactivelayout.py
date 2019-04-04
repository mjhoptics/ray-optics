#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
""" Interactive paraxial layout figure

.. Created on Thu Mar 14 10:20:33 2019

.. codeauthor: Michael J. Hayford
"""

import logging
import warnings
import matplotlib.cbook

from matplotlib.figure import Figure
from matplotlib.lines import Line2D
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
        self.hilited_artist = None
        self.selected_artist = None

        Figure.__init__(self, **kwargs)

        self.set_facecolor(backgrnd_color)

        self.update_data()

    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.canvas.mpl_connect('button_press_event',
                                                self.on_press)
        self.cidrelease = self.canvas.mpl_connect('button_release_event',
                                                  self.on_release)
        self.cidmotion = self.canvas.mpl_connect('motion_notify_event',
                                                 self.on_motion)

    def disconnect(self):
        'disconnect all the stored connection ids'
        self.canvas.mpl_disconnect(self.cidpress)
        self.canvas.mpl_disconnect(self.cidrelease)
        self.canvas.mpl_disconnect(self.cidmotion)

    def refresh(self):
        self.update_data()
        self.plot()

    def update_data(self):
        self.artists = []

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
            handles = shape.update_shape(self)
            for key, value in handles.items():
                poly, bbox = value
                poly.shape = (shape, key)
                self.artists.append(poly)
                if len(bbox_list) == 0:
                    bbox_list = bbox
                else:
                    bbox_list = np.vstack((bbox_list, bbox))
        bbox = layout.bbox_from_poly(bbox_list)
        return bbox

    def create_polygon(self, poly, rgb_color, **kwargs):
        def highlight(p):
            fc = p.get_facecolor()
            ec = p.get_edgecolor()
            lw = p.get_linewidth()
            p.unhilite = (fc, ec, lw)
            alpha = fc[3]+0.5
            p.set_facecolor((fc[0], fc[1], fc[2], alpha))

        def unhighlight(p):
            fc, ec, lw = p.unhilite
            p.set_facecolor(fc)
            p.set_edgecolor(ec)
            p.set_linewidth(lw)
            p.unhilite = None
        p = Polygon(poly, closed=True, fc=rgb2mpl(rgb_color),
                    ec='black', **kwargs)
        p.set_linewidth(self.linewidth)
        p.highlight = highlight
        p.unhighlight = unhighlight
        return p

    def create_polyline(self, poly, **kwargs):
        def highlight(p):
            lw = p.get_linewidth()
            c = p.get_color()
            p.unhilite = (c, lw)
            p.set_linewidth(2)
            p.set_color('red')

        def unhighlight(p):
            c, lw = p.unhilite
            p.set_linewidth(lw)
            p.set_color(c)
            p.unhilite = None
        x = poly.T[0]
        y = poly.T[1]
        p = Line2D(x, y, linewidth=self.linewidth)
        p.highlight = highlight
        p.unhighlight = unhighlight
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

        for a in self.artists:
            a.set_picker(5)
            self.ax.add_artist(a)

        self.scale_bounds(self.oversize_factor)
        self.draw_frame(self.do_draw_frame)
        self.ax.set_facecolor(backgrnd_color)

        self.connect()
        self.canvas.draw()

        return self

    def on_press(self, event):
        artists = self.find_artists_at_location(event)
        obj = artists[0] if len(artists) > 0 else None
        self.selected_artist = self.hilited_artist
        if obj is None:
            logging.debug("on_press: no object found")
        else:
            shape, handle = obj.shape
            if id(obj) != id(self.hilited_artist):
                logging.debug('press event: different than hilite object')
            else:
                logging.debug('press event: same as hilite object')
                logging.debug("on_press:", shape.get_label(), handle,
                              obj.get_zorder())
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

    def find_artists_at_location(self, event):
        artists = []
        for artist in self.ax.get_children():
            if hasattr(artist, 'shape'):
                inside, _ = artist.contains(event)
                if inside:
                    shape, handle = artist.shape
                    artists.append(artist)
                    logging.debug("on_motion:", len(artists),
                                  shape.get_label(), handle,
                                  artist.get_zorder())
        return sorted(artists, key=lambda a: a.get_zorder(), reverse=True)

    def on_motion(self, event):
        if self.selected_artist is None:
            artists = self.find_artists_at_location(event)
            next_hilited_artist = artists[0] if len(artists) > 0 else None

            if id(next_hilited_artist) != id(self.hilited_artist):
                if self.hilited_artist:
                    self.hilited_artist.unhighlight(self.hilited_artist)
                    self.hilited_artist.figure.canvas.draw()
                if next_hilited_artist:
                    next_hilited_artist.highlight(next_hilited_artist)
                    next_hilited_artist.figure.canvas.draw()
                self.hilited_artist = next_hilited_artist
                if next_hilited_artist is None:
                    logging.debug("hilite_change: no object found")
                else:
                    shape, handle = self.hilited_artist.shape
                    logging.debug("hilite_change:", shape.get_label(), handle,
                                  self.hilited_artist.get_zorder())
        else:
            shape, handle = self.selected_artist.shape
            logging.debug("on_drag:", shape.get_label(), handle,
                          self.selected_artist.get_zorder())

    def on_release(self, event):
        'on release we reset the press data'
        logging.debug("on_release")
        self.selected_artist = None
