#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
""" Interactive paraxial layout figure

.. Created on Thu Mar 14 10:20:33 2019

.. codeauthor: Michael J. Hayford
"""

import logging
import warnings
from collections import namedtuple

import matplotlib.cbook

from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon

import numpy as np

import rayoptics.gui.layout as layout

Fit_All, User_Scale = range(2)


warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)

SelectInfo = namedtuple('SelectInfo', ['artist', 'info'])
""" tuple grouping together an artist and info returned from contains(event)

    Attributes:
        artist: the artist
        info: a dictionary of artist specific details of selection
"""


def rgb2mpl(rgb):
    if len(rgb) == 3:
        return [rgb[0]/255., rgb[1]/255., rgb[2]/255.]
    elif len(rgb) == 4:
        return [rgb[0]/255., rgb[1]/255., rgb[2]/255., rgb[3]/255.]


backgrnd_color = rgb2mpl([237, 243, 254])  # light blue


class InteractiveLayout(Figure):

    def __init__(self, opt_model, refresh_gui,
                 do_draw_frame=False,
                 scale_type=Fit_All,
                 user_scale_value=1.0,
                 oversize_factor=0.05,
                 offset_factor=0.05,
                 do_draw_rays=False,
                 do_paraxial_layout=True,
                 **kwargs):
        self.refresh_gui = refresh_gui
        self.layout = layout.LensLayout(opt_model)
        self.user_scale_value = user_scale_value
        self.scale_type = scale_type
        self.linewidth = 0.5
        self.do_draw_frame = do_draw_frame
        self.do_draw_rays = do_draw_rays
        self.do_paraxial_layout = do_paraxial_layout
        self.oversize_factor = oversize_factor
        self.offset_factor = offset_factor
        self.hilited = None
        self.selected = None
        self.do_scale_bounds = True
        self.do_action = self.do_shape_action

        Figure.__init__(self, **kwargs)

        self.set_facecolor(backgrnd_color)

        self.update_data()

    def connect(self, action_dict):
        'connect to all the events we need'
        self.callback_ids = []
        for event, action in action_dict.items():
            self.callback_ids.append(self.canvas.mpl_connect(event, action))

    def disconnect(self):
        'disconnect all the stored connection ids'
        for clbk in self.callback_ids:
            self.canvas.mpl_disconnect(clbk)
        self.callback_ids = None

    def refresh(self):
        self.update_data()
        self.plot()

    def update_data(self):
        self.artists = []
        concat_bbox = []

        self.ele_shapes = self.layout.create_element_model(self)
        self.ele_bbox = self.update_patches(self.ele_shapes)
        concat_bbox.append(self.ele_bbox)

        if self.do_draw_rays:
            system_length = self.layout.system_length(self.ele_bbox)
            start_offset = self.offset_factor*system_length
            self.ray_shapes = self.layout.create_ray_model(self, start_offset)
            self.ray_bbox = self.update_patches(self.ray_shapes)
            concat_bbox.append(self.ray_bbox)

        if self.do_paraxial_layout:
            self.parax_shapes = self.layout.create_paraxial_layout(self)
            self.parax_bbox = self.update_patches(self.parax_shapes)
            concat_bbox.append(self.parax_bbox)

        sys_bbox = np.concatenate(concat_bbox)
        self.sys_bbox = layout.bbox_from_poly(sys_bbox)

        return self

    def update_patches(self, shapes):
        bbox_list = []
        for shape in shapes:
            handles = shape.update_shape(self)
            for key, gui_handle in handles.items():
                poly, bbox = gui_handle
                # add shape and handle key as attribute on artist
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
            if alpha > 1.0:
                alpha -= 1.0  # subtract 0.5 instead of adding
            p.set_facecolor((fc[0], fc[1], fc[2], alpha))

        def unhighlight(p):
            fc, ec, lw = p.unhilite
            p.set_facecolor(fc)
            p.set_edgecolor(ec)
            p.set_linewidth(lw)
            p.unhilite = None

        if 'linewidth' not in kwargs:
            kwargs['linewidth'] = self.linewidth
        p = Polygon(poly, closed=True, fc=rgb2mpl(rgb_color),
                    ec='black', **kwargs)
        p.highlight = highlight
        p.unhighlight = unhighlight
        return p

    def create_polyline(self, poly, hilite='red', **kwargs):
        def highlight(p):
            lw = p.get_linewidth()
            c = p.get_color()
            p.unhilite = (c, lw)
            p.set_linewidth(2)
            p.set_color(hilite)

        def unhighlight(p):
            c, lw = p.unhilite
            p.set_linewidth(lw)
            p.set_color(c)
            p.unhilite = None

        x = poly.T[0]
        y = poly.T[1]
        if 'linewidth' not in kwargs:
            kwargs['linewidth'] = self.linewidth
        p = Line2D(x, y, **kwargs)
        p.highlight = highlight
        p.unhighlight = unhighlight
        return p

    def update_axis_limits(self):
        self.ax.set_xlim(self.view_bbox[0][0], self.view_bbox[1][0])
        self.ax.set_ylim(self.view_bbox[0][1], self.view_bbox[1][1])

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
            if isinstance(a, Line2D):
                self.ax.add_line(a)
            elif isinstance(a, Patch):
                self.ax.add_patch(a)
            else:
                self.ax.add_artist(a)

        if self.do_scale_bounds:
            self.view_bbox = layout.scale_bounds(self.sys_bbox,
                                                 self.oversize_factor)
        self.update_axis_limits()

        self.draw_frame(self.do_draw_frame)
        self.ax.set_facecolor(backgrnd_color)

        default_actions = {'button_press_event': self.on_press,
                           'button_release_event': self.on_release,
                           'motion_notify_event': self.on_motion}

        self.connect(default_actions)
        self.canvas.draw()

        return self

    def find_artists_at_location(self, event):
        artists = []
        for artist in self.ax.get_children():
            if hasattr(artist, 'shape'):
                inside, info = artist.contains(event)
                if inside:
                    shape, handle = artist.shape
                    artists.append(SelectInfo(artist, info))
                    if 'ind' in info:
                        logging.debug("on motion, artist {}: {}.{}, z={}, \
                                      hits={}".format(len(artists),
                                      shape.get_label(), handle,
                                      artist.get_zorder(), info['ind']))
                    else:
                        logging.debug("on motion, artist {}: {}.{}, z={}"
                                      .format(len(artists), shape.get_label(),
                                              handle, artist.get_zorder()))

        return sorted(artists, key=lambda a: a.artist.get_zorder(),
                      reverse=True)

    def do_shape_action(self, event, target, event_key):
        if target is not None:
            shape, handle = target.artist.shape
            try:
                action = shape.actions[event_key]
                action(self, handle, event, target.info)
            except KeyError:
                pass

    def on_press(self, event):
        self.do_scale_bounds = False
        target_artist = self.selected = self.hilited
        self.do_action(event, target_artist, 'press')

    def on_motion(self, event):
        if self.selected is None:
            artists = self.find_artists_at_location(event)
            next_hilited = artists[0] if len(artists) > 0 else None

            if next_hilited is not self.hilited:
                if self.hilited:
                    self.hilited.artist.unhighlight(self.hilited.artist)
                    self.hilited.artist.figure.canvas.draw()
                if next_hilited:
                    next_hilited.artist.highlight(next_hilited.artist)
                    next_hilited.artist.figure.canvas.draw()
                self.hilited = next_hilited
                if next_hilited is None:
                    logging.debug("hilite_change: no object found")
                else:
                    shape, handle = self.hilited.artist.shape
                    logging.debug("hilite_change:", shape.get_label(), handle,
                                  self.hilited.artist.get_zorder())
        else:
            self.do_action(event, self.selected, 'drag')
            shape, handle = self.selected.artist.shape
            logging.debug("on_drag:", shape.get_label(), handle,
                          self.selected.artist.get_zorder())

    def on_release(self, event):
        'on release we reset the press data'
        logging.debug("on_release")

        self.do_action(event, self.selected, 'release')
        self.do_scale_bounds = True
        self.selected = None
        self.do_action = self.do_shape_action
