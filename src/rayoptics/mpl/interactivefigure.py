#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
"""

.. Created on Thu Oct 10 21:57:25 2019

.. codeauthor: Michael J. Hayford
"""

import logging
from collections import namedtuple

import numpy as np
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon

from rayoptics.gui.util import bbox_from_poly, scale_bounds

from rayoptics.util.rgb2mpl import rgb2mpl, backgrnd_color


SelectInfo = namedtuple('SelectInfo', ['artist', 'info'])
""" tuple grouping together an artist and info returned from contains(event)

    Attributes:
        artist: the artist
        info: a dictionary of artist specific details of selection
"""


class InteractiveFigure(Figure):
    """ Editable version of optical system layout, aka Live Layout

    Attributes:
        opt_model: parent optical model
        refresh_gui: function to be called on refresh_gui event
        do_draw_frame: if True, draw frame around system layout
        oversize_factor: what fraction to oversize the system bounding box
        do_draw_rays: if True, draw edge rays
        do_paraxial_layout: if True, draw editable paraxial axial and chief ray
    """
    def __init__(self,
                 do_draw_frame=False,
                 do_draw_axes=False,
                 oversize_factor=0.05,
                 aspect='equal',
                 **kwargs):
        self.linewidth = 0.5
        self.do_draw_frame = do_draw_frame
        self.do_draw_axes = do_draw_axes
        self.oversize_factor = oversize_factor
        self.aspect = aspect
        self.hilited = None
        self.selected = None
        self.do_scale_bounds = True
        self.do_action = self.do_shape_action

        super().__init__(**kwargs)

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
        pass

    def action_complete(self):
        pass

    def update_patches(self, shapes):
        """ loop over the input shapes, fetching their current geometry
        and attaching it to the corresponding ``Artist``
        """
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
        bbox = bbox_from_poly(bbox_list)
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
        fill_color = rgb2mpl(kwargs.pop('fill_color', rgb_color))
        p = Polygon(poly, closed=True, fc=fill_color,
                    ec='black', **kwargs)
        p.highlight = highlight
        p.unhighlight = unhighlight
        return p

    def create_polyline(self, poly, **kwargs):
        def highlight(p):
            lw = p.get_linewidth()
            c = p.get_color()
            p.unhilite = (c, lw)
            p.set_linewidth(2)
            p.set_color(hilite_color)

        def unhighlight(p):
            c, lw = p.unhilite
            p.set_linewidth(lw)
            p.set_color(c)
            p.unhilite = None

        x = poly.T[0]
        y = poly.T[1]
        hilite_color = kwargs.pop('hilite', 'red')
        if 'linewidth' not in kwargs:
            kwargs['linewidth'] = self.linewidth
        p = Line2D(x, y, **kwargs)
        p.highlight = highlight
        p.unhighlight = unhighlight
        return p

    def create_vertex(self, vertex, **kwargs):
        def highlight(p):
            lw = p.get_linewidth()
            c = p.get_color()
            p.unhilite = (c, lw)
            p.set_linewidth(2)
            p.set_color(hilite_color)

        def unhighlight(p):
            c, lw = p.unhilite
            p.set_linewidth(lw)
            p.set_color(c)
            p.unhilite = None

        x = [vertex[0]]
        y = [vertex[1]]
        hilite_color = kwargs.pop('hilite', 'red')
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
        else:
            self.ax.set_axis_off()
            self.tight_layout(pad=0.)

    def draw_axes(self, do_draw_axes):
        if do_draw_axes:
            self.ax.grid(True)
            self.ax.axvline(0, c='black', lw=1)
            self.ax.axhline(0, c='black', lw=1)
            if hasattr(self, 'header'):
                self.ax.set_title(self.header, pad=10.0, fontsize=18)
            if hasattr(self, 'x_label'):
                self.ax.set_xlabel(self.x_label)
            if hasattr(self, 'y_label'):
                self.ax.set_ylabel(self.y_label)
        else:
            self.ax.grid(False)

    def plot(self):
        try:
            self.ax.cla()
        except AttributeError:
            self.ax = self.add_subplot(1, 1, 1, aspect=self.aspect)

        for a in self.artists:
            a.set_picker(5)
            if isinstance(a, Line2D):
                self.ax.add_line(a)
            elif isinstance(a, Patch):
                self.ax.add_patch(a)
            else:
                self.ax.add_artist(a)

        if self.do_scale_bounds:
            self.view_bbox = scale_bounds(self.sys_bbox,
                                          self.oversize_factor)
        self.update_axis_limits()

        self.draw_frame(self.do_draw_frame)
        self.ax.set_facecolor(backgrnd_color)

        self.draw_axes(self.do_draw_axes)

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

            cur_art = self.hilited.artist if self.hilited is not None else None
            nxt_art = next_hilited.artist if next_hilited is not None else None
            if nxt_art is not cur_art:
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
                    logging.debug("hilite_change: %s %s %d",
                                  shape.get_label(), handle,
                                  self.hilited.artist.get_zorder())
        else:
            self.do_action(event, self.selected, 'drag')
            shape, handle = self.selected.artist.shape
            logging.debug("on_drag: %s %s %d", shape.get_label(), handle,
                          self.selected.artist.get_zorder())

    def on_release(self, event):
        'on release we reset the press data'
        logging.debug("on_release")

        self.do_action(event, self.selected, 'release')
        self.do_scale_bounds = True
        self.selected = None
        self.action_complete()
