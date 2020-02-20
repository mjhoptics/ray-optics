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
from matplotlib import lines
from matplotlib import patches
from matplotlib import widgets

from rayoptics.mpl.styledfigure import StyledFigure

from rayoptics.gui import util
from rayoptics.util import rgb2mpl


SelectInfo = namedtuple('SelectInfo', ['artist', 'info'])
""" tuple grouping together an artist and info returned from contains(event)

    Attributes:
        artist: the artist
        info: a dictionary of artist specific details of selection
"""


def display_artist_and_event(callback_str, event, artist):
    if artist and hasattr(artist, 'shape'):
        shape, handle = artist.shape
        print('{} {}: shape: {}.{}'.format(callback_str, event.name,
                                           shape.get_label(), handle))
    else:
        print('{} {}: shape: None'.format(callback_str, event.name))


class InteractiveFigure(StyledFigure):
    """ Editable version of optical system layout, aka Live Layout

    Attributes:
        opt_model: parent optical model
        do_draw_frame: if True, draw frame around the figure
        do_draw_axes: if True, draw coordinate axes for the figure
        oversize_factor: what fraction to oversize the view bounding box
        aspect: 'equal' for 1:1 aspect ratio, 'auto' for best ratio
    """

    def __init__(self,
                 do_draw_frame=False,
                 do_draw_axes=False,
                 oversize_factor=0.05,
                 aspect='equal',
                 view_bbox=None,
                 do_scale_bounds=False,
                 **kwargs):
        self.linewidth = 0.5
        self.do_draw_frame = do_draw_frame
        self.do_draw_axes = do_draw_axes
        self.oversize_factor = oversize_factor
        self.aspect = aspect
        self.hilited = None
        self.selected = None
        self.do_scale_bounds = do_scale_bounds

        self.do_action = self.do_shape_action
        self.event_dict = {}

        self.on_finished = None

        super().__init__(**kwargs)

        self.update_data()
        self.view_bbox = view_bbox if view_bbox else self.fit_axis_limits()

    def connect_events(self, action_dict=None):
        'connect to all the events we need'
        if action_dict is None:
            action_dict = {'motion_notify_event': self.on_motion,
                           'button_press_event': self.on_press,
                           'button_release_event': self.on_release,
                           # 'button_press_event': self.display_event,
                           # 'button_release_event': self.display_event,
                           # 'pick_event': self.on_select,
                           }
        self.callback_ids = []
        for event, action in action_dict.items():
            self.event_dict[event] = action
            cid = self.canvas.mpl_connect(event, action)
            self.callback_ids.append(cid)

    def disconnect_events(self):
        'disconnect all the stored connection ids'
        for clbk in self.callback_ids:
            self.canvas.mpl_disconnect(clbk)
        self.callback_ids = None
        event_dict, self.event_dict = self.event_dict, {}
        return event_dict

    @property
    def is_unit_aspect_ratio(self):
        return self.aspect == 'equal'

    @is_unit_aspect_ratio.setter
    def is_unit_aspect_ratio(self, value):
        self.aspect = 'equal' if value else 'auto'

    def refresh(self):
        self.update_data()
        self.plot()

    def update_data(self):
        pass

    def action_complete(self):
        if self.on_finished:
            self.on_finished()
            self.on_finished = None

    def register_action(self, *args, **kwargs):
        action_obj = args[0]
        fig = kwargs.pop('figure')

        def do_command_action(event, target, event_key):
            nonlocal action_obj, fig
            try:
                action_obj.actions[event_key](fig, event)
            except KeyError:
                pass
        self.do_action = do_command_action

    def register_pan(self, on_finished):
        action_obj = PanAction()
        self.register_action(action_obj, figure=self)
        self.on_finished = on_finished

    def register_zoom_box(self, on_finished):
        self.zoom_box_action = ZoomBoxAction(self)

        def do_command_action(event, target, event_key):
            pass

        self.do_action = do_command_action
        self.on_finished = on_finished

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
        bbox = util.bbox_from_poly(bbox_list)
        return bbox

    def create_patches(self, handles):
        gui_handles = {}
        for key, graphics_handle in handles.items():
            poly_data, poly_type, kwargs = graphics_handle
            poly = np.array(poly_data)
            if poly_type == 'vertex':
                p = self.create_vertex(poly, **kwargs)
            elif poly_type == 'polyline':
                p = self.create_polyline(poly, **kwargs)
            elif poly_type == 'polygon':
                p = self.create_polygon(poly, **kwargs)
            else:
                break
            if len(poly.shape) > 1:
                bbox = util.bbox_from_poly(poly)
            else:
                x = poly[0]
                y = poly[1]
                bbox = np.array([[x, y], [x, y]])
            gui_handles[key] = util.GUIHandle(p, bbox)
        return gui_handles

    def create_polygon(self, poly, **kwargs):
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
        fill_color = kwargs.pop('fill_color',
                                self._rgb['background1'] + '40')  # 25%
        if not isinstance(fill_color, str):
            fill_color = rgb2mpl.rgb2mpl(fill_color)
        p = patches.Polygon(poly, closed=True,
                            facecolor=fill_color,
                            edgecolor=self._rgb['foreground'],
                            **kwargs)
        p.highlight = highlight
        p.unhighlight = unhighlight
        return p

    def create_polyline(self, poly, **kwargs):
        def highlight(p):
            lw = p.get_linewidth()
            ls = p.get_linestyle()
            c = p.get_color()
            p.unhilite = (c, lw, ls)

            p.update(hilite_kwargs)

        def unhighlight(p):
            c, lw, ls = p.unhilite
            p.set_linewidth(lw)
            p.set_linestyle(ls)
            p.set_color(c)
            p.unhilite = None

        hilite_kwargs = {}
        if 'hilite' in kwargs:
            if isinstance(kwargs['hilite'], dict):
                hilite_kwargs = kwargs.pop('hilite')
                if 'linewidth' not in hilite_kwargs:
                    hilite_kwargs['linewidth'] = 5
            else:
                hilite_kwargs['color'] = kwargs.pop('hilite')
                hilite_kwargs['linewidth'] = hilite_kwargs.get('linewidth', 5)
        else:
            hilite_kwargs['color'] = self._rgb['hilite']

        kwargs['linewidth'] = kwargs.get('linewidth', self.linewidth)
        x = poly.T[0]
        y = poly.T[1]
        p = lines.Line2D(x, y, **kwargs)

        p.highlight = highlight
        p.unhighlight = unhighlight
        return p

    def create_vertex(self, vertex, **kwargs):
        def highlight(p):
            mkr = p.get_marker()
            ms = p.get_markersize()
            mfc = p.get_markerfacecolor()
            mec = p.get_markeredgecolor()
            mew = p.get_markeredgewidth()
            p.unhilite = (mkr, ms, mfc, mec, mew)

            p.update(hilite_kwargs)

        def unhighlight(p):
            mkr, ms, mfc, mec, mew = p.unhilite
            p.set_marker(mkr)
            p.set_markersize(ms)
            p.set_markerfacecolor(mfc)
            p.set_markeredgecolor(mec)
            p.set_markeredgewidth(mew)
            p.unhilite = None

        hilite_kwargs = {}
        hilite_kwargs['linewidth'] = hilite_kwargs.get('linewidth', 5)
        if 'hilite' in kwargs:
            if isinstance(kwargs['hilite'], dict):
                hilite_kwargs = kwargs.pop('hilite')
                if 'color' in hilite_kwargs:
                    hc = hilite_kwargs['color']
                    hilite_kwargs['markerfacecolor'] = hc
                    hilite_kwargs['markeredgecolor'] = hc

            else:
                hc = kwargs.pop('hilite')
                hilite_kwargs['markerfacecolor'] = hc
                hilite_kwargs['markeredgecolor'] = hc
                hilite_kwargs['markersize'] = 10
        else:
            hilite_kwargs['markerfacecolor'] = self._rgb['hilite']
            hilite_kwargs['markeredgecolor'] = self._rgb['hilite']
            hilite_kwargs['markersize'] = 10

        kwargs['linewidth'] = kwargs.get('linewidth', self.linewidth)
        x = [vertex[0]]
        y = [vertex[1]]
        p = lines.Line2D(x, y, **kwargs)

        p.highlight = highlight
        p.unhighlight = unhighlight
        return p

    def update_axis_limits(self, bbox):
        self.ax.set_xlim(bbox[0][0], bbox[1][0])
        self.ax.set_ylim(bbox[0][1], bbox[1][1])

    def fit_axis_limits(self):
        """ returns a numpy bounding box that fits the current data """
        pass

    def set_view_bbox(self, bbox):
        self.view_bbox = bbox
        self.update_axis_limits(bbox=self.view_bbox)

    def fit(self):
        self.set_view_bbox(self.fit_axis_limits())
        self.plot()

    def zoom(self, factor):
        bbox = self.view_bbox
        # calculate the bbox half-widths
        hlf_x, hlf_y = (bbox[1][0] - bbox[0][0])/2, (bbox[1][1] - bbox[0][1])/2
        # calculate the center of the bbox
        cen_x, cen_y = (bbox[1][0] + bbox[0][0])/2, (bbox[1][1] + bbox[0][1])/2
        # scale the bbox dimensions by the requested factor
        hlf_x *= factor
        hlf_y *= factor
        # rebuild the scaled bbox
        view_bbox = np.array([[cen_x-hlf_x, cen_y-hlf_y],
                              [cen_x+hlf_x, cen_y+hlf_y]])

        self.set_view_bbox(view_bbox)
        self.plot()

    def zoom_in(self):
        self.zoom(factor=0.8)

    def zoom_out(self):
        self.zoom(factor=1.2)

    def draw_frame(self, do_draw_frame):
        if do_draw_frame:
            self.ax.set_axis_on()
        else:
            self.ax.set_axis_off()
            self.tight_layout(pad=0.)

    def draw_axes(self, do_draw_axes):
        if do_draw_axes:
            self.ax.grid(True)
            self.ax.axvline(0, c=self._rgb['foreground'], lw=1)
            self.ax.axhline(0, c=self._rgb['foreground'], lw=1)
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
            if isinstance(a, lines.Line2D):
                self.ax.add_line(a)
            elif isinstance(a, patches.Patch):
                self.ax.add_patch(a)
            else:
                self.ax.add_artist(a)

        if self.do_scale_bounds:
            self.view_bbox = util.scale_bounds(self.sys_bbox,
                                               self.oversize_factor)

        self.ax.set_aspect(self.aspect, adjustable='datalim')
        self.update_axis_limits(bbox=self.view_bbox)

        self.draw_frame(self.do_draw_frame)
        self.ax.set_facecolor(self._rgb['background1'])

        self.draw_axes(self.do_draw_axes)

        self.connect_events()
        self.canvas.draw_idle()

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
                        logging.debug("on motion, artist {}: {}.{}, z={}, "
                                      "hits={}".format(len(artists),
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

    def on_select(self, event):
        artist = event.artist
        display_artist_and_event('on_select', event, artist)

    def display_event(self, event):
        artist = self.hilited.artist if self.hilited else None
        display_artist_and_event('display_event', event, artist)

    def on_press(self, event):
        self.save_do_scale_bounds = self.do_scale_bounds
        self.do_scale_bounds = False
        target_artist = self.selected = self.hilited
        # if target_artist is not None:
        #     display_artist_and_event('on_press', event, self.selected.artist)
        self.do_action(event, target_artist, 'press')

    def on_motion(self, event):
        if self.selected is None:
            artists = self.find_artists_at_location(event)
            next_hilited = artists[0] if len(artists) > 0 else None
            # if next_hilited is not None:
            #     display_artist_and_event('on_motion ({})'.format(len(artists)),
            #                              event, next_hilited.artist)

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
            # display_artist_and_event('on_drag', event, self.selected.artist)
            self.do_action(event, self.selected, 'drag')
            shape, handle = self.selected.artist.shape
            logging.debug("on_drag: %s %s %d", shape.get_label(), handle,
                          self.selected.artist.get_zorder())

    def on_release(self, event):
        'on release we reset the press data'
        logging.debug("on_release")
        # display_artist_and_event('on_release', event, self.selected.artist)

        self.do_action(event, self.selected, 'release')
        self.do_scale_bounds = self.save_do_scale_bounds
        self.selected = None
        self.action_complete()


class PanAction():
    ''' wrapper class to handle pan action, handing off to Axes '''

    def __init__(self, **kwargs):

        def on_press(fig, event):
            self._button_pressed = event.button
            fig.ax.start_pan(event.x, event.y, event.button)

        def on_drag(fig, event):
            fig.ax.drag_pan(self._button_pressed, event.key, event.x, event.y)
            fig.canvas.draw_idle()

        def on_release(fig, event):
            fig.ax.end_pan()
            # update figure view_bbox with result of pan action
            x_min, x_max = fig.ax.get_xbound()
            y_min, y_max = fig.ax.get_ybound()
            fig.view_bbox = np.array([[x_min, y_min], [x_max, y_max]])

        self.actions = {}
        self.actions['press'] = on_press
        self.actions['drag'] = on_drag
        self.actions['release'] = on_release


class ZoomBoxAction():
    """ handle zoom box action by using a RectangleSelector widget """

    def __init__(self, fig, **kwargs):
        def on_release(press_event, release_event):
            bbox = np.array([[press_event.xdata, press_event.ydata],
                             [release_event.xdata, release_event.ydata]])
            fig.set_view_bbox(bbox)
            fig.canvas.draw_idle()
            self.rubber_box.disconnect_events()
            fig.connect_events(self.saved_events)
            fig.action_complete()

        self.saved_events = fig.disconnect_events()
        rectprops = dict(edgecolor=fig._rgb['foreground'], fill=False)
        self.rubber_box = widgets.RectangleSelector(
            fig.ax, on_release, drawtype='box', useblit=False,
            button=[1, 3],  # don't use middle button
            minspanx=5, minspany=5, spancoords='pixels', rectprops=rectprops,
            interactive=False)
