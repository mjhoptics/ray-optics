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
from matplotlib import collections
from matplotlib import lines
from matplotlib import patches
from matplotlib import widgets
from matplotlib.path import Path
from matplotlib.patches import PathPatch

from rayoptics.mpl.styledfigure import StyledFigure

from rayoptics.gui import util
from rayoptics.util import rgb2mpl
from rayoptics.util import misc_math

logger = logging.getLogger(__name__)

SelectInfo = namedtuple('SelectInfo', ['artist', 'info'])
SelectInfo.artist.__doc__ = "the artist"
SelectInfo.info.__doc__ = "a dictionary of artist specific details of selection"


""" tuple grouping together an artist and info returned from contains(event)

    Attributes:
        artist: the artist
        info: a dictionary of artist specific details of selection
"""

def artist_str(a):
    shape, handle = a.shape
    return f"{a}<{id(a)}>: {shape.listobj_str()}, {handle}, <{id(shape)}>"

def display_artist_and_event(callback_str, event, artist):
    if isinstance(artist, list):
        artist_list = ""
        spacer = ''
        for a in artist:
            if hasattr(a, 'shape'):
                artist_list += spacer + artist_str(a)
                spacer = ', '
            else:
                artist_list += spacer + f"{type(a).__name__}"
                spacer = ', '
    elif artist and hasattr(artist, 'shape'):
        print(f"{callback_str} {event.name}: "
              +artist_str(artist))
    else:
        print(f"{callback_str} {event.name}: shape: None")


class InteractiveFigure(StyledFigure):
    """Base class for domain specific figures with support for input events

    The **InteractiveFigure** class supplies common implementations for:

        - polyline and polygon 2D graphics
        - selection support for mpl graphic objects
        - mouse/touch event handling
        - interface commands for zooming and panning the display area

    Attributes:
        do_draw_frame: if True, draw frame around the figure
        do_draw_axes: if True, draw coordinate axes for the figure
        oversize_factor: what fraction to oversize the view bounding box
        aspect: 'equal' for 1:1 aspect ratio, 'auto' for best ratio
        artist_filter: an (optional) callable applied in
                       find_artists_at_location(), returns True if rejected
    """

    def __init__(self,
                 do_draw_frame=False,
                 do_draw_axes=False,
                 oversize_factor=0.05,
                 aspect='equal',
                 view_bbox=None,
                 do_scale_bounds=False,
                 **kwargs):
        self.do_draw_frame = do_draw_frame
        self.do_draw_axes = do_draw_axes
        self.oversize_factor = oversize_factor
        self.aspect = aspect

        self.hilited_artists: list[SelectInfo] = []
        self.selected_shape = None
        self.do_scale_bounds = do_scale_bounds

        self.artist_filter = None
        self.do_action = self.do_shape_action
        self.event_dict = {}

        self.is_mouse_down = False
        # self.mouse_down_count = 0
        self.on_finished = None

        super().__init__(**kwargs)

        self.update_data()
        self.view_bbox = view_bbox if view_bbox else self.fit_axis_limits()

        self.connect_events()


    def connect_events(self, action_dict=None):
        'connect to all the events we need'
        if action_dict is None:
            action_dict = {'motion_notify_event': self.on_motion,
                           'button_press_event': self.on_press,
                           'button_release_event': self.on_release,
                           'key_press_event': self.on_key_press,
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
        self.callback_ids = []
        event_dict, self.event_dict = self.event_dict, {}
        return event_dict

    @property
    def is_unit_aspect_ratio(self):
        return self.aspect == 'equal'

    @is_unit_aspect_ratio.setter
    def is_unit_aspect_ratio(self, value):
        self.aspect = 'equal' if value else 'auto'

    def refresh(self, **kwargs):
        """Call update_data() followed by plot(), return self.

        Args:
            kwargs: keyword arguments are passed to update_data

        Returns:
            self (class Figure) so scripting envs will auto display results
        """
        self.update_data(**kwargs)
        self.plot()
        return self

    def update_data(self, **kwargs):
        pass

    # --- command support
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

    # --- graphics element creation
    def update_patches(self, shapes):
        """ loop over the input shapes, fetching their current geometry
        and attaching it to the corresponding ``Artist``
        """
        bbox_list = []
        # print('update_patches')
        for shape in shapes:
            handles = shape.update_shape(self)
            for handle, gui_handle in handles.items():
                artist, bbox = gui_handle
                # add shape and handle key as attribute on artist
                if ~np.isnan(bbox).all():
                    # this is the tie between artist and shape
                    artist.shape = (shape, handle)
                    # print(f'shape: {type(shape).__name__} {key}: '
                    #       f'{type(poly).__name__}')
                    self.artists.append(artist)
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
            if isinstance(poly_data, tuple):
                poly = []
                for poly_seg in poly_data:
                    poly.append(np.array(poly_seg))
                poly = tuple(poly)
            else:
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

    def create_polygon(self, poly, linewidth=0.5, **kwargs):
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
            kwargs['linewidth'] = linewidth
        fill_color = kwargs.pop('fill_color',
                                self._rgb['background1'] + '40')  # 25%
        if not isinstance(fill_color, str):
            fill_color = rgb2mpl.rgb2mpl(fill_color)

        if isinstance(poly, tuple):
            p = collections.PolyCollection(poly, closed=True,
                                           facecolor=fill_color, 
                                           edgecolor=self._rgb['foreground'], **kwargs)
        else:
            p = patches.Polygon(poly, closed=True,
                                facecolor=fill_color,
                                edgecolor=self._rgb['foreground'],
                                **kwargs)
        p.highlight = highlight
        p.unhighlight = unhighlight
        # print(f"create_polygon: {p} <{id(p)}>")
        return p

    def create_polyline(self, poly,
                        linewidth=0.5, hilite_linewidth=2,
                        **kwargs):
        def highlight(p):
            # print("highlight: "+artist_str(p))
            lw = p.get_linewidth()
            ls = p.get_linestyle()
            c = p.get_color()
            p.unhilite = (c, lw, ls)

            p.update(hilite_kwargs)

        def unhighlight(p):
            # print("unhighlight: "+artist_str(p))
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
                    hilite_kwargs['linewidth'] = hilite_linewidth
            else:
                hilite_kwargs['color'] = kwargs.pop('hilite')
                hilite_kwargs['linewidth'] = hilite_linewidth
        else:
            hilite_kwargs['color'] = self._rgb['hilite']

        kwargs['linewidth'] = linewidth
        if isinstance(poly, tuple):
            p = collections.LineCollection(poly, **kwargs)
        else:
            x = poly.T[0]
            y = poly.T[1]
            p = lines.Line2D(x, y, **kwargs)

        p.highlight = highlight
        p.unhighlight = unhighlight
        # print(f"create_polyline: {p} <{id(p)}>")
        return p

    def create_vertex(self, vertex,
                      markersize=5, hilite_markersize=7,
                      **kwargs):
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
        # hilite_kwargs['linewidth'] = hilite_kwargs.get('linewidth', 5)
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
                hilite_kwargs['markersize'] = hilite_markersize
        else:
            hilite_kwargs['markerfacecolor'] = self._rgb['hilite']
            hilite_kwargs['markeredgecolor'] = self._rgb['hilite']
            hilite_kwargs['markersize'] = hilite_markersize

        kwargs['markersize'] = markersize
        x = [vertex[0]]
        y = [vertex[1]]
        p = lines.Line2D(x, y, **kwargs)

        p.highlight = highlight
        p.unhighlight = unhighlight
        # print(f"create_vertex: {p} <{id(p)}>")
        return p

    def get_artist_for_handle(self, handle) -> SelectInfo:
        (shape, hdl), info = handle
        for a in self.artists:
            if shape == a.shape:
                return SelectInfo(a, info)
        return None

    def update_axis_limits(self, bbox):
        self.ax.set_xlim(bbox[0][0], bbox[1][0])
        self.ax.set_ylim(bbox[0][1], bbox[1][1])

    def fit_axis_limits(self):
        """ returns a numpy bounding box that fits the current data """
        pass

    def set_view_bbox(self, bbox):
        self.view_bbox = bbox
        self.update_axis_limits(bbox=self.view_bbox)

    # --- 2D view controls
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

    # --- drawing
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
        """Draw the actual figure."""
        try:
            self.ax.cla()
        except AttributeError:
            self.ax = self.add_subplot(1, 1, 1, aspect=self.aspect)

        for a in self.artists:
            if isinstance(a, lines.Line2D):
                a.set_pickradius(5)
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

        self.canvas.draw_idle()

        return self

    # --- interactive actions
    def find_artists_at_location(self, event) -> list[SelectInfo]:
        """Returns a list of shapes in zorder at the event location."""
        artists = []
        for artist in self.ax.get_children():
            if hasattr(artist, 'shape'):
                inside, info = artist.contains(event)
                if inside:
                    shape, handle = artist.shape
                    if self.artist_filter:
                        if self.artist_filter(artist):
                            continue

                    artists.append(SelectInfo(artist, info))
                    if 'ind' in info:
                        logger.debug("on motion, artist {}: {}.{}, z={}, "
                                     "hits={}".format(len(artists),
                                     shape.get_label(), handle,
                                     artist.get_zorder(), info['ind']))
                    else:
                        logger.debug("on motion, artist {}: {}.{}, z={}"
                                     .format(len(artists), shape.get_label(),
                                             handle, artist.get_zorder()))

        return sorted(artists, key=lambda a: a.artist.get_zorder(),
                      reverse=True)

    def do_shape_action(self, event, target:SelectInfo, event_key):
        """Execute the target shape's action for the event_key.
        
        This is the default function that the do_action callable attribute is
        initialized to.
        """
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

    def on_key_press(self, event):
        pass
    
    def on_motion(self, event):
        # print(f"ifig.on_motion: mouse_down={self.is_mouse_down}")
        # if self.selected_shape is None:
        #     # print("selected_shape: None")
        #     pass
        # else:
        #     print(f"selected_shape: {self.selected_shape[0]}, "
        #           f"handle={self.selected_shape[1]}")
        # if self.is_mouse_down:
        #     self.mouse_down_count += 1
        # if self.mouse_down_count == 1:
        #     pass
        if self.selected_shape is None:
            artist_infos = self.find_artists_and_hilite(event)
            self.artist_infos = artist_infos
            # display_artist_and_event(f"on_motion ({len(artist_infos)})", 
            #                          event, artist_infos)
        else:
            if self.is_mouse_down:
                selection: SelectInfo = self.get_artist_for_handle(
                    self.selected_shape)
                # display_artist_and_event('on_drag', event, 
                #                          selection.artist)
                if selection is not None:
                    self.do_action(event, selection, 'drag')
            else:
                artist_infos = self.find_artists_and_hilite(event)
                self.artist_infos = artist_infos

    def on_press(self, event):
        self.save_do_scale_bounds = self.do_scale_bounds
        self.do_scale_bounds = False

        # print("on_press")
        # for i, a in enumerate(self.artists):
        #     if isinstance(a, SelectInfo):
        #         print(f"{i:2d} {artist_str(a.artist)}")
        #     else:
        #         print(f"{i:2d} {artist_str(a)}")

        artist_infos = self.find_artists_and_hilite(event)
        self.artist_infos = artist_infos
        selection: SelectInfo = (artist_infos[0]
                           if len(artist_infos) > 0 else None)

        self.is_mouse_down = True
        if selection is not None:
            self.selected_shape = (selection.artist.shape,
                                   selection.info)
            # print(f"selected_shape: {self.selected_shape[0]}, "
            #       f"handle={self.selected_shape[1]}")
            self.do_action(event, selection, 'press')
        else:
            self.selected_shape = None
        # print(f"selected_shape: {self.selected_shape[0]}, "
        #       f"handle={self.selected_shape[1]}")

        # print("selection.artist: "+artist_str(selection.artist))
        # for i, a in enumerate(self.artists):
        #     if isinstance(a, SelectInfo):
        #         print(f"{i:2d} {artist_str(a.artist)}")
        #     else:
        #         print(f"{i:2d} {artist_str(a)}")
        # if selection is not None:
        #     display_artist_and_event('on_press', event, selection.artist)

    def on_release(self, event):
        'on release we reset the press data'
        if self.selected_shape is not None:
            selection: SelectInfo = self.get_artist_for_handle(
                self.selected_shape)
            logger.debug("on_release")
            # print("on_release")
            if selection is not None:
                # display_artist_and_event('on_release', event, 
                #                         selection.artist)
                self.do_action(event, selection, 'release')
            self.selected_shape = None
        self.do_scale_bounds = self.save_do_scale_bounds
        # print(f"selected_shape: {self.selected_shape[0]}, "
        #       f"handle={self.selected_shape[1]}")
        self.is_mouse_down = False
        self.action_complete()

    def find_artists_and_hilite(self, event) -> list[SelectInfo]:
        """ identify and hilite artists at the mouse event location"""
        artist_infos = self.find_artists_at_location(event)
        self.hilited_artists = update_artist_hiliting(
            self.hilited_artists, artist_infos)
        return artist_infos


def update_artist_hiliting(hilited_artists, new_artists):
    """ manage artist hiliting in response to mouse event"""
    hilited_artist_set = {a.artist for a in hilited_artists}
    new_artists_set = {a.artist for a in new_artists}
    to_unhilite = list(hilited_artist_set - new_artists_set)
    for a in to_unhilite:
        try:
            a.unhighlight(a)
            a.figure.canvas.draw()
        except Exception as e:
            pass
    to_hilite = list(new_artists_set - hilited_artist_set)
    for a in to_hilite:
        a.highlight(a)
        a.figure.canvas.draw()
    return new_artists


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
            fig.ax, on_release, useblit=False,
            button=[1, 3],  # don't use middle button
            minspanx=5, minspany=5, spancoords='pixels', props=rectprops,
            interactive=False)


def enter_polyline(point_filter=None, constrain_to_wedge=True, **inputs):
    """ Graphical input of a 2d polyline """
    fig = inputs.pop('figure')
    canvas = fig.canvas
    ax = fig.ax
    polyline = np.array([])
    line = None
    cur_node = None
    pt0 = None
    pt2 = None
    do_on_finished = inputs.pop('do_on_finished', None)

    useblit = inputs.pop('useblit', True)
    background = fig.canvas.copy_from_bbox(fig.bbox)

    def init_poly(init_pt):
        nonlocal pt0
        pt0 = init_pt
        poly0 = np.array([init_pt, np.copy(init_pt)])
        poly0_data = poly0.T
        x, y = poly0_data[0], poly0_data[1]
        line0 = lines.Line2D(x, y, color=fig._rgb['foreground'],
                            marker='o',
                            animated=True)
        ax.add_line(line0)
        return poly0, line0

    def do_constrain_to_wedge(input_pt):
        """ keep the input point inside the wedge of adjacent points """
        projected_point_on_radial_line = misc_math.projected_point_on_radial_line
        if pt0 is not None:
            x_prod0 = input_pt[0]*pt0[1] - pt0[0]*input_pt[1]
            if x_prod0 < 0:
                # pin to boundary
                output_pt = projected_point_on_radial_line(input_pt, pt0)
                return output_pt

        if pt2 is not None:
            x_prod2 = input_pt[0]*pt2[1] - pt2[0]*input_pt[1]
            if x_prod2 > 0:
                # pin to boundary
                output_pt = projected_point_on_radial_line(input_pt, pt2)
                return output_pt

        return input_pt

    def update():
        nonlocal line, background
        if useblit:
            canvas.restore_region(background)
            ax.draw_artist(line)
            canvas.blit(fig.bbox)
        else:
            canvas.draw_idle()

    def draw_polyline(poly):
        """ draw the polyline """
        line.set_data(poly.T)
        ax.draw_artist(line)
        update()

    def add_point(pt):
        nonlocal polyline, pt0
        pt0 = polyline[-1]
        polyline = np.append(polyline, np.array([pt]), axis=0)
        return polyline

    def constrain_point(pt):
        if point_filter is not None:
            pt = point_filter(pt)
        if constrain_to_wedge:
            pt = do_constrain_to_wedge(pt)
        return pt

    def remove_last_point():
        nonlocal polyline, pt0
        polyline = np.delete(polyline, -1, axis=0)
        pt0 = polyline[-1]
        return polyline

    def on_mouse_move(event):
        """Callback for mouse movements."""
        nonlocal polyline, line, background

        if (event.inaxes is None):
            return

        if len(polyline) > 0:
            if event.xdata is not None and event.ydata is not None:
                event_data = np.array([event.xdata, event.ydata])
                event_data = constrain_point(event_data)
                polyline[-1] = event_data
    
                draw_polyline(polyline)

    def on_button_press(event):
        """Callback for mouse button presses."""
        nonlocal polyline, line

        if (event.inaxes is None or event.button != 1):
            return

        if event.xdata is not None and event.ydata is not None:
            event_data = np.array([event.xdata, event.ydata])
            event_data = constrain_point(event_data)

            if len(polyline) == 0:
                polyline, line = init_poly(event_data)
            else:
                polyline = add_point(event_data)
            draw_polyline(polyline)

    def on_button_release(event):
        """Callback for mouse button releases."""
        if (event.button != 1):
            return

    def on_key_press(event):
        """Callback for key presses."""
        if not event.inaxes:
            return

        event_key = event.key
        if event_key == 'escape':
            # we're done
            remove_last_point()  # potential new point not needed, delete
            fig.disconnect_events()
            fig.connect_events(saved_events)
            fig.action_complete()

        elif event_key == 'backspace' or event_key == 'backspace':
            # remove the most recently entered point
            remove_last_point()
            draw_polyline(polyline)

        elif event_key == 'shift':
            # reserved for add to selection
            pass
        elif event_key == 'control':
            # perhaps some sort of snap mode?
            pass
        canvas.draw()

    def on_finished_clb():
        """Callback for mouse movements."""
        nonlocal polyline, line, pt0
        pt0 = None
        draw_polyline(polyline)
        if do_on_finished is not None:
            do_on_finished(polyline, line)

    action_dict = {'button_press_event': on_button_press,
                   'button_release_event': on_button_release,
                   'key_press_event': on_key_press,
                   'motion_notify_event': on_mouse_move}

    saved_events = fig.disconnect_events()
    fig.connect_events(action_dict=action_dict)
    fig.on_finished = on_finished_clb


def snap_to_grid_fct(grid_space: float):
    def snap_to_grid(pt):
        return grid_space * np.around(pt/grid_space)
    return snap_to_grid
