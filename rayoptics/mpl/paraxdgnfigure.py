#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
"""
Created on Mon Apr  2 19:20:27 2018

@author: Michael J. Hayford
"""

from enum import Enum, auto

from matplotlib.figure import Figure

import rayoptics.optical.paraxialdesign as pd
from rayoptics.optical.model_constants import ax, pr, lns
from rayoptics.optical.model_constants import ht, slp, aoi
from rayoptics.util.misc_math import distance_sqr_2d


class Dgm(Enum):
    ht = auto()
    slp = auto()


class EditableLine:
    def __init__(self, line):
        self.line = line
        self.press = None
        self.is_hilited = False
        self.markers, = self.line.axes.plot([], [], 'rs')
        self.pick_radius_sqr = self.line.get_pickradius()**2

    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.line.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.line.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.line.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

    def process_hit_location(self, event, hit_list):
        xdata, ydata = self.line.get_data()
        line_hits = [[xdata[i], ydata[i]] for i in hit_list]
        dsp_hits = self.line.axes.transData.transform(line_hits)
        hit_pt = [event.x, event.y]

        hit_vertex = None
        min_hit_dist = 1e10
        for i, pt in enumerate(dsp_hits):
            hit_dist = distance_sqr_2d(pt, hit_pt)
            if hit_dist < min_hit_dist:
                min_hit_dist = hit_dist
                if hit_dist < self.pick_radius_sqr:
                    hit_vertex = hit_list[i]
        if hit_vertex is None:
            h = hit_list[0]
#            print("edge selected", hit_list, min_hit_dist,
#                  xdata[h:h+2], ydata[h:h+2])
            return xdata[h:h+2], ydata[h:h+2], ''
        else:
            return xdata[hit_vertex], ydata[hit_vertex], 's'
#            print("vertex selected %d: event x=%f y=%f data x=%f y=%f" %
#                  (hit_vertex,
#                   event.xdata, event.ydata,
#                   xdata[hit_vertex], ydata[hit_vertex]))

    def on_press(self, event):
        # on button press we will see if the mouse is over us and store
        #  some data
        if event.inaxes != self.line.axes:
            return

        contains, attrd = self.line.contains(event)
        if not contains:
            return

        hit_list = attrd['ind']
        xd, yd, mkr = self.process_hit_location(event, hit_list)

        self.press = xd, yd, event.xdata, event.ydata

    def on_motion(self, event):
        'on motion we will highlight a vertex of edge if the mouse is over it '
        if event.inaxes != self.line.axes:
            return

        contains, props = self.line.contains(event)
        if not contains:
            self.markers.set_data([], [])
            self.line.figure.canvas.draw()
            return

        hit_list = props['ind']
        xd, yd, mkr = self.process_hit_location(event, hit_list)
        self.markers.set_xdata(xd)
        self.markers.set_ydata(yd)
        self.markers.set_marker(mkr)
        self.markers.set_linestyle('solid')
        self.line.figure.canvas.draw()

    def on_release(self, event):
        'on release we reset the press data'
        self.press = None
        self.line.figure.canvas.draw()

    def disconnect(self):
        'disconnect all the stored connection ids'
        self.line.figure.canvas.mpl_disconnect(self.cidpress)
        self.line.figure.canvas.mpl_disconnect(self.cidrelease)
        self.line.figure.canvas.mpl_disconnect(self.cidmotion)


class ParaxialDesignFigure(Figure):

    def __init__(self, seq_model, refresh_gui, dgm_type, **kwargs):
        self.seq_model = seq_model
        self.refresh_gui = refresh_gui
        self.setup_dgm_type(dgm_type)
        self.skip_build = False

        Figure.__init__(self, **kwargs)

        self.vertex = None
        self.update_data()

    def setup_dgm_type(self, dgm_type):
        if dgm_type == Dgm.ht:
            self.type_sel = ht
            self.data_slice = slice(1, None)
            self.x_label = r'$\overline{y}$'
            self.y_label = 'y'
            self.apply_data = pd.apply_ht_dgm_data
            self.header = r'$y-\overline{y}$ Diagram'
        elif dgm_type == Dgm.slp:
            self.type_sel = slp
            self.data_slice = slice(0, -1)
            self.x_label = r'$\overline{\omega}$'
            self.y_label = r'$\omega$'
            self.apply_data = pd.apply_slope_dgm_data
            self.header = r'$\omega-\overline{\omega}$ Diagram'

    def refresh(self):
        self.update_data()
        self.plot()

    def update_data(self):
        if not self.skip_build:
            self.lens = pd.build_lens(self.seq_model)
        self.skip_build = False
        return self

    def plot(self):
        self.clf()
        self.ax = self.add_subplot(1, 1, 1)
        self.ax.axvline(0, c='black', lw=1)
        self.ax.axhline(0, c='black', lw=1)
        self.ax.set_title(self.header, pad=10.0, fontsize=18)

        x_data = self.lens[pr][self.type_sel][self.data_slice]
        y_data = self.lens[ax][self.type_sel][self.data_slice]
        self.line, = self.ax.plot(x_data, y_data, marker='s', picker=6)
        self.ax.set_xlabel(self.x_label)
        self.ax.set_ylabel(self.y_label)

        self.eline = EditableLine(self.line)
        self.eline.connect()
        self.canvas.mpl_connect('pick_event', self.on_pick)
        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('button_release_event', self.on_release)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.canvas.draw()

        return self

    def on_press(self, event):
        hit, props = self.line.contains(event)
        if hit:
            hit_list = props['ind']
#            print("event.inaxes", event.inaxes)
            x_data = self.lens[pr][self.type_sel][self.data_slice]
            y_data = self.lens[ax][self.type_sel][self.data_slice]
            line_hits = [[x_data[i], y_data[i]] for i in hit_list]
            dsp_hits = self.ax.transData.transform(line_hits)
            hit_pt = [event.x, event.y]
            pick_radius_sqr = self.line.get_pickradius()**2
            hit_vertex = None
            min_hit_dist = 1e10
            for i, pt in enumerate(dsp_hits):
                hit_dist = distance_sqr_2d(pt, hit_pt)
#                print("distance_sqr", hit_list[i], hit_dist)
                if hit_dist < min_hit_dist:
                    min_hit_dist = hit_dist
                    if hit_dist < pick_radius_sqr:
                        hit_vertex = hit_list[i]
            if hit_vertex is None:
                pass
#                print("edge selected", hit_list, min_hit_dist)
            else:
                self.vertex = hit_vertex + self.data_slice.start
#                print("vertex selected", hit_vertex, min_hit_dist)
#            print('on_press', event.button, event.x, event.y,
#                  event.xdata, event.ydata, event.key,
#                  len(hit_list), hit_list)

    def on_motion(self, event):
        if self.vertex:
            self.lens[pr][self.type_sel][self.vertex] = event.xdata
            self.lens[ax][self.type_sel][self.vertex] = event.ydata
#            print("on_motion", self.vertex, event.xdata, event.ydata)
            self.apply_data(self.lens, self.vertex)
            self.seq_model.paraxial_lens_to_seq_model(self.lens)
            self.skip_build = True
            self.refresh_gui()

    def on_release(self, event):
        'on release we reset the press data'
        if self.vertex:
            self.lens[pr][self.type_sel][self.vertex] = event.xdata
            self.lens[ax][self.type_sel][self.vertex] = event.ydata
#            print("on_release", self.vertex, event.xdata, event.ydata)
            self.apply_data(self.lens, self.vertex)
            self.seq_model.paraxial_lens_to_seq_model(self.lens)
            self.skip_build = True
            self.refresh_gui()
            self.vertex = None

    def on_pick(self, event):
        line = event.artist
        me = event.mouseevent
        xdata, ydata = line.get_data()
        ind = event.ind
        id = line.get_gid()
#        print("on_pick", id, ind, xdata[ind], ydata[ind], me.name, me.x, me.y,
#              me.button, me.key, me.xdata, me.ydata, me.dblclick)
