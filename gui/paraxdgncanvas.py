#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
"""
Created on Mon Apr  2 19:20:27 2018

@author: Michael J. Hayford
"""

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QSizePolicy

from matplotlib.backends.backend_qt5agg \
     import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import optical.paraxialdesign as pd


def distance_sqr(pt0, pt1):
    return (pt0[0] - pt1[0])**2 + (pt0[1] - pt1[1])**2


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
            hit_dist = distance_sqr(pt, hit_pt)
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
#        print('event contains', self.line.xy)
        x0, y0 = self.line.xy
        self.press = x0, y0, event.xdata, event.ydata

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


class ParaxialDesignCanvas(FigureCanvas):

    def __init__(self, parent, seq_model, width=5, height=4, dpi=100):
        self.seq_model = seq_model
        self.lens = pd.build_lens(self.seq_model)

        self.fig = Figure(figsize=(width, height), dpi=dpi)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        # Next 2 lines are needed so that key press events are correctly
        #  passed with mouse events
        # https://github.com/matplotlib/matplotlib/issues/707/
        self.setFocusPolicy(Qt.ClickFocus)
        self.setFocus()

        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.plot()

    def update_plot(self):
        self.lens = pd.build_lens(self.seq_model)
        self.plot()

    def plot(self):
        self.fig.clf()
        self.axes = self.fig.add_subplot(1, 1, 1)
        self.axes.axvline(0, c='black', lw=1)
        self.axes.axhline(0, c='black', lw=1)

        ax_ht = self.lens[pd.ax][pd.ht][1:]
        pr_ht = self.lens[pd.pr][pd.ht][1:]
        self.line, = self.axes.plot(pr_ht, ax_ht, marker='s', picker=6)
        self.axes.set_xlabel(r'$\overline{y}$')
        self.axes.set_ylabel('y')
#        ax_nu = self.lens[pd.ax][pd.slp][:-1]
#        pr_nu = self.lens[pd.pr][pd.slp][:-1]
#        self.line, = self.axes.plot(pr_nu, ax_nu, marker='s', picker=5)
#        self.axes.set_xlabel(r'$\overline{\omega}$')
#        self.axes.set_ylabel(r'$\omega$')

        self.eline = EditableLine(self.line)
        self.eline.connect()
        self.fig.canvas.mpl_connect('pick_event', self.on_pick)
        self.fig.canvas.mpl_connect('button_press_event', self.on_press)

        self.draw()

    def on_press(self, event):
        hit, props = self.line.contains(event)
        if hit:
            hit_list = props['ind']
            print("event.inaxes", event.inaxes)
            ax_ht = self.lens[pd.ax][pd.ht][1:]
            pr_ht = self.lens[pd.pr][pd.ht][1:]
            line_hits = [[pr_ht[i], ax_ht[i]] for i in hit_list]
            dsp_hits = self.axes.transData.transform(line_hits)
            hit_pt = [event.x, event.y]
            pick_radius_sqr = self.line.get_pickradius()**2
            hit_vertex = None
            min_hit_dist = 1e10
            for i, pt in enumerate(dsp_hits):
                hit_dist = distance_sqr(pt, hit_pt)
#                print("distance_sqr", hit_list[i], hit_dist)
                if hit_dist < min_hit_dist:
                    min_hit_dist = hit_dist
                    if hit_dist < pick_radius_sqr:
                        hit_vertex = hit_list[i]
            if hit_vertex is None:
                print("edge selected", hit_list, min_hit_dist)
            else:
                print("vertex selected", hit_vertex, min_hit_dist)
#            print('on_press', event.button, event.x, event.y,
#                  event.xdata, event.ydata, event.key, len(hit_list), hit_list)

    def on_pick(self, event):
        line = event.artist
        me = event.mouseevent
        xdata, ydata = line.get_data()
        ind = event.ind
        id = line.get_gid()
#        print("on_pick", id, ind, xdata[ind], ydata[ind], me.name, me.x, me.y,
#              me.button, me.key, me.xdata, me.ydata, me.dblclick)
