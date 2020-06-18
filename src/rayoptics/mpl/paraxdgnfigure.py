#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
"""

.. Created on Mon Apr  2 19:20:27 2018

.. codeauthor: Michael J. Hayford

.. deprecated:: 0.4.5
"""

import numpy as np
from matplotlib.figure import Figure

from rayoptics.optical.model_constants import ht, slp, aoi
from rayoptics.optical.elements import (create_thinlens, create_mirror,
                                        create_lens)

from rayoptics.util.misc_math import distance_sqr_2d, perpendicular_distance_2d


def create_parax_design_commands(fig):
    cmds = []
    # Select an existing point
    cmds.append(('Select', (fig.select_point, (), {})))
    # Add thin lens
    cmds.append(('Add Thin Lens',
                 (fig.add_point, (), {'factory': create_thinlens})))
    # Add lens
    cmds.append(('Add Lens', (fig.add_point, (), {'factory': create_lens})))
    # Add mirror
    cmds.append(('Add Mirror',
                 (fig.add_point, (), {'factory': create_mirror})))

    return cmds


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

    def process_hit_location(self, event, hit_list, prt=False):
        xdata, ydata = self.line.get_data()
        # create a list of vertex coordinates in the model coordinate system
        line_hits = [[xdata[i], ydata[i]] for i in hit_list]
        # transform the list into display coordinates
        dsp_hits = self.line.axes.transData.transform(line_hits)
        hit_pt = [event.x, event.y]

        hit_vertex = None
        hit_edge = None
        min_hit_dist = 1e10
        min_perp_dist = 1e10
        for i, pt in enumerate(dsp_hits):
            # consider hit point as a vertex
            hit_dist = distance_sqr_2d(pt, hit_pt)
            if hit_dist < min_hit_dist:
                min_hit_dist = hit_dist
                if hit_dist < self.pick_radius_sqr:
                    hit_vertex = hit_list[i]

            # consider hit point on an edge
            nxt = hit_list[i]+1
            try:
                next_pt = np.array([xdata[nxt], ydata[nxt]])
            except IndexError:
                continue
            else:
                next_dsp_pt = self.line.axes.transData.transform(next_pt)
                perp_dist = perpendicular_distance_2d(pt, hit_pt, next_dsp_pt)
                if perp_dist < min_perp_dist:
                    min_perp_dist = perp_dist
#                    if perp_dist < self.pick_radius_sqr:
#                    print("perp_dist", hit_list[i], perp_dist)
                    hit_edge = hit_list[i]
        v = None
        if hit_vertex is not None:
            v = (hit_vertex, event.xdata, event.ydata,
                 xdata[hit_vertex], ydata[hit_vertex], 's')
#            print("vertex selected %d: event x=%f y=%f data x=%f y=%f" %
#                  (hit_vertex,
#                   event.xdata, event.ydata,
#                   xdata[hit_vertex], ydata[hit_vertex]))
        e = None
        if hit_edge is not None:
            h = hit_edge
#            print("edge selected", hit_list, min_hit_dist,
#                  xdata[h:h+2], ydata[h:h+2])
            e = (hit_edge, event.xdata, event.ydata,
                 xdata[h:h+2], ydata[h:h+2], '')

        return (v, e)

    def on_press(self, event):
        # on button press we will see if the mouse is over us and store
        #  some data
#        print("eline.on_press")
        if event.inaxes != self.line.axes:
            return

        contains, attrd = self.line.contains(event)
        if not contains:
            return

        hit_list = attrd['ind']
        v, e = self.process_hit_location(event, hit_list)

        self.press = v, e

    def on_motion(self, event):
        'on motion we will highlight a vertex of edge if the mouse is over it '
        if event.inaxes != self.line.axes:
            return

        contains, props = self.line.contains(event)
        if not contains:
            self.markers.set_data([], [])
            self.line.figure.canvas.draw()
            return

#        print("eline.on_motion")
        hit_list = props['ind']
        v, e = self.process_hit_location(event, hit_list)

        if v:
            hit_vertex, x_hit, y_hit, x_data, y_data, mkr = v
        elif e:
            hit_vertex, x_hit, y_hit, x_data, y_data, mkr = e

        if v or e:
            self.markers.set_xdata(x_data)
            self.markers.set_ydata(y_data)
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

    def __init__(self, opt_model, refresh_gui, dgm_type, **kwargs):
        self.actions = {}
        self.opt_model = opt_model
        self.parax_model = opt_model.parax_model
        self.refresh_gui = refresh_gui
        self.setup_dgm_type(dgm_type)
        self.skip_build = False

        Figure.__init__(self, **kwargs)

        self.vertex = None
        self.update_data()
        self.select_point()

    def setup_dgm_type(self, dgm_type):
        if dgm_type == 'ht':
            self.type_sel = ht
            self.data_slice = slice(0, None)
#            self.data_slice = slice(1, None)
            self.x_label = r'$\overline{y}$'
            self.y_label = 'y'
            self.apply_data = self.parax_model.apply_ht_dgm_data
            self.header = r'$y-\overline{y}$ Diagram'
        elif dgm_type == 'slp':
            self.type_sel = slp
            self.data_slice = slice(0, -1)
            self.x_label = r'$\overline{\omega}$'
            self.y_label = r'$\omega$'
            self.apply_data = self.parax_model.apply_slope_dgm_data
            self.header = r'$\omega-\overline{\omega}$ Diagram'

    def refresh(self, **kwargs):
        self.update_data(**kwargs)
        self.plot()

    def update_data(self, **kwargs):
        if not self.skip_build:
            self.parax_model.build_lens()
        self.skip_build = False
        return self

    def fit_data_range(self, x_data, margin=0.05, range_trunc=0.25):
        x_min = min(0., min(x_data))
        x_max = max(0., max(x_data))
        x_range = x_max - x_min
        if x_range != 0.0 and len(x_data) > 2:
            x1_min = min(0., min(x_data[1:]))
            x1_max = max(0., max(x_data[1:]))
            x1_range = x1_max - x1_min
            if abs(x1_range/x_range) < range_trunc:
                x_min = x1_min
                x_max = x1_max
                x_range = x1_range

        if x_range > 0.:
            x_margin = margin*x_range
        else:
            x_margin = 0.01
        return x_min-x_margin, x_max+x_margin

    def plot(self):
        self.clf()
        self.ax = self.add_subplot(1, 1, 1)
        self.ax.axvline(0, c='black', lw=1)
        self.ax.axhline(0, c='black', lw=1)
        self.ax.set_title(self.header, pad=10.0, fontsize=18)

        x_data = []
        y_data = []
        for x, y in zip(self.parax_model.pr,
                        self.parax_model.ax):
            x_data.append(x[self.type_sel])
            y_data.append(y[self.type_sel])

        x_min, x_max = self.fit_data_range(x_data)
        y_min, y_max = self.fit_data_range(y_data)

        self.line, = self.ax.plot(x_data,
                                  y_data,
                                  marker='s', picker=6)
        self.ax.set_xlabel(self.x_label)
        self.ax.set_ylabel(self.y_label)
        self.ax.set_xlim(x_min, x_max)
        self.ax.set_ylim(y_min, y_max)

        xl, xr = self.ax.get_xlim()
        yl, yu = self.ax.get_ylim()
        xlb, xrb = self.ax.get_xbound()
        ylb, yub = self.ax.get_ybound()

        self.eline = EditableLine(self.line)
        self.eline.connect()
        self.canvas.mpl_connect('pick_event', self.on_pick)
        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('button_release_event', self.on_release)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.canvas.draw()

        return self

    def select_point(self, **kwargs):
        def on_press_select_point(pdfig, press):
            hit_vertex, x_hit, y_hit, x_data, y_data, mkr = press
            pdfig.vertex = hit_vertex + pdfig.data_slice.start
        self.actions['press'] = on_press_select_point

        def on_drag_select_point(pdfig, event):
            pdfig.apply_data(pdfig.vertex, (event.xdata, event.ydata))
            pdfig.parax_model.paraxial_lens_to_seq_model()
            pdfig.skip_build = True
            pdfig.refresh_gui()
        self.actions['drag'] = on_drag_select_point

        def on_release_select_point(pdfig, event):
            pdfig.apply_data(pdfig.vertex, (event.xdata, event.ydata))
            pdfig.parax_model.paraxial_lens_to_seq_model()
            pdfig.skip_build = True
            pdfig.refresh_gui()
            pdfig.vertex = None
        self.actions['release'] = on_release_select_point

        self.move_action = None

    def add_point(self, factory=None, **kwargs):
        def on_press_add_point(pdfig, press):
            hit_vertex, x_hit, y_hit, x_data, y_data, mkr = press
            new_vertex = hit_vertex + pdfig.data_slice.start
#            print("add_point:", new_vertex)
            pdfig.parax_model.add_object(new_vertex, (x_hit, y_hit),
                                         pdfig.type_sel, factory, 'transmit')
            new_vertex += 1
            pdfig.parax_model.paraxial_lens_to_seq_model()
            pdfig.skip_build = True
            pdfig.refresh_gui()
        self.actions['press'] = on_press_add_point

    def on_press(self, event):
#        print("on_press")
        hit, props = self.line.contains(event)
        if hit:
            if self.eline.press is None:
                self.eline.on_press(event)
            if self.eline.press:
                v, e = self.eline.press
                if v:
                    self.actions['press'](self, v)
                    hit_vertex, x_hit, y_hit, x_data, y_data, mkr = v
                    self.vertex = hit_vertex + self.data_slice.start
#                    print("vertex selected", hit_vertex)
                elif e:
                    hit_edge, x_hit, y_hit, x_data, y_data, mkr = e
                    self.actions['press'](self, e)
#                    print("edge selected", e[0])

#            print('on_press', event.button, event.x, event.y,
#                  event.xdata, event.ydata, event.key,
#                  len(hit_list), hit_list)

    def on_motion(self, event):
        if self.vertex:
            self.actions['drag'](self, event)

    def on_release(self, event):
        'on release we reset the press data'
        if self.vertex:
            self.actions['release'](self, event)

    def on_pick(self, event):
        line = event.artist
        me = event.mouseevent
        xdata, ydata = line.get_data()
        ind = event.ind
        id = line.get_gid()
#        print("on_pick", id, ind, xdata[ind], ydata[ind], me.name, me.x, me.y,
#              me.button, me.key, me.xdata, me.ydata, me.dblclick)
