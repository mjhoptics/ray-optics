#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
"""
.. Created on Wed Oct 16 14:20:49 2019

.. codeauthor: Michael J. Hayford
"""

import numpy as np

from rayoptics.gui.util import GUIHandle, bbox_from_poly

from rayoptics.optical.model_constants import ht, slp
from rayoptics.util.rgb2mpl import rgb2mpl


class Diagram():
    """ class for paraxial ray rendering/editing """
    def __init__(self, opt_model, dgm_type, seq_start=1,
                 label='paraxial'):
        self.label = label
        self.opt_model = opt_model

        self.setup_dgm_type(dgm_type)

        if opt_model.specsheet.conjugate_type == 'finite':
            self.seq_start = 0
        elif opt_model.specsheet.conjugate_type == 'infinite':
            self.seq_start = 1
        else:
            self.seq_start = seq_start

        self.cur_node = None

    def setup_dgm_type(self, dgm_type):
        parax_model = self.opt_model.parax_model
        if dgm_type == 'ht':
            self.type_sel = ht
            self.apply_data = parax_model.apply_ht_dgm_data
        elif dgm_type == 'slp':
            self.type_sel = slp
            self.apply_data = parax_model.apply_slope_dgm_data

    def get_label(self):
        return self.label

    def register_commands(self, *args, **kwargs):
        fig = kwargs.pop('figure')
        cmd_actions = kwargs.pop('cmd_actions')
        actions = cmd_actions(**kwargs)

        def do_command_action(event, target, event_key):
            nonlocal fig, actions
            if target is not None:
                shape, handle = target.artist.shape
                try:
                    action = actions[event_key]
                    action(fig, handle, event, target.info)
                except KeyError:
                    pass
        fig.do_action = do_command_action

    def render_shape(self):
        # render the diagram into a shape list
        parax_model = self.opt_model.parax_model
        shape = []
        for x, y in zip(parax_model.pr, parax_model.ax):
            shape.append([x[self.type_sel], y[self.type_sel]])
        return shape

    def edit_diagram_actions(self):
        parax_model = self.opt_model.parax_model
        actions = {}

        def on_select_point(fig, handle, event, info):
            self.cur_node = self.seq_start
            if 'ind' in info:
                self.cur_node += info['ind'][0]
        actions['press'] = on_select_point

        def on_edit_point(fig, handle, event, info):
            event_data = np.array([event.xdata, event.ydata])
            self.apply_data(self.cur_node, event_data)
            parax_model.paraxial_lens_to_seq_model()
            fig.refresh_gui()
        actions['drag'] = on_edit_point

        def on_release_point(fig, handle, event, info):
            event_data = np.array([event.xdata, event.ydata])
            self.apply_data(self.cur_node, event_data)
            parax_model.paraxial_lens_to_seq_model()
            fig.refresh_gui()
            self.cur_node = None
        actions['release'] = on_release_point

        return actions

    def add_element_actions(self, factory=None, node_init=None, **kwargs):
        parax_model = self.opt_model.parax_model
        actions = {}

        def on_press_add_point(fig, handle, event, info):
            self.cur_node = self.seq_start
            if 'ind' in info:
                self.cur_node += info['ind'][0]
#            print("add_point: press", self.cur_node)
            event_data = np.array([event.xdata, event.ydata])
            parax_model.add_node(self.cur_node, event_data, self.type_sel)
            self.cur_node += 1
            parax_model.assign_object_to_node(self.cur_node, node_init)
            parax_model.paraxial_lens_to_seq_model()
            fig.skip_build = True
            fig.refresh_gui()
        actions['press'] = on_press_add_point

        def on_drag_add_point(fig, handle, event, info):
#            print("add_point: drag", self.cur_node)
            event_data = np.array([event.xdata, event.ydata])
            self.apply_data(self.cur_node, event_data)
            parax_model.paraxial_lens_to_seq_model()
            fig.refresh_gui()
        actions['drag'] = on_drag_add_point

        def on_release_add_point(fig, handle, event, info):
#            print("add_point: release", self.cur_node)
            parax_model.assign_object_to_node(self.cur_node, factory)
            parax_model.paraxial_lens_to_seq_model()
            fig.skip_build = True
            fig.refresh_gui()
            self.cur_node = None
        actions['release'] = on_release_add_point

        return actions


class DiagramNode():
    def __init__(self, diagram, idx):
        self.diagram = diagram
        self.node = idx
        self.select_pt = None
        self.move_direction = None
        self.handles = {}
#        self.actions = self.edit_diagram_actions()

    def update_shape(self, view):
        n_color = rgb2mpl([138, 43, 226])  # blueviolet
        shape = self.diagram.shape
        self.handles['shape'] = shape[self.node], 'vertex', {'linestyle': '',
                                                             'marker': 's',
                                                             'picker': 6,
                                                             'color': n_color,
                                                             'hilite': 'red'}
        gui_handles = {}
        for key, graphics_handle in self.handles.items():
            poly_data, poly_type, kwargs = graphics_handle
            poly = np.array(poly_data)
            if poly_type == 'vertex':
                priority = 2.5
                p = view.create_vertex(poly, zorder=priority, **kwargs)
            elif poly_type == 'polyline':
                priority = 2.
                p = view.create_polyline(poly, zorder=priority, **kwargs)
            elif poly_type == 'polygon':
                p = view.create_polygon(poly, self.render_color(),
                                        zorder=2.5, **kwargs)
            else:
                break
            if len(poly.shape) > 1:
                bbox = bbox_from_poly(poly)
            else:
                x = poly[0]
                y = poly[1]
                bbox = np.array([[x, y], [x, y]])
            gui_handles[key] = GUIHandle(p, bbox)
        return gui_handles

    def render_color(self):
        e = self.diagram.opt_model.ele_model.elements[self.node]
        return e.render_color

    def get_label(self):
        return 'node' + str(self.node)

    @property
    def actions(self):
        return self.diagram.actions

    def edit_diagram_actions(self):
        parax_model = self.diagram.opt_model.parax_model
        actions = {}

        def on_select_point(fig, handle, event, info):
            self.cur_node = self.node
        actions['press'] = on_select_point

        def on_edit_point(fig, handle, event, info):
            event_data = np.array([event.xdata, event.ydata])
            self.diagram.apply_data(self.cur_node, event_data)
            parax_model.paraxial_lens_to_seq_model()
            fig.refresh_gui()
        actions['drag'] = on_edit_point

        def on_release_point(fig, handle, event, info):
            event_data = np.array([event.xdata, event.ydata])
            self.diagram.apply_data(self.cur_node, event_data)
            parax_model.paraxial_lens_to_seq_model()
            fig.refresh_gui()
            self.cur_node = None
        actions['release'] = on_release_point

        return actions


class DiagramEdge():
    def __init__(self, diagram, idx):
        self.diagram = diagram
        self.edge = idx
        self.select_pt = None
        self.move_direction = None
        self.handles = {}
#        self.actions = self.add_element_actions()

    def update_shape(self, view):
        shape = self.diagram.shape
        edge_poly = shape[self.edge:self.edge+2]
        self.handles['shape'] = edge_poly, 'polyline', {'picker': 6,
                                                        'hilite': 'red'}
        area_poly = [[0, 0]]
        area_poly.extend(edge_poly)
        fill_color = self.render_color()
        self.handles['area'] = area_poly, 'polygon', {'fill_color': fill_color}

        gui_handles = {}
        for key, graphics_handle in self.handles.items():
            poly_data, poly_type, kwargs = graphics_handle
            poly = np.array(poly_data)
            if poly_type == 'polygon':
                p = view.create_polygon(poly, self.render_color(),
                                        zorder=2.5, **kwargs)
            elif poly_type == 'polyline':
                priority = 2.
                p = view.create_polyline(poly, zorder=priority, **kwargs)
            else:
                break
            gui_handles[key] = GUIHandle(p, bbox_from_poly(poly))
        return gui_handles

    def render_color(self):
        if self.edge == 0:
            return (237, 243, 254, 64)  # light blue
        else:
            e = self.diagram.opt_model.ele_model.elements[self.edge-1]
            return e.render_color

    def get_label(self):
        return 'edge' + str(self.edge)

    @property
    def actions(self):
        return self.diagram.actions

    def add_element_actions(self, factory=None, node_init=None, **kwargs):
        parax_model = self.diagram.opt_model.parax_model
        actions = {}

        def on_press_add_point(fig, handle, event, info):
            self.cur_node = self.edge
#            print("add_point: press", self.cur_node)
            event_data = np.array([event.xdata, event.ydata])
            parax_model.add_node(self.cur_node, event_data,
                                 self.diagram.type_sel)
            self.cur_node += 1
            parax_model.assign_object_to_node(self.cur_node, node_init)
            fig.parax_model.paraxial_lens_to_seq_model()
            fig.skip_build = True
            fig.refresh_gui()
        actions['press'] = on_press_add_point

        def on_drag_add_point(fig, handle, event, info):
#            print("add_point: drag", self.cur_node)
            event_data = np.array([event.xdata, event.ydata])
            self.diagram.apply_data(self.cur_node, event_data)
            parax_model.paraxial_lens_to_seq_model()
            fig.refresh_gui()
        actions['drag'] = on_drag_add_point

        def on_release_add_point(fig, handle, event, info):
#            print("add_point: release", self.cur_node)
            parax_model.assign_object_to_node(self.cur_node, factory)
            parax_model.paraxial_lens_to_seq_model()
            fig.skip_build = True
            fig.refresh_gui()
            self.cur_node = None
        actions['release'] = on_release_add_point

        return actions


class EditNodeAction():
    """ Action to move a diagram node, using an input pt """
    def __init__(self, dgm, parax_model, dgm_type):
        self.parax_model = parax_model
        self.select_pt_lcl = None
        self.cur_value = None
        self.new_value = None
        self.actions = {}

        def on_select(fig, event):
            self.cur_node = self.dgm.seq_start
            if 'ind' in info:
                self.cur_node += info['ind'][0]
        self.actions['press'] = on_select

        def on_edit(fig, event, pt):
            event_data = np.array([event.xdata, event.ydata])
            self.diagram.apply_data(self.cur_node, event_data)
            self.parax_model.paraxial_lens_to_seq_model()
            fig.refresh_gui()
        self.actions['drag'] = on_edit

        def on_release(fig, event):
            self.cv_new = self.ele.reference_interface().profile_cv
            xsag = sag(self.cv_new, event.lcl_pt[1])
#            print('on_release: {:.3f} {:.3f} {:.5f}'
#                  .format(event.lcl_pt[0], xsag, self.cv_new))
            fig.refresh_gui()
        self.actions['release'] = on_release
