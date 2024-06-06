#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
"""
.. Created on Wed Oct 16 14:20:49 2019

.. codeauthor: Michael J. Hayford
"""
import logging
import math
import numpy as np
from copy import deepcopy

from rayoptics.gui.util import bbox_from_poly, fit_data_range
from rayoptics.gui.actions import ReplaceGlassAction

import rayoptics.optical.model_constants as mc
from rayoptics.util.rgb2mpl import rgb2mpl
from rayoptics.util.misc_math import (normalize, distance_sqr_2d,
                                      projected_point_on_radial_line)
from rayoptics.util.line_intersection import get_intersect
from rayoptics.util import misc_math
from rayoptics.util import colors

from rayoptics.gui.util import calc_render_color_for_material

from rayoptics.parax import paraxialdesign


# --- color management
dgm_lw = {
    'data': 2,
    'data_hilite': 3,
    'line': 0.5,
    'hilite': 2,
    'guide': 1,
    'edge': 0.5,
    'node': 6,
    'node_hilite': 8,
    'node_edge': 0.75,
    }


def light_or_dark(is_dark=True):
    accent = colors.accent_colors(is_dark)
    fb = colors.foreground_background(is_dark)
    rgb = {
        'node': fb['foreground'],
        'edge': fb['foreground'],
        # 'node': accent['cyan'],
        # 'edge': accent['cyan'],
        'slide': accent['blue'],
        'object_image': accent['magenta'],
        'stop': accent['magenta'],
        'conj_line': accent['orange'],
        'shift': fb['foreground'],
        'barrel': accent['green'],
        }
    return {**rgb, **fb}


class Diagram():
    """ class for paraxial ray rendering/editing """

    def __init__(self, opt_model, parax_model, parax_model_key, dgm_type,
                 seq_start=1, do_barrel_constraint=False, barrel_constraint=1.0,
                 label='paraxial', bend_or_gap='bend', do_node_annotation=False,
                 is_dark=True):
        self.label = label
        self.opt_model = opt_model
        self.parax_model = parax_model

        self.dgm_rgb = light_or_dark(is_dark=is_dark)

        self.dgm_type = dgm_type
        self.setup_dgm_type(dgm_type)

        self.do_barrel_constraint = do_barrel_constraint
        self.barrel_constraint_radius = barrel_constraint

        self.do_node_annotation = do_node_annotation

        self.bend_or_gap = bend_or_gap

        self.active_layer = parax_model_key

        self.command_inputs = {}

    def setup_dgm_type(self, dgm_type):
        if dgm_type == 'ht':
            self.type_sel = mc.ht
        elif dgm_type == 'slp':
            self.type_sel = mc.slp

        pm = self.parax_model
        self._apply_data = pm.update_composite_node_fct(self.type_sel)

    def get_label(self):
        return self.label

    def sync_light_or_dark(self, is_dark):
        self.dgm_rgb = light_or_dark(is_dark)

    def set_active_layer(self, layer_key):
        opm = self.opt_model
        prx = paraxialdesign.update_diagram_for_key(opm, layer_key,
                                                    self.type_sel)
        prx_key, prx_model = prx
        self.active_layer = prx_key
        self.parax_model = prx_model
        self._apply_data = prx_model.update_composite_node_fct(self.type_sel)
        return self

    def update_data(self, fig, **kwargs):
        parax_model = self.parax_model

        self.shape = self.render_shape()
        if len(self.shape) > 0:
            self.shape_bbox = bbox_from_poly(self.shape)
            build = kwargs.get('build', 'rebuild')

            num_nodes = len(parax_model.sys)
            if self.dgm_type == 'slp':
                num_nodes -= 1

            if build == 'rebuild':
                self.node_list = []
                self.edge_list = []

            node_list_len = len(self.node_list)
            if num_nodes < node_list_len:
                del self.node_list[num_nodes:]
                del self.edge_list[num_nodes-1:]
            elif num_nodes > node_list_len:
                for i in range(node_list_len, num_nodes):
                    self.node_list.append(DiagramNode(self, i))
                for i in range(node_list_len, num_nodes-1):
                    self.edge_list.append(DiagramEdge(self, i))

            self.node_bbox = fig.update_patches(self.node_list)
            self.edge_bbox = fig.update_patches(self.edge_list)

            if build == 'rebuild':
                self.object_shift = ConjugateLine(self, 'object_image')
            self.object_shift_bbox = fig.update_patches([self.object_shift])

            if self.opt_model.seq_model.stop_surface is None:
                # stop shift conjugate line is only editable for floating stop
                if build == 'rebuild':
                    self.stop_shift = ConjugateLine(self, 'stop')
                self.stop_shift_bbox = fig.update_patches([self.stop_shift])

            if self.do_barrel_constraint:
                if build == 'rebuild':
                    self.barrel_constraint = BarrelConstraint(self)
                self.barrel_bbox = fig.update_patches([self.barrel_constraint])
            
            if self.do_node_annotation:
                pass
        else:
            self.shape_bbox = np.array([[-1, -1], [1, 1]])

        return self.shape_bbox

    def apply_data(self, node, vertex):
        self._apply_data(node, vertex)
        self.opt_model.parax_model.paraxial_lens_to_seq_model()

    def assign_object_to_node(self, node, new_node, factory, **kwargs):
        parax_model = self.parax_model
        inputs = parax_model.assign_object_to_node(node, new_node, 
                                                   self.type_sel, 
                                                   factory, **kwargs)
        return inputs

    def register_commands(self, *args, **inputs):
        fig = inputs['figure']
        gui_fct = inputs.pop('gui_fct', None)
        if gui_fct is not None:
            gui_fct(*args, **inputs)
        inputs.pop('figure')
        self.command_inputs = deepcopy(inputs)

        def do_command_action(event, target, event_key):
            nonlocal fig
            if target is not None:
                shape, handle = target.artist.shape
                try:
                    # print(type(shape).__name__, shape.node, handle, event_key)
                    handle_action_obj = shape.actions[handle]
                    if isinstance(handle_action_obj, dict):
                        handle_action_obj[event_key](fig, event)
                    else:
                        handle_action_obj.actions[event_key](fig, event)
                except KeyError:
                    pass
        fig.do_action = do_command_action

    def register_add_replace_element(self, *args, **inputs):
        fig = inputs.pop('figure')
        self.command_inputs = deepcopy(inputs)
        action_obj = AddReplaceElementAction(self, **inputs)

        def do_command_action(event, target, event_key):
            nonlocal fig
            shape, handle = target.artist.shape
            if isinstance(shape, DiagramNode) or \
               isinstance(shape, DiagramEdge):
                try:
                    action_obj.actions[event_key](fig, event, shape)
                except KeyError:
                    pass
        fig.do_action = do_command_action

    def render_shape(self):
        """ render the diagram into a shape list """
        parax_model = self.parax_model
        ts = self.type_sel
        shape = []
        for x, y in zip(parax_model.pr, parax_model.ax):
            x_guard = misc_math.infinity_guard(x[ts])
            y_guard = misc_math.infinity_guard(y[ts])
            shape.append([x_guard, y_guard])
        shape = np.array(shape)
        return shape

    def update_diagram_from_shape(self, shape):
        ''' use the shape list to update the paraxial model '''
        parax_model = self.parax_model
        for x, y, shp in zip(parax_model.pr, parax_model.ax, shape):
            x[self.type_sel] = shp[0]
            y[self.type_sel] = shp[1]

    def fit_axis_limits(self):
        ''' define diagram axis limits as the extent of the shape polygon '''
        if len(self.shape) > 0:
            x_min, x_max = fit_data_range([x[0] for x in self.shape])
            y_min, y_max = fit_data_range([x[1] for x in self.shape])
        else:
            x_min, x_max = -1, 1
            y_min, y_max = -1, 1
        return np.array([[x_min, y_min], [x_max, y_max]])


def compute_slide_line(shape, node, imode):
    """ compute a constraint line to keep the overall length of the airspaces
        surrounding `node` constant
    """
    def distance_sqr(pt):
        return pt[0]**2 + pt[1]**2
    num_nodes = len(shape)
    if node > 0 and node < num_nodes-1:
        pt0 = shape[node-1]
        pt1 = shape[node]
        pt2 = shape[node+1]
        if imode == 'transmit':
            # if transmitting, constrain movement of node to a line parallel
            #  to the line between the previous and next nodes.
            origin2line = misc_math.perpendicular_from_origin(pt0, pt2)
            pt2line = misc_math.perpendicular_to_line(pt1, pt0, pt2)
            scale_factor = (origin2line - pt2line)/origin2line
            new_pt0 = scale_factor*pt0
            new_pt2 = scale_factor*pt2
            return new_pt0, new_pt2
        elif imode == 'reflect':
            # if reflecting, constrain movement to a radial line from the
            #  origin through the selected node
            origin = np.array([0., 0.])
            dist_pt0 = distance_sqr(pt0)
            dist_pt1 = distance_sqr(pt1)
            dist_pt2 = distance_sqr(pt2)
            # figure out how long a line to draw
            dist = dist_pt0 if dist_pt0 > dist_pt1 else dist_pt1
            dist = dist_pt2 if dist_pt2 > dist else dist
            scale_factor = math.sqrt(dist/dist_pt1)
            new_pt1 = scale_factor*pt1
            return origin, new_pt1
    return None


def constrain_to_line_action(pt0, pt2):

    def constrain_to_line(input_pt):
        """ constrain the input point to the line pt0 to pt2 """
        output_pt = misc_math.projected_point_on_line(input_pt, pt0, pt2)
        return output_pt

    return constrain_to_line


class DiagramNode():
    def __init__(self, diagram, idx):
        self.diagram = diagram
        self.node = idx
        self.select_pt = None
        self.move_direction = None
        self.handles = {}
        self.actions = self.handle_actions()

    def listobj_str(self):
        o_str = f"node{str(self.node):>3}"
        return o_str

    def update_shape(self, view):
        diagram = self.diagram
        # set alpha to 25% -> #40
        bkgrnd_rbga = diagram.dgm_rgb['background1'] + '40'
        opt_model = diagram.opt_model
        shape = diagram.shape
        dgm_rgb = diagram.dgm_rgb
        if self.diagram.dgm_type == 'ht':
            hilite_kwargs = {
                'color': dgm_rgb['node'],
                'markersize': dgm_lw['node_hilite'],
                }
            self.handles['shape'] = (shape[self.node], 'vertex',
                                     {'linestyle': '',
                                      'linewidth': dgm_lw['data'],
                                      'marker': 'o',
                                      'markersize': dgm_lw['node'],
                                      'pickradius': 6,
                                      'color': dgm_rgb['node'],
                                      'hilite': dgm_rgb['hilite'],
                                      'zorder': 3.})
    
            # define the "constant spacing" or "slide" constraint
            if view.enable_slide:
                sys = self.diagram.parax_model.sys
                slide_pts = compute_slide_line(shape, self.node,
                                               sys[self.node][mc.rmd])
                if slide_pts is not None:
                    seg = [*slide_pts]
    
                    hilite_kwargs = {
                        'color': dgm_rgb['slide'],
                        'linewidth': dgm_lw['guide'],
                        'linestyle': '-'
                        }
                    self.handles['slide'] = (seg, 'polyline',
                                             {'linestyle': '--',
                                              'linewidth': dgm_lw['guide'],
                                              'pickradius': 6,
                                              'color': dgm_rgb['slide'],
                                              'hilite': hilite_kwargs,
                                              'zorder': 2.5})
        elif self.diagram.dgm_type == 'slp':
            hilite_kwargs = {
                'color': dgm_rgb['node'],
                'markersize': dgm_lw['node_hilite'],
                }

            gap, z_dir = diagram.parax_model.get_gap_for_node(self.node, 'slp')
            e_node = opt_model.part_tree.parent_node((gap, z_dir),
                                                     '#element#airgap')
            e = e_node.id if e_node else None
            marker_color = bkgrnd_rbga
            if e and '#airgap' not in e_node.tag:
                marker_color = calc_render_color_for_material(gap.medium)
                marker_color = rgb2mpl(marker_color)
            self.handles['shape'] = (shape[self.node], 'vertex',
                                     {'linestyle': '',
                                      'linewidth': dgm_lw['data'],
                                      'marker': 'o',
                                      'markersize': dgm_lw['node'],
                                      'markeredgewidth': dgm_lw['node_edge'],
                                      'markeredgecolor': dgm_rgb['node'],
                                      'pickradius': 6,
                                      'markerfacecolor': marker_color,
                                      'hilite': dgm_rgb['hilite'],
                                      'zorder': 3.})
        return view.create_patches(self.handles)

    def get_label(self):
        return 'node' + str(self.node)

    def handle_actions(self):
        actions = {}
        wedge_constraint = True if self.diagram.dgm_type == 'ht' else False
        actions['shape'] = EditNodeAction(self,
                                          constrain_to_wedge=wedge_constraint)
        slide_filter = None
        sys = self.diagram.parax_model.sys
        slide_pts = compute_slide_line(self.diagram.shape, self.node,
                                       sys[self.node][mc.rmd])
        if slide_pts is not None:
            slide_filter = constrain_to_line_action(*slide_pts)
        actions['slide'] = EditNodeAction(self, filter=slide_filter,
                                          constrain_to_wedge=wedge_constraint)
        return actions


class DiagramEdge():
    def __init__(self, diagram, idx):
        self.diagram = diagram
        self.node = idx
        self.select_pt = None
        self.move_direction = None
        self.handles = {}
        self.actions = self.handle_actions()

    def listobj_str(self):
        o_str = f"edge{str(self.node):>3}"
        return o_str

    def update_shape(self, view):
        shape = self.diagram.shape
        dgm_rgb = self.diagram.dgm_rgb
        edge_poly = shape[self.node:self.node+2]
        self.handles['shape'] = (edge_poly, 'polyline',
                                 {'pickradius': 6,
                                  'linewidth': dgm_lw['data'],
                                  'hilite_linewidth': dgm_lw['data_hilite'],
                                  'color': dgm_rgb['edge'],
                                  'hilite': dgm_rgb['hilite'],
                                  'zorder': 2.})
        area_poly = [[0, 0]]
        area_poly.extend(edge_poly)
        fill_color = self.render_color()
        self.handles['area'] = (area_poly, 'polygon',
                                {'fill_color': fill_color,
                                 'zorder': 1.})

        return view.create_patches(self.handles)

    def render_color(self):
        diagram = self.diagram
        opt_model = diagram.opt_model
        # set alpha to 25% -> #40
        bkgrnd_rbga = diagram.dgm_rgb['background1'] + '40'
        if diagram.dgm_type == 'ht':
            gap, z_dir = diagram.parax_model.get_gap_for_node(self.node, 'ht')
            e_node = opt_model.part_tree.parent_node((gap, z_dir),
                                                     '#element#airgap')
            e = e_node.id if e_node else None
    
            if e:
                if '#airgap' in e_node.tag:
                    return bkgrnd_rbga
                else:
                    return calc_render_color_for_material(gap.medium)
            else:
                # single surface element, like mirror or thinlens, use airgap
                return bkgrnd_rbga
        elif diagram.dgm_type == 'slp':
            return bkgrnd_rbga

    def get_label(self):
        return 'edge' + str(self.node)

    def handle_actions(self):
        actions = {}
        actions['shape'] = EditLensAction(self)
        actions['area'] = EditAreaAction(self)
        return actions


# --- Editing actions
class BarrelConstraint():
    def __init__(self, diagram):
        self.diagram = diagram
        self.handles = {}
        self.actions = self.handle_actions()

    def listobj_str(self):
        o_str = "BarrelConstraint"
        return o_str

    def update_shape(self, view):
        dgm_rgb = self.diagram.dgm_rgb
        barrel_radius = self.diagram.barrel_constraint_radius
        diamond = []
        diamond.append([0., barrel_radius])
        diamond.append([barrel_radius, 0.])
        diamond.append([0., -barrel_radius])
        diamond.append([-barrel_radius, 0.])
        diamond.append([0., barrel_radius])
        diamond = np.array(diamond)
        self.handles['shape'] = (diamond, 'polyline',
                                 {'color': dgm_rgb['barrel'],
                                  'linewidth': dgm_lw['guide'],
                                  'zorder': 1.})
        square = []
        square.append([ barrel_radius,  barrel_radius])
        square.append([-barrel_radius,  barrel_radius])
        square.append([-barrel_radius, -barrel_radius])
        square.append([ barrel_radius, -barrel_radius])
        square.append([ barrel_radius,  barrel_radius])
        square = np.array(square)
        self.handles['square'] = (square, 'polyline',
                                  {'color': dgm_rgb['barrel'],
                                   'linewidth': dgm_lw['guide'],
                                   'zorder': 1.})

        return view.create_patches(self.handles)

    def get_label(self):
        return 'barrel constraint'

    def handle_actions(self):
        actions = {}
        return actions


class ConjugateLine():
    def __init__(self, diagram, line_type):
        self.diagram = diagram
        self.line_type = line_type
        self.sys_orig = []
        self.shape_orig = []
        self.handles = {}
        self.actions = self.handle_actions()

    def listobj_str(self):
        o_str = f"ConjugateLine: {self.line_type}"
        return o_str

    def update_shape(self, view):
        dgm_rgb = self.diagram.dgm_rgb
        shape_bbox = self.diagram.shape_bbox
        line = []
        if self.line_type == 'stop':
            ht = shape_bbox[1][1] - shape_bbox[0][1]
            line.append([0., -2*ht])
            line.append([0., 2*ht])
            color = dgm_rgb['stop']
        elif self.line_type == 'object_image':
            wid = shape_bbox[1][0] - shape_bbox[0][0]
            line.append([-2*wid, 0.])
            line.append([2*wid, 0.])
            color = dgm_rgb['object_image']
        self.handles['shape'] = (line, 'polyline',
                                 {'color': dgm_rgb['foreground'],
                                  'hilite': color,
                                  'zorder': 1.})

        if len(self.shape_orig) > 0:
            conj_line = []
            if self.line_type == 'stop':
                lwr, upr = view.ax.get_ybound()
                ht = upr - lwr
                conj_line.append([-self.k*ht, -ht])
                conj_line.append([self.k*ht, ht])
            elif self.line_type == 'object_image':
                lwr, upr = view.ax.get_xbound()
                wid = upr - lwr
                conj_line.append([-wid, -self.k*wid])
                conj_line.append([wid, self.k*wid])
            self.handles['conj_line'] = (conj_line, 'polyline',
                                         {'color': dgm_rgb['conj_line'],
                                          'linewidth': dgm_lw['guide'],
                                          'zorder': 1.})
            self.handles['shift'] = (self.shape_orig, 'polyline',
                                     {'color': dgm_rgb['shift'],
                                      'linewidth': 1.5,
                                      'zorder': 1.})

        return view.create_patches(self.handles)

    def get_label(self):
        if self.line_type == 'stop':
            return 'stop shift line'
        elif self.line_type == 'object_image':
            return 'object shift line'
        else:
            return ''

    def edit_conjugate_line_actions(self):
        dgm = self.diagram
        pm = dgm.parax_model

        def calculate_slope(x, y):
            """ x=ybar, y=y  """
            if self.line_type == 'stop':
                k = x/y
                return k, np.array([[1, -k], [0, 1]]).T
            elif self.line_type == 'object_image':
                k = y/x
                return k, np.array([[1, 0], [-k, 1]]).T
            else:
                k = 0
                return 0, np.array([[1, 0], [0, 1]])

        def apply_data(event_data):
            nonlocal self
            k, mat = calculate_slope(event_data[0], event_data[1])
            pm.apply_conjugate_shift(self.shape_orig, k, mat,
                                     self.line_type)
            self.k = k

        def on_select(fig, event):
            nonlocal self
            self.sys_orig = deepcopy(pm.sys)
            self.shape_orig = deepcopy(dgm.shape)

        def on_edit(fig, event):
            nonlocal self
            if event.xdata is not None and event.ydata is not None:
                event_data = np.array([event.xdata, event.ydata])
                apply_data(event_data)
                fig.refresh_gui(build='update', src_model=pm)

        def on_release(fig, event):
            nonlocal self
            event_data = np.array([event.xdata, event.ydata])
            apply_data(event_data)
            self.sys_orig = []
            self.shape_orig = []
            fig.refresh_gui(build='rebuild', src_model=pm)

        actions = {}
        actions['press'] = on_select
        actions['drag'] = on_edit
        actions['release'] = on_release
        return actions

    def handle_actions(self):
        actions = {}
        actions['shape'] = self.edit_conjugate_line_actions()
        return actions


class EditNodeAction():
    """ Action to move a diagram node, using an input pt """

    def __init__(self, dgm_node: DiagramNode, 
                 filter=None, constrain_to_wedge=True):
        if not isinstance(dgm_node, DiagramNode):
            raise TypeError
        diagram = dgm_node.diagram
        parax_model = diagram.parax_model
        self.cur_node = None
        self.pt0 = None
        self.pt2 = None
        self.filter = filter
        self.constrain_to_wedge = constrain_to_wedge

        def point_on_line(pt1, pt2, t):
            d = pt2 - pt1
            return pt1 + t*d

        def do_constrain_to_wedge(input_pt):
            """ keep the input point inside the wedge of adjacent points """

            if self.pt0 is not None:
                x_prod0 = input_pt[0]*self.pt0[1] - self.pt0[0]*input_pt[1]
                if x_prod0 < 0:
                    # pin to boundary
                    output_pt = projected_point_on_radial_line(input_pt,
                                                               self.pt0)
                    return output_pt

            if self.pt2 is not None:
                x_prod2 = input_pt[0]*self.pt2[1] - self.pt2[0]*input_pt[1]
                if x_prod2 > 0:
                    # pin to boundary
                    output_pt = projected_point_on_radial_line(input_pt,
                                                               self.pt2)
                    return output_pt

            return input_pt

        def on_select(fig, event):
            nonlocal diagram
            # we don't allow points to move onto their adjacent neighbors. Use
            #  a buffer amount when constraining to the wedge
            buffer_fraction = 0.0025
            self.cur_node = dgm_node.node
            pt1 = diagram.shape[self.cur_node]
            if self.cur_node == 0:
                pt2 = diagram.shape[self.cur_node+1]
                self.pt2 = point_on_line(pt1, pt2, 1-buffer_fraction)
                self.pt0 = None
            elif self.cur_node == len(diagram.shape)-1:
                pt0 = diagram.shape[self.cur_node-1]
                self.pt0 = point_on_line(pt0, pt1, buffer_fraction)
                self.pt2 = None
            else:
                pt0 = diagram.shape[self.cur_node-1]
                self.pt0 = point_on_line(pt0, pt1, buffer_fraction)
                pt2 = diagram.shape[self.cur_node+1]
                self.pt2 = point_on_line(pt1, pt2, 1-buffer_fraction)

        def on_edit(fig, event):
            if event.xdata is not None and event.ydata is not None:
                event_data = np.array([event.xdata, event.ydata])
                if self.filter:
                    event_data = self.filter(event_data)
                if self.constrain_to_wedge:
                    event_data = do_constrain_to_wedge(event_data)
                diagram.apply_data(self.cur_node, event_data)
                fig.refresh_gui(build='update', src_model=parax_model)

        def on_release(fig, event):
            if event.xdata is not None and event.ydata is not None:
                event_data = np.array([event.xdata, event.ydata])
                if self.filter:
                    event_data = self.filter(event_data)
                if self.constrain_to_wedge:
                    event_data = do_constrain_to_wedge(event_data)
                diagram.apply_data(self.cur_node, event_data)
                fig.refresh_gui(build='update', src_model=parax_model)
                self.cur_node = None

        self.actions = {}
        self.actions['drag'] = on_edit
        self.actions['press'] = on_select
        self.actions['release'] = on_release


class EditLensAction():
    """ Action for diagram edge, using an input pt

    This is a simple wrapper class to choose the correct action, i.e. bending
    or thickness change, depending on the UI setting.
    """

    def __init__(self, dgm_edge: DiagramEdge):
        if not isinstance(dgm_edge, DiagramEdge):
            raise TypeError
        diagram = dgm_edge.diagram
        actions = {}
        actions['gap'] = EditThicknessAction(dgm_edge)
        actions['bend'] = EditBendingAction(dgm_edge)

        def create_dispatch_action(event_key):
            def dispatch_action(fig, event):
                actions[diagram.bend_or_gap].actions[event_key](fig, event)
            return dispatch_action

        self.actions = {}
        self.actions['press'] = create_dispatch_action('press')
        self.actions['drag'] = create_dispatch_action('drag')
        self.actions['release'] = create_dispatch_action('release')


class EditAreaAction():
    """ Action for diagram area, placeholder for now

    This is a simple wrapper class to choose the correct action, i.e. bending
    or thickness change, depending on the UI setting.
    """

    def __init__(self, dgm_edge: DiagramEdge):
        # actions['power'] = EditPowerAction(dgm_edge)
        # actions['bend'] = EditBendingAction(dgm_edge)

        def create_dispatch_action(event_key):
            def dispatch_action(fig, event):
                pass
                # actions[diagram.bend_or_gap].actions[event_key](fig, event)
            return dispatch_action

        self.actions = {}
        self.actions['press'] = create_dispatch_action('press')
        self.actions['drag'] = create_dispatch_action('drag')
        self.actions['release'] = create_dispatch_action('release')


class EditThicknessAction():
    """ Action to move a diagram edge, using an input pt

    The movement is constrained to be parallel to the original edge. By doing
    this the power and bending of the element remains constant, while the
    element thickness changes. Movement of the edge is limited to keep
    the thickness greater than zero and not to interfere with adjacent spaces.
    """

    def __init__(self, dgm_edge: DiagramEdge):
        if not isinstance(dgm_edge, DiagramEdge):
            raise TypeError
        diagram = dgm_edge.diagram
        parax_model = diagram.parax_model
        self.node = None
        self.bundle = None

        def on_select(fig, event):
            nonlocal diagram
            shape = diagram.shape
            self.node = node = dgm_edge.node
            if node > 0 and node < (len(shape)-2):
                # get the virtual vertex of the combined element surfaces
                vertex = np.array(get_intersect(shape[node-1], shape[node],
                                                shape[node+1], shape[node+2]))
                # get direction cosines for the edge
                edge = normalize(shape[node+1] - shape[node])
                # construct perpendicular to the edge. use this to define a
                # range for allowed inputs
                perp_edge = np.array([edge[1], -edge[0]])
                self.bundle = vertex, edge, perp_edge
                # measure distances along the perpendicular thru the vertex
                #  and project the 2 outer nodes and the first vertex of the
                #  edge
                pp0 = np.dot(shape[node-1]-vertex, perp_edge)
                pp1 = np.dot(shape[node]-vertex, perp_edge)
                pp3 = np.dot(shape[node+2]-vertex, perp_edge)
                # use the first edge vertex to calculate an allowed input range
                if pp1 > 0:
                    self.pp_lim = 0, min(i for i in (pp0, pp3) if i > 0)
                else:
                    self.pp_lim = max(i for i in (pp0, pp3) if i < 0), 0

        def on_edit(fig, event):
            buffer = 0.0025
            nonlocal diagram
            shape = diagram.shape
            node = self.node
            if node > 0 and node < (len(shape)-2):
                if event.xdata is not None and event.ydata is not None:
                    inpt = np.array([event.xdata, event.ydata])
                    vertex, edge, perp_edge = self.bundle
                    # project input pt onto perp_edge
                    ppi = np.dot(inpt - vertex, perp_edge)
                    # constrain the input point within the range, if needed
                    if ppi < self.pp_lim[0]:
                        inpt = vertex + (1+buffer)*self.pp_lim[0]*perp_edge
                    elif ppi > self.pp_lim[1]:
                        inpt = vertex + (1-buffer)*self.pp_lim[1]*perp_edge

                    # compute new edge vertices from intersection of adjacent
                    #  edges and line shifted parallel to the initial edge
                    edge_pt = inpt + edge
                    pt1 = np.array(get_intersect(shape[node-1], vertex,
                                                 inpt, edge_pt))
                    diagram.apply_data(self.node, pt1)
                    pt2 = np.array(get_intersect(vertex, shape[node+2],
                                                 inpt, edge_pt))
                    diagram.apply_data(self.node+1, pt2)
                    fig.refresh_gui(build='update', src_model=parax_model)

        def on_release(fig, event):
            if event.xdata is not None and event.ydata is not None:
                event_data = np.array([event.xdata, event.ydata])
                fig.refresh_gui(build='rebuild', src_model=parax_model)
                self.node = None
                self.bundle = None

        self.actions = {}
        self.actions['drag'] = on_edit
        self.actions['press'] = on_select
        self.actions['release'] = on_release


class EditBendingAction():
    """ Action to bend the lens element for diagram edge, using an input pt.

    The movement is constrained to be along the object ray for the lens if the
    input point is closer to the leading node of the edge. Otherwise the
    movement is constrained to be along the image ray. The unconstrained point
    is solved to keep the element thickness constant and maintain the
    object-image properties of the lens.
    """

    def __init__(self, dgm_edge: DiagramEdge):
        if not isinstance(dgm_edge, DiagramEdge):
            raise TypeError
        diagram = dgm_edge.diagram
        pm = diagram.parax_model

        self.node = None
        self.bundle = None
        self.filter = None

        def cross_prod(pt1, pt2):
            return pt1[1]*pt2[0] - pt2[1]*pt1[0]

        def calc_coef_fct(vertex, iNode, dir_inpt, oNode, dir_out):
            nonlocal pm
            tau_factor = pm.sys[self.node][mc.tau]*pm.opt_inv
            constrain_to_line = constrain_to_line_action(vertex,
                                                         vertex+dir_inpt)
            def calc_t(inpt):
                pt = constrain_to_line(inpt)
                if iNode < oNode:
                    arg1 = pt, vertex
                    arg2 = pt, dir_out
                else:
                    arg1 = vertex, pt
                    arg2 = dir_out, pt
                t = (tau_factor - cross_prod(*arg1))/(cross_prod(*arg2))
                return (iNode, pt), (oNode, vertex + t*dir_out)
            return calc_t

        def on_select(fig, event):
            nonlocal diagram
            if event.xdata is None or event.ydata is None:
                return
            shape = diagram.shape
            self.node = node = dgm_edge.node
            if node > 0 and node < (len(shape)-2):
                inpt = np.array([event.xdata, event.ydata])
                # get the virtual vertex of the combined element surfaces
                vertex = np.array(get_intersect(shape[node-1], shape[node],
                                                shape[node+1], shape[node+2]))
                edge_dir_01 = normalize(shape[node] - shape[node-1])
                edge_dir_23 = normalize(shape[node+2] - shape[node+1])
                # which node is closer to the input point?
                pt1_dist = distance_sqr_2d(shape[node], inpt)
                pt2_dist = distance_sqr_2d(shape[node+1], inpt)
                if pt1_dist < pt2_dist:
                    self.filter = calc_coef_fct(vertex, node, edge_dir_01,
                                                node+1, edge_dir_23)
                else:
                    self.filter = calc_coef_fct(vertex, node+1, edge_dir_23,
                                                node, edge_dir_01)
                # get direction cosines for the edge
                edge = normalize(shape[node+1] - shape[node])
                # construct perpendicular to the edge. use this to define a
                # range for allowed inputs
                perp_edge = np.array([edge[1], -edge[0]])
                self.bundle = (vertex, edge, perp_edge, edge_dir_01,
                               edge_dir_23)

        def on_edit(fig, event):
            nonlocal diagram
            shape = diagram.shape
            node = self.node
            if node > 0 and node < (len(shape)-2):
                if event.xdata is not None and event.ydata is not None:
                    inpt = np.array([event.xdata, event.ydata])
                    inp, out = self.filter(inpt)

                    diagram.apply_data(inp[0], inp[1])
                    diagram.apply_data(out[0], out[1])
                    fig.refresh_gui(build='update', src_model=pm)

        def on_release(fig, event):
            if event.xdata is not None and event.ydata is not None:
                event_data = np.array([event.xdata, event.ydata])
                fig.refresh_gui(build='rebuild', src_model=pm)
                self.node = None
                self.bundle = None
                self.filter = None

        self.actions = {}
        self.actions['drag'] = on_edit
        self.actions['press'] = on_select
        self.actions['release'] = on_release


class AddReplaceElementAction():
    ''' insert or replace a node with a chunk from a factory fct

    The do_command_action fct registered for this operation passes the shape
    being operated upon; these can be:

        - DiagramEdge: insert/add the chunk returned by the factory fct
        - DiagramNode: replace the selected node with the factory fct output

    Inserting is done by splitting the corresponding gap in two. A new gap
    and an AirGap element are tacked on to the chunk returned from the factory
    fct.
    Replacing is done when a DiagramNode is selected. The gaps surrounding the
    node are retained, and modified as needed to accomodate the chunk.
    '''
    from rayoptics.mpl.interactivefigure import SelectInfo

    def __init__(self, diagram, **kwargs):
        seq_model = diagram.opt_model.seq_model
        parax_model = diagram.parax_model
        self.command_inputs = dict(kwargs)
        self.cur_node = None
        self.init_inputs = None

        def on_press_add_point(fig, event, shape):
            # if we don't have factory functions, skip the command
            nonlocal self, diagram
            if isinstance(shape, DiagramEdge):
                if 'node_init' in self.command_inputs and \
                   'factory' in self.command_inputs:

                    node = shape.node
                    event_data = np.array([event.xdata, event.ydata])
                    interact = self.command_inputs['interact_mode']
                    sys_data = [1.0, interact]
                    new_node = parax_model.add_node(node, event_data, 
                                                    diagram.type_sel, sys_data)
                    self.cur_node = cur_node = node, new_node
                    do_insert = True

            elif isinstance(shape, DiagramNode):
                if 'factory' in self.command_inputs:
                    # replacing a node with a chunk only requires recording
                    # what chunk corresponds to the current node. There is
                    # no drag action
                    node = shape.node
                    # self.init_inputs = parax_model.get_object_for_node(node)
                    # node_init = self.command_inputs['node_init']
                    self.cur_node = cur_node = node, node
                    do_insert = False

            else:
                return

            factory = self.command_inputs['factory']
            self.init_inputs = diagram.assign_object_to_node(
                *cur_node, factory, insert=do_insert, do_update=False)
            parax_model.paraxial_lens_to_seq_model()
            parax_model.update_parax_to_dgms()
            fig.refresh_gui(build='rebuild', src_model=parax_model)

        def on_drag_add_point(fig, event, shape):
            nonlocal self, diagram
            if self.cur_node is not None and isinstance(shape, DiagramEdge):
                event_data = np.array([event.xdata, event.ydata])
                # node, new_node = self.cur_node
                # diagram.apply_data(new_node, event_data)
                # fig.refresh_gui(build='update', src_model=parax_model)

        def on_release_add_point(fig, event, shape):
            nonlocal self, diagram
            if self.cur_node is not None:
                factory = self.command_inputs['factory']
                # if factory and node_init fcts are the same, we're done;
                # always call factory fct for a node
                # if factory != self.command_inputs['node_init'] or \
                #               isinstance(shape, DiagramNode):
                #     inputs = diagram.assign_object_to_node(
                #                   *self.cur_node, 
                #                   factory)
                # fig.refresh_gui(build='rebuild', src_model=parax_model)
            self.cur_node = None

        self.actions = {}
        self.actions['press'] = on_press_add_point
        self.actions['drag'] = on_drag_add_point
        self.actions['release'] = on_release_add_point


class AddElementAction():
    ''' insert or replace a node with a chunk from a factory fct

    The do_command_action fct registered for this operation passes the shape
    being operated upon; these can be:

        - DiagramEdge: insert/add the chunk returned by the factory fct
        - DiagramNode: replace the selected node with the factory fct output

    Inserting is done by splitting the corresponding gap in two. A new gap
    and an AirGap element are tacked on to the chunk returned from the factory
    fct.
    Replacing is done when a DiagramNode is selected. The gaps surrounding the
    node are retained, and modified as needed to accomodate the chunk.
    '''

    def __init__(self, dgm_edge: DiagramEdge, **kwargs):
        if not isinstance(dgm_edge, DiagramEdge):
            raise TypeError
        diagram = dgm_edge.diagram
        parax_model = diagram.parax_model

        self.command_inputs = dict(kwargs)
        self.cur_node = None
        self.init_inputs = None

        def on_press_add_point(fig, event):
            # if we don't have factory functions, skip the command
            nonlocal self, diagram
            if 'node_init' in self.command_inputs and \
                'factory' in self.command_inputs:

                node = dgm_edge.node
                event_data = np.array([event.xdata, event.ydata])
                interact = self.command_inputs['interact_mode']
                sys_data = [1.0, interact]
                new_node = parax_model.add_node(node, event_data, 
                                                diagram.type_sel, sys_data)
                self.cur_node = node, new_node

                # create a node for editing during the drag action
                #  'node_init' will currently be a thinlens or a mirror
                node_init = self.command_inputs['node_init']
                self.init_inputs = diagram.assign_object_to_node(
                                        *self.cur_node, node_init,
                                        insert=True, do_update=False)
                parax_model.paraxial_lens_to_seq_model()
                fig.refresh_gui(build='update', src_model=parax_model)

        def on_drag_add_point(fig, event):
            nonlocal self, diagram
            if self.cur_node is not None:
                event_data = np.array([event.xdata, event.ydata])
                diagram.apply_data(self.cur_node, event_data)
                fig.refresh_gui(build='update', src_model=parax_model)

        def on_release_add_point(fig, event):
            nonlocal self, diagram
            if self.cur_node is not None:
                factory = self.command_inputs['factory']
                # if factory and node_init fcts are the same, we're done;
                # always call factory fct for a node
                if factory != self.command_inputs['node_init']:
                    inputs = diagram.assign_object_to_node(
                        *self.cur_node, factory)
                fig.refresh_gui(build='update', src_model=parax_model)
            self.cur_node = None

        self.actions = {}
        self.actions['press'] = on_press_add_point
        self.actions['drag'] = on_drag_add_point
        self.actions['release'] = on_release_add_point


class ReplaceElementAction():
    ''' insert or replace a node with a chunk from a factory fct

    The do_command_action fct registered for this operation passes the shape
    being operated upon; these can be:

        - DiagramEdge: insert/add the chunk returned by the factory fct
        - DiagramNode: replace the selected node with the factory fct output

    Inserting is done by splitting the corresponding gap in two. A new gap
    and an AirGap element are tacked on to the chunk returned from the factory
    fct.
    Replacing is done when a DiagramNode is selected. The gaps surrounding the
    node are retained, and modified as needed to accomodate the chunk.
    '''

    def __init__(self, dgm_node: DiagramNode, **kwargs):
        if not isinstance(dgm_node, DiagramNode):
            raise TypeError
        diagram = dgm_node.diagram
        parax_model = diagram.parax_model
        self.command_inputs = dict(kwargs)
        self.cur_node = None
        self.init_inputs = None

        def on_press_add_point(fig, event, shape):
            # if we don't have factory functions, skip the command
            nonlocal self, diagram
            if isinstance(shape, DiagramNode):
                if 'factory' in self.command_inputs:
                    # replacing a node with a chunk only requires recording
                    # what chunk corresponds to the current node. There is
                    # no drag action
                    node = shape.node
                    self.cur_node = node, None
                    self.init_inputs = parax_model.get_object_for_node(node)

        def on_drag_add_point(fig, event, shape):
            nonlocal self, diagram
            pass

        def on_release_add_point(fig, event, shape):
            nonlocal self, diagram
            if self.cur_node is not None:
                factory = self.command_inputs['factory']
                # if factory and node_init fcts are the same, we're done;
                # always call factory fct for a node
                if factory != self.command_inputs['node_init'] or \
                              isinstance(shape, DiagramNode):
                    inputs = diagram.assign_object_to_node(
                                  *self.cur_node, 
                                  factory)
                fig.refresh_gui(build='update', src_model=parax_model)
            self.cur_node = None

        self.actions = {}
        self.actions['press'] = on_press_add_point
        self.actions['drag'] = on_drag_add_point
        self.actions['release'] = on_release_add_point


class GlassDropAction():

    def dragEnterEvent(self, canvas_fig, event):
        def glass_target_filter(artist):
            shape, handle = artist.shape
            if handle == 'area' and 'area' in shape.actions:
                sm = shape.diagram.opt_model.seq_model
                gap = sm.gaps[shape.node]
                if gap.medium.name().lower() != 'air':
                    return False
            return True
        canvas_fig.figure.artist_filter = glass_target_filter

    def dragMoveEvent(self, canvas_fig, event):
        canvas_fig.mouseMoveEvent(event)
        
    def dragLeaveEvent(self, canvas_fig, event):
        canvas_fig.figure.artist_filter = None

    def dropEvent(self, canvas_fig, event):
        dropped_it = False
        fig = canvas_fig.figure
        if fig.artist_infos is not None:
            for selection in fig.artist_infos:
                shape, handle = selection.artist.shape
                if handle == 'area' and 'area' in shape.actions:
                    sm = shape.diagram.opt_model.seq_model
                    gap = sm.gaps[shape.node]
                    action = ReplaceGlassAction(gap, update=False)
                    action(fig, event)
                    pm = shape.diagram.parax_model
                    pm.update_rindex(shape.node)
                    pm.paraxial_lens_to_seq_model()
                    fig.refresh_gui()
                    dropped_it = True
        canvas_fig.figure.artist_filter = None
        return dropped_it


from rayoptics.mpl.interactivefigure import SelectInfo, display_artist_and_event


def add_shape(point_filter=None, constrain_to_wedge=True, **inputs):
    """ Graphical input of a 2d polyline """
    fig = inputs.pop('figure')
    do_on_finished = inputs.pop('do_on_finished', None)
    
    def do_shape_action(fig, event, target:"SelectInfo", event_key):
        """Execute the target shape's action for the event_key.
        """
        if target is not None:
            shape, handle = target.artist.shape
            try:
                handle_action_obj = shape.actions[handle]
                if isinstance(handle_action_obj, dict):
                    handle_action_obj[event_key](fig, event)
                else:
                    handle_action_obj.actions[event_key](fig, event)
            except KeyError:
                pass
    
    def do_selected_shape_action(fig, event, selected_shape, event_key):
        """Execute the target shape's action for the event_key.
        """
        shape, handle = selected_shape[0]
        try:
            handle_action_obj = shape.actions[handle]
            handle_action_obj.actions[event_key](fig, event)
        except KeyError:
            pass

    def do_command_action(fig, event, target, event_key):
        shape, handle = target.artist.shape
        if isinstance(shape, DiagramNode) or \
            isinstance(shape, DiagramEdge):
            try:
                action_obj.actions[event_key](fig, event, shape)
            except KeyError:
                pass

    def on_mouse_move(event):
        if fig.selected_shape is None:
            artist_infos = fig.find_artists_and_hilite(event)
            fig.artist_infos = artist_infos
        else:
            if fig.is_mouse_down:
                selected_artist = fig.get_artist_for_handle(fig.selected_shape)
                do_shape_action(fig, event, selected_artist, 'drag')
            else:
                artist_infos = fig.find_artists_and_hilite(event)
                fig.artist_infos = artist_infos

    def on_button_press(event):
        fig.save_do_scale_bounds = fig.do_scale_bounds
        fig.do_scale_bounds = False

        artist_infos = fig.find_artists_and_hilite(event)
        fig.artist_infos = artist_infos
        selected_artist = (artist_infos[0]
                           if len(artist_infos) > 0 else None)

        fig.is_mouse_down = True
        if fig.selected_shape is None and selected_artist is not None:
            fig.selected_shape = (selected_artist.artist.shape, 
                                  selected_artist.info)
            do_shape_action(fig, event, selected_artist, 'press')
        elif fig.selected_shape is not None:
            do_selected_shape_action(fig, event, fig.selected_shape, 'press')
        else:
            fig.selected_shape = None

    def on_button_release(event):
        'on release we reset the press data'
        if fig.selected_shape is not None:
            do_selected_shape_action(fig, event, fig.selected_shape, 'release')
            # selected_artist = fig.get_artist_for_handle(fig.selected_shape)
            # do_shape_action(fig, event, selected_artist, 'release')
            fig.selected_shape = None
        fig.do_scale_bounds = fig.save_do_scale_bounds
        fig.is_mouse_down = False
        fig.action_complete()

    def on_key_press(event):
        """Callback for key presses."""
        if not event.inaxes:
            return

    def on_finished_clb():
        """Callback for mouse movements."""
        nonlocal saved_events
        if do_on_finished is not None:
            do_on_finished()
        fig.disconnect_events()
        fig.connect_events(saved_events)

    action_dict = {'button_press_event': on_button_press,
                   'button_release_event': on_button_release,
                   'key_press_event': on_key_press,
                   'motion_notify_event': on_mouse_move}

    saved_events = fig.disconnect_events()
    fig.connect_events(action_dict=action_dict)
    fig.on_finished = on_finished_clb

def context_menu_actions():
    """ return list of action_desc, action_handle, action """
    menu_defs = []
    menu_defs.append(("Add node", 'add', AddElementAction))
    menu_defs.append(("Bend element", 'bend', EditBendingAction))
    menu_defs.append(("Thicken element", 'gap', EditThicknessAction))
    return menu_defs
