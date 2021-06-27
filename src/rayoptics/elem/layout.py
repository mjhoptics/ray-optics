#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Interactive 2D lens picture

    The Lens Layout capability provides a 2D display of the optical model
    represented by **shapes**. Shapes contain `dict` attributes to manage the
    graphical rendering and editing actions associated with their parent
    objects. These attributes must include:

    **handles**: named graphic handles for different aspects of the parent object

    **actions**: functions for **press**, **drag**, and **release** actions

.. Created on Tue Sep 18 14:23:28 2018

.. codeauthor: Michael J. Hayford
"""

import numpy as np
from anytree import Node

from rayoptics.gui.util import (GUIHandle, transform_ray_seg, bbox_from_poly,
                                transform_poly, inv_transform_poly)

from rayoptics.elem import parttree
from rayoptics.elem import transform
from rayoptics.raytr.analyses import retrieve_ray, RayFan
from rayoptics.raytr.trace import (trace_boundary_rays_at_field,
                                   boundary_ray_dict, RaySeg)
from rayoptics.elem import elements as ele
import rayoptics.seq.medium as medium

import rayoptics.optical.model_constants as mc
from rayoptics.util.rgb2mpl import rgb2mpl
import rayoptics.gui.appcmds as cmds
from rayoptics.util import colors


def light_or_dark(is_dark=True):
    accent = colors.accent_colors(is_dark)
    fb = colors.foreground_background(is_dark)
    rgb = {
        'profile': fb['foreground'],
        'edge': fb['foreground'],
        'ct': accent['blue'],
        'ray': fb['foreground'],
        'rayfan_fill': [254, 197, 254, 64],  # magenta, 25%
        'ax_ray': rgb2mpl([138,  43, 226]),  # blueviolet
        'pr_ray': rgb2mpl([198, 113, 113]),  # sgi salmon
        }
    return {**rgb, **accent, **fb}


lo_rgb = {}

lo_lw = {
    'line': 1,
    'hilite': 2,
    'guide': 1.0,
    'edge': 0.5,
    }


def setup_shift_of_ray_bundle(seq_model, start_offset):
    """ compute transformation for rays "start_offset" from 1st surface

    Args:
        seq_model: the sequential model
        start_offset: z distance rays should start wrt first surface.
                      positive if to left of first surface
    Returns:
        transformation rotation and translation::
            (rot, t)
    """

    s1 = seq_model.ifcs[1]
    s0 = seq_model.ifcs[0]
    rot, t = transform.reverse_transform(s0, start_offset, s1)
    return rot, t


def shift_start_of_ray_bundle(start_bundle, ray_bundle, rot, t, cr_indx=0):
    """ modify ray_bundle so that rays begin "start_offset" from 1st surface

    Args:
        ray_bundle: list of rays in a bundle, i.e. all for one field.
                    ray_bundle[cr_indx] is the chief/central ray
        start_offset: z distance rays should start wrt first surface.
                      positive if to left of first surface
        rot: transformation rotation
        t: transformation translation
        cr_indx: index of the central ray in the bundle
    """

    # For the chief ray, use the input offset.
    ray, op_delta, wvl = ray_bundle[cr_indx]
    pt1_t = rot.dot(ray[1].p) - t
    dir0 = rot.dot(ray[0].d)
    dst = -pt1_t[2]/dir0[2]
    pt0 = pt1_t + dst*dir0
    ray0 = RaySeg(pt0, dir0, dst, ray[0].nrml)
    start_bundle[cr_indx] = ray0
    
    for ri, ray_pkg in enumerate(ray_bundle):
        ray, op_delta, wvl = ray_pkg
        b4_pt = rot.dot(ray[1].p) - t
        b4_dir = rot.dot(ray[0].d)
        if ri != cr_indx:
            # Calculate distance along ray to plane perpendicular to
            #  the chief ray.
            dst = -(b4_pt - pt0).dot(dir0)/b4_dir.dot(dir0)
            pt = b4_pt + dst*b4_dir
            ray0 = RaySeg(pt, b4_dir, dst, ray[0].nrml)
            start_bundle[ri] = ray0


def create_optical_element(opt_model, e):
    # if isinstance(e, ele.CementedElement):
    #     e_list = []
    #     for elemnt in e.ele_list:
    #         e_list.append(OpticalElement(opt_model, elemnt))
    # else:
    #     e_list = (OpticalElement(opt_model, e))
    # return e_list
    return OpticalElement(opt_model, e)


# --- Graphics elements
class OpticalElement():
    """ mediator class for rendering and editing optical elements """

    def __init__(self, opt_model, e):
        self.opt_model = opt_model
        self.e = e
        self.select_pt = None
        self.move_direction = None
        self.handles = {}
        self.actions = self.edit_shape_actions()
        self.handle_actions = self.e.handle_actions()

    def update_shape(self, view):
        self.e.render_handles(self.opt_model)
        for key, graphics_handle in self.e.handles.items():
            poly = np.array(graphics_handle.polydata)
            poly_gbl, bbox = transform_poly(graphics_handle.tfrm, poly)
            if graphics_handle.polytype == 'polygon':
                p = view.create_polygon(poly_gbl,
                                        fill_color=graphics_handle.color,
                                        zorder=2.5)
            elif graphics_handle.polytype == 'polyline':
                hilite_kwargs = {
                    'color': lo_rgb['profile'],
                    'linewidth': lo_lw['hilite'],
                    'linestyle': '-'
                    }
                priority = 2.

                if key == 'ct':
                    priority = 3.
                    hilite_kwargs['color'] = lo_rgb['ct']
                p = view.create_polyline(poly_gbl,
                                         hilite=hilite_kwargs,
                                         zorder=priority)
            else:
                break
            self.handles[key] = GUIHandle(p, bbox)
        return self.handles

    def render_color(self):
        return self.e.render_color

    def get_label(self):
        if not hasattr(self.e, 'label'):
            self.e.label = 'element'
        return self.e.label

    def edit_shape_actions(self):
        def add_event_data(self, event, handle):
            gbl_pt = np.array([event.xdata, event.ydata])
            lcl_pt = inv_transform_poly(self.e.handles[handle].tfrm, gbl_pt)
            event.lcl_pt = lcl_pt
            if self.select_pt is not None:
                xdata, ydata = self.select_pt[1]
            else:
                xdata, ydata = 0., 0.
            dxdy = event.xdata - xdata, event.ydata - ydata
            event.dxdy = dxdy

        def on_select_shape(fig, handle, event, info):
            handle_actions = self.handle_actions[handle]
            add_event_data(self, event, handle)
            for key, action_obj in handle_actions.items():
                action_obj.actions['press'](fig, event)
            self.select_pt = ((event.x, event.y), (event.xdata, event.ydata))
#            print('select pt:', self.select_pt)
#            print('select pt:', event.x, event.y)

        def on_edit_shape(fig, handle, event, info):
            handle_actions = self.handle_actions[handle]
            x, y = self.select_pt[0]
            xdata, ydata = self.select_pt[1]
            delta_x, delta_y = abs(x - event.x), abs(y - event.y)
            delta_xdata, delta_ydata = (abs(xdata - event.xdata),
                                        abs(ydata - event.ydata))
            if self.move_direction is None:
                if delta_x > delta_y:
                    self.move_direction = 'x'
#                    print('move horizontal: delta x, y:', delta_x, delta_y,
#                          delta_xdata, delta_ydata)
                elif delta_x < delta_y:
                    self.move_direction = 'y'
#                    print('move vertical: delta x, y:', delta_x, delta_y,
#                          delta_xdata, delta_ydata)
                else:
                    self.move_direction = None
#                    print('move same: delta x, y:', delta_x, delta_y)

            add_event_data(self, event, handle)
#            print('move pt:', event.xdata, event.ydata, event.lcl_pt)
            if 'pt' in handle_actions:
                if 'drag' in handle_actions['pt'].actions:
                    handle_actions['pt'].actions['drag'](fig, event,
                                                         event.lcl_pt)
            elif self.move_direction in handle_actions:
                if 'drag' in handle_actions[self.move_direction].actions:
                    if self.move_direction == 'x':
                        handle_actions['x'].actions['drag'](fig, event,
                                                            event.dxdy[0])
                    elif self.move_direction == 'y':
                        handle_actions['y'].actions['drag'](fig, event,
                                                            event.dxdy[1])

            fig.refresh_gui(build='update')

        def on_release_shape(fig, handle, event, info):
            # print('release pt:', event.x, event.y)
            handle_actions = self.handle_actions[handle]
            add_event_data(self, event, handle)
            for key, action_obj in handle_actions.items():
                action_obj.actions['release'](fig, event)
            self.select_pt = None
            self.move_direction = None
            fig.refresh_gui(build='rebuild')

        actions = {}
        actions['press'] = on_select_shape
        actions['drag'] = on_edit_shape
        actions['release'] = on_release_shape
        return actions


class RayBundle():
    """ class for ray bundle from a single field point """

    def __init__(self, opt_model, fld, fld_label, wvl, start_offset,
                 ray_table_callback=None):
        self.opt_model = opt_model
        self.fld = fld
        self.fld_label = fld_label
        self.wvl = wvl
        self.start_offset = start_offset
        self.handles = {}
        self.actions = self.edit_ray_bundle_actions()
        self.handle_actions = {}
        self.ray_table_callback = ray_table_callback

    def get_label(self):
        return self.fld_label

    def render_ray(self, ray, start_seg, tfrms):
        poly = []
        transform_ray_seg(poly, start_seg, tfrms[0])
        for i, r in enumerate(ray[1:], 1):
            transform_ray_seg(poly, r, tfrms[i])
        return np.array(poly)

    def render_shape(self, rayset, start_bundle, tfrms):
        poly1 = []
        transform_ray_seg(poly1, start_bundle[3], tfrms[0])
        for i, r in enumerate(rayset['+Y'].ray[1:], 1):
            transform_ray_seg(poly1, r, tfrms[i])

        poly2 = []
        transform_ray_seg(poly2, start_bundle[4], tfrms[0])
        for i, r in enumerate(rayset['-Y'].ray[1:], 1):
            transform_ray_seg(poly2, r, tfrms[i])

        poly2.reverse()
        poly1.extend(poly2)
        bbox = bbox_from_poly(poly1)
        return poly1, bbox

    def update_shape(self, view):
        wvl = self.opt_model['optical_spec']['wvls'].central_wvl
        self.wvl = wvl
        rayset = trace_boundary_rays_at_field(self.opt_model,
                                              self.fld, wvl,
                                              use_named_tuples=True)

        self.rayset = boundary_ray_dict(self.opt_model, rayset)
        # If the object distance (tfrms[0][1][2]) is greater than the
        #  start_offset, then modify rayset start to match start_offset.
        # Remember object transformation for resetting at the end.
        seq_model = self.opt_model.seq_model
        tfrms = seq_model.gbl_tfrms
        tfrtm0 = tfrms[0]

        start_bundle = [r.ray[0] for r in self.rayset.values()]
        ray_list = [r for r in self.rayset.values()]
        if abs(tfrtm0[1][2]) > self.start_offset:
            rot, t = setup_shift_of_ray_bundle(seq_model, self.start_offset)
            tfrms[0] = (rot, t)
            shift_start_of_ray_bundle(start_bundle, ray_list, rot, t)

        try:
            if view.do_draw_beams:
                poly, bbox = self.render_shape(self.rayset,
                                               start_bundle, tfrms)

                p = view.create_polygon(poly, fill_color=lo_rgb['rayfan_fill'])
                self.handles['shape'] = GUIHandle(p, bbox)

            if view.do_draw_edge_rays:
                cr = self.render_ray(self.rayset['00'].ray,
                                     start_bundle[0], tfrms)
                upr = self.render_ray(self.rayset['+Y'].ray,
                                      start_bundle[3], tfrms)
                lwr = self.render_ray(self.rayset['-Y'].ray,
                                      start_bundle[4], tfrms)
                kwargs = {
                    'linewidth': lo_lw['line'],
                    'color': lo_rgb['ray'],
                    'hilite_linewidth': lo_lw['hilite'],
                    'hilite': lo_rgb['ray'],
                    }
                cr_poly = view.create_polyline(cr, **kwargs)
                self.handles['00'] = GUIHandle(cr_poly, bbox_from_poly(cr))
        
                upr_poly = view.create_polyline(upr, **kwargs)
                self.handles['+Y'] = GUIHandle(upr_poly, bbox_from_poly(upr))
        
                lwr_poly = view.create_polyline(lwr, **kwargs)
                self.handles['-Y'] = GUIHandle(lwr_poly, bbox_from_poly(lwr))

        finally:
            tfrms[0] = tfrtm0

        return self.handles

    def edit_ray_bundle_actions(self):
        def on_select_ray(fig, handle, event, info):
            if handle != 'shape':
                ray_table = self.ray_table_callback()
                ray_table.root = self.rayset[handle].ray
                fig.refresh_gui(build='update')

        actions = {}
        actions['press'] = on_select_ray
        return actions


class RayFanBundle():
    """ class for a RayFan from a single field point """

    def __init__(self, opt_model, ray_fan, start_offset, label='ray fan'):
        self.opt_model = opt_model

        rayerr_filter = ray_fan.rayerr_filter
        ray_fan.rayerr_filter = ('full' if rayerr_filter is None
                                 else rayerr_filter)
        self.ray_fan = ray_fan
        self.start_offset = start_offset
        self.label = label
        self.handles = {}
        self.actions = {}
        self.handle_actions = {}

    def get_label(self):
        return self.label

    def render_ray(self, ray_pkg, start_seg, tfrms):
        poly = []
        ray, op_delta, wvl = ray_pkg
        transform_ray_seg(poly, start_seg, tfrms[0])
        for i, r in enumerate(ray[1:], 1):
            transform_ray_seg(poly, r, tfrms[i])
        return np.array(poly)

    def update_shape(self, view):
        self.ray_fan.update_data()
        fan = self.ray_fan.fan_pkg[0]

        # Remember object transformation for resetting at the end.
        seq_model = self.opt_model.seq_model
        tfrms = seq_model.gbl_tfrms
        tfrtm0 = tfrms[0]

        start_bundle = []
        ray_list = []
        for ray_item in fan:
            ray_pkg = retrieve_ray(ray_item)
            ray = ray_pkg[0]
            start_bundle.append(ray[0])
            ray_list.append(ray_pkg)
            
        # If the object distance (tfrms[0][1][2]) is greater than the
        #  start_offset, then modify rayset start to match start_offset.
        if abs(tfrtm0[1][2]) > self.start_offset:
            rot, t = setup_shift_of_ray_bundle(seq_model, self.start_offset)
            tfrms[0] = (rot, t)
            cr_index = len(ray_list)//2
            shift_start_of_ray_bundle(start_bundle, ray_list, rot, t,
                                      cr_indx=cr_index)

        kwargs = {
            'linewidth': lo_lw['line'],
            'color': lo_rgb['ray'],
            'hilite_linewidth': lo_lw['hilite'],
            'hilite': lo_rgb['ray'],
            }

        for i, rs in enumerate(zip(ray_list, start_bundle)):
            ray_pkg, start_seg = rs
            global_ray = self.render_ray(ray_pkg, start_seg, tfrms)
            ray_poly = view.create_polyline(global_ray, **kwargs)
            self.handles[i] = GUIHandle(ray_poly, bbox_from_poly(global_ray))
            
        tfrms[0] = tfrtm0

        return self.handles


class ParaxialRay():
    """ class for paraxial ray rendering/editing """

    def __init__(self, opt_model, ray, color, seq_start=1, label='paraxial'):
        self.label = label
        self.opt_model = opt_model
        self.ray = ray
        self.seq_start = seq_start
        self.color = color
        self.handles = {}
        self.actions = self.edit_paraxial_layout_actions()
        self.handle_actions = {}
        self.vertex = None

    def get_label(self):
        return self.label

    def render_ray(self, ray, tfrms):
        poly = []
        for i, r in enumerate(ray[self.seq_start:], self.seq_start):
            rot, trns = tfrms[i]
            ps = np.array([0., r[mc.ht], 0.])
            p = rot.dot(ps) + trns
            poly.append([p[2], p[1]])
        return np.array(poly)

    def update_shape(self, view):
        seq_model = self.opt_model.seq_model
        tfrms = seq_model.gbl_tfrms

        ray_poly = self.render_ray(self.ray, tfrms)

        hilite_kwargs = {
            'color': self.color,
            'linewidth': lo_lw['hilite'],
            'linestyle': '-'
            }
        priority = 3
        rp = view.create_polyline(ray_poly,
                                  color=self.color,
                                  linewidth=lo_lw['line'],
                                  hilite=hilite_kwargs,
                                  zorder=priority)
        self.handles['shape'] = GUIHandle(rp, bbox_from_poly(ray_poly))

        return self.handles

    def apply_data(self, vertex, lcl_pt):
        ray = self.ray
        p = ray[vertex-1]
        c = ray[vertex]
        c_slp_new = (lcl_pt[1] - c[mc.ht])/lcl_pt[0]
        pwr = (p[mc.slp] - c_slp_new)/c[mc.ht]
        self.opt_model.seq_model.ifcs[vertex].optical_power = pwr

    def edit_paraxial_layout_actions(self):
        def add_event_data(shape, event, handle, info):
            gbl_pt = np.array([event.xdata, event.ydata])
            lcl_pt = inv_transform_poly(shape.tfrm, gbl_pt)
            event.lcl_pt = lcl_pt

        def on_select_point(fig, handle, event, info):
            self.vertex = self.seq_start
            if 'ind' in info:
                self.vertex += info['ind'][0]
            seq_model = self.opt_model.seq_model
            self.tfrm = seq_model.gbl_tfrms[self.vertex]

        def on_edit_point(fig, handle, event, info):
            add_event_data(self, event, handle, info)
            self.apply_data(self.vertex, event.lcl_pt)
            fig.refresh_gui(build='update')

        def on_release_point(fig, handle, event, info):
            add_event_data(self, event, handle, info)
            self.apply_data(self.vertex, event.lcl_pt)
            fig.refresh_gui(build='rebuild')
            self.vertex = None
            self.tfrm = None

        actions = {}
        actions['press'] = on_select_point
        actions['drag'] = on_edit_point
        actions['release'] = on_release_point
        return actions


# --- Layout manager
class LensLayout():
    """ manager for live layout graphics entities """

    def __init__(self, opt_model, is_dark=True, **kwargs):
        self.opt_model = opt_model
        self.ray_table = None

        light_or_dark(is_dark=is_dark)

        seq_model = self.opt_model['seq_model']
        ele_model = self.opt_model['ele_model']
        part_tree = self.opt_model['part_tree']
        if len(part_tree.nodes_with_tag(tag='#element')) == 0:
            parttree.elements_from_sequence(ele_model, seq_model, part_tree)

    def sync_light_or_dark(self, is_dark):
        global lo_rgb
        lo_rgb = light_or_dark(is_dark)

    def system_length(self, ele_bbox, offset_factor=0.05):
        """ returns system length and ray start offset """
        specsheet = self.opt_model['specsheet']
        fod = self.opt_model['optical_spec'].parax_data[2]
        ele_length = ele_bbox[1][0] - ele_bbox[0][0]
        image_thi = abs(self.opt_model['seq_model'].gaps[-1].thi)
        if specsheet.imager_defined():
            if specsheet.conjugate_type == 'finite':
                return specsheet.imager.tt, (2/3)*specsheet.imager.sp
            elif specsheet.conjugate_type == 'infinite':
                if fod.efl == 0:
                    estimated_length = ele_length
                else:
                    # img_dist = abs(fod.img_dist)
                    # estimated_length = ele_length + img_dist
                    estimated_length = ele_length + image_thi
                return estimated_length, offset_factor*estimated_length
                # return specsheet.imager.sp, offset_factor*specsheet.imager.sp
            elif specsheet.conjugate_type == 'afocal':
                return ele_length, offset_factor*ele_length
        else:
            if fod.efl == 0:
                estimated_length = ele_length
            else:
                estimated_length = ele_length + image_thi
            return estimated_length, offset_factor*estimated_length

    def create_element_entities(self, view):
        opm = self.opt_model
        pt = opm.part_tree
        if opm['ss'].conjugate_type == 'infinite':
            e_nodes = pt.nodes_with_tag(tag='#element#dummyifc#airgap',
                                        not_tag='#object')
        else:
            e_nodes = pt.nodes_with_tag(tag='#element#dummyifc#airgap')
        elements = [create_optical_element(opm, e_node.id)
                    for e_node in e_nodes]
        return elements

    def create_ray_entities(self, view, start_offset):
        ray_bundles = []
        fov = self.opt_model['optical_spec']['fov']
        wvl = self.opt_model['seq_model'].central_wavelength()
        for i, fld in enumerate(fov.fields):
            fld_label = fov.index_labels[i]
            rb = RayBundle(self.opt_model, fld, fld_label, wvl, start_offset,
                           ray_table_callback=self.get_ray_table)
            ray_bundles.append(rb)
        return ray_bundles

    def create_ray_fan_entities(self, view, start_offset, num_rays=21):
        ray_fan_bundles = []
        opt_model = self.opt_model
        fov = opt_model['optical_spec']['fov']
        for i, fld in enumerate(fov.fields):
            fld_label = fov.index_labels[i]
            rayfan = RayFan(opt_model, f=fld, xyfan='y', num_rays=num_rays,
                            label=fld_label)
            rb = RayFanBundle(opt_model, rayfan, start_offset)
            ray_fan_bundles.append(rb)
        return ray_fan_bundles

    def create_paraxial_ray_entities(self, view):
        parax_model = self.opt_model['parax_model']

        ax_poly = ParaxialRay(self.opt_model, parax_model.ax,
                              color=lo_rgb['ax_ray'],
                              label='axial ray')

        pr_poly = ParaxialRay(self.opt_model, parax_model.pr,
                              color=lo_rgb['pr_ray'],
                              label='chief ray')
        return [ax_poly, pr_poly]

    def get_ray_table(self):
        if self.ray_table is None:
            self.ray_table = cmds.create_ray_table_model(self.opt_model, None)
            gui_parent = self.opt_model.app_manager.gui_parent
            gui_parent.create_table_view(self.ray_table, "Ray Table")

        return self.ray_table

    def register_commands(self, *args, **kwargs):
        self.apply_fct = kwargs.pop('apply_fct')
        fig = kwargs.pop('figure')
        actions = self.add_element_cmd_actions(**kwargs)

        def do_command_action(event, target, event_key):
            nonlocal fig, actions
            if target is not None:
                shape, handle = target.artist.shape
                if handle == 'ct' and isinstance(shape.e, ele.AirGap):
                    try:
                        action = actions[event_key]
                        action(fig, shape, event, target.info)
                    except KeyError:
                        pass
        fig.do_action = do_command_action

    def add_element_cmd_actions(self, **kwargs):
        idx = None

        def add_event_data(tfrm, event):
            gbl_pt = np.array([event.xdata, event.ydata])
            lcl_pt = inv_transform_poly(tfrm, gbl_pt)
            event.lcl_pt = lcl_pt

        def on_select_cmd(fig, shape, event, info):
            nonlocal idx
            idx = shape.e.idx
            tfrm = self.opt_model.seq_model.gbl_tfrms[idx]
            add_event_data(tfrm, event)
            self.apply_fct(self.opt_model, idx, event.lcl_pt, **kwargs)
            fig.refresh_gui(build='update')

#        def on_edit_cmd(fig, shape, event, info):
#            add_event_data(self, event)
#            self.apply_fct(idx, event.lcl_pt)
#            fig.refresh_gui()
#        actions['drag'] = on_edit_cmd

        def on_release_cmd(fig, shape, event, info):
            nonlocal idx
#            add_event_data(self, event)
#            self.apply_fct(idx, event.lcl_pt)
#            fig.refresh_gui(build='rebuild')
            idx = None

        actions = {}
        actions['press'] = on_select_cmd
        actions['release'] = on_release_cmd
        return actions


# --- Command functions
def split_gap(opt_model, idx, lcl_pt):
    """ split g=gap[idx] into t_old = t_0 + t_k using t_0 = lcl_pt.x """
    g = opt_model.seq_model.gaps[idx]
    x, y = lcl_pt
    t_k = g.thi - x
    t_0 = x
    return g, t_0, t_k


def add_elements(opt_model, idx, lcl_pt, create, **kwargs):
    g, t_0, t_k = split_gap(opt_model, idx, lcl_pt)
    g.thi = t_0
    opt_model.insert_ifc_gp_ele(*create(**kwargs), idx=idx, t=t_k,
                                insert=True)


def add_reflector(opt_model, idx, lcl_pt, create, **kwargs):
    g, t_0, t_k = split_gap(opt_model, idx, lcl_pt)
    g.thi = t_0
    opt_model.insert_ifc_gp_ele(*create(**kwargs), idx=idx, t=-t_k,
                                insert=True)


def add_thinlens(opt_model, idx, lcl_pt, **kwargs):
    add_elements(opt_model, idx, lcl_pt, ele.create_thinlens, **kwargs)


def add_lens(opt_model, idx, lcl_pt, **kwargs):
    g, t_0, t_k = split_gap(opt_model, idx, lcl_pt)
    t_old = g.thi
    if 'th' in kwargs:
        t_ct = kwargs['th']
    else:
        t_ct = 0.05*t_old
    t0_new = t_0 - t_ct/2
    tk_new = t_k - t_ct/2
    g.thi = t0_new

    descriptor = ele.create_lens(th=t_ct)
    opt_model.insert_ifc_gp_ele(*descriptor, idx=idx, t=tk_new, insert=True)


def add_mirror(opt_model, idx, lcl_pt, **kwargs):
    add_reflector(opt_model, idx, lcl_pt, ele.create_mirror, **kwargs)


def add_conic(opt_model, idx, lcl_pt, **kwargs):
    add_reflector(opt_model, idx, lcl_pt, ele.create_mirror, **kwargs)


def add_doublet(opt_model, idx, lcl_pt, **kwargs):
    add_elements(opt_model, idx, lcl_pt, ele.create_cemented_doublet, **kwargs)


class GlassDropAction():

    def dragEnterEvent(self, view, event):
        def glass_target_filter(artist):
            shape, handle = artist.shape
            if handle == 'shape' and 'shape' in shape.handle_actions:
                if 'glass' in shape.handle_actions['shape']:
                    return False
            return True
        view.figure.artist_filter = glass_target_filter

    def dragMoveEvent(self, view, event):
        x, y = view.mouseEventCoords(event.pos())
        view.motion_notify_event(x, y, guiEvent=event)

    def dragLeaveEvent(self, view, event):
        view.figure.artist_filter = None

    def dropEvent(self, view, event):
        dropped_it = False
        fig = view.figure
        if fig.hilited is not None:
            target = fig.hilited
            shape, handle = target.artist.shape
            if 'glass' in shape.handle_actions['shape']:
                action = shape.handle_actions['shape']['glass']
                action(fig, event)
                dropped_it = True
        view.figure.artist_filter = None
        return dropped_it
