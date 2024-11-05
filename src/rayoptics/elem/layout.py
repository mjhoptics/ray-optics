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

from rayoptics.gui.util import (GUIHandle, transform_ray_seg, bbox_from_poly,
                                transform_poly, inv_transform_poly)

from rayoptics.raytr.analyses import RayFan
from rayoptics.raytr.trace import (trace_boundary_rays_at_field,
                                   boundary_ray_dict)
from rayoptics.elem import elements as ele

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

    def listobj_str(self):
        o_str = f"{type(self).__name__}: {self.get_label()}\n"
        return o_str

    def update_shape(self, view):
        self.e.render_handles(self.opt_model)
        for key, graphics_handle in self.e.handles.items():
            poly_data = graphics_handle.polydata 
            # print(self.listobj_str()
            #       +f"{type(poly_data).__name__}, # polys={len(poly_data)}, {key} {graphics_handle.polytype}")
            if isinstance(poly_data, tuple):
                polys = []
                for poly_list in poly_data:
                    for poly_seg in poly_list:
                        poly = np.array(poly_seg)
                        poly_gbl = transform_poly(graphics_handle.tfrm, poly)
                        bbox = bbox_from_poly(poly_gbl)
                        polys.append(np.array(poly_gbl))
                poly_gbl = tuple(polys)
            else:
                poly = np.array(poly_data)
                poly_gbl = transform_poly(graphics_handle.tfrm, poly)
                bbox = bbox_from_poly(poly_gbl)

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

            elif graphics_handle.polytype == 'vertex':
                hilite_kwargs = {
                    'color': lo_rgb['edge'],
                    }
                priority = 2.
                p = view.create_vertex(poly_gbl,
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
            fig.refresh_gui(build='update')

        actions = {}
        actions['press'] = on_select_shape
        actions['drag'] = on_edit_shape
        actions['release'] = on_release_shape
        return actions


class RayBundle():
    """ class for ray bundle from a single field point """

    def __init__(self, opt_model, fld, fld_label, wvl, start_offset,
                 ray_table_callback=None, **kwargs):
        self.opt_model = opt_model
        self.fld = fld
        self.fld_label = fld_label
        self.wvl = wvl
        self.start_offset = start_offset
        self.handles = {}
        self.actions = self.edit_ray_bundle_actions()
        self.handle_actions = {}
        self.ray_table_callback = ray_table_callback

    def listobj_str(self):
        o_str = f"RayBundle: {self.fld_label}\n"
        o_str += f"wavelength={self.wvl:10.4f} nm\n"
        return o_str

    def get_label(self):
        return self.fld_label

    def render_ray(self, ray, tfrms):
        poly = []
        for i, r in enumerate(ray):
            transform_ray_seg(poly, r, tfrms[i])
        return np.array(poly)

    def render_shape(self, rayset, tfrms):
        poly1 = []
        for i, r in enumerate(rayset['+Y'].ray):
            transform_ray_seg(poly1, r, tfrms[i])

        poly2 = []
        for i, r in enumerate(rayset['-Y'].ray):
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
                                              use_named_tuples=True,
                                              rayerr_filter='full',
                                              check_apertures=view.clip_rays)

        self.rayset = boundary_ray_dict(self.opt_model, rayset)

        seq_model = self.opt_model.seq_model
        tfrms = seq_model.gbl_tfrms        #     shift_start_of_ray_bundle(start_bundle, ray_list, rot, t)

        if view.do_draw_beams:
            poly, bbox = self.render_shape(self.rayset, tfrms)

            p = view.create_polygon(poly, fill_color=lo_rgb['rayfan_fill'])
            self.handles['shape'] = GUIHandle(p, bbox)

        if view.do_draw_edge_rays:
            cr = self.render_ray(self.rayset['00'].ray, tfrms)
            upr = self.render_ray(self.rayset['+Y'].ray, tfrms)
            lwr = self.render_ray(self.rayset['-Y'].ray, tfrms)
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

        rayerr_filter = ray_fan.rt_kwargs['rayerr_filter']
        ray_fan.rayerr_filter = ('full' if rayerr_filter is None
                                 else rayerr_filter)
        self.ray_fan = ray_fan
        self.start_offset = start_offset
        self.label = label
        self.handles = {}
        self.actions = {}
        self.handle_actions = {}

    def listobj_str(self):
        o_str = f"{type(self).__name__}: {self.get_label()}\n"
        return o_str

    def get_label(self):
        return self.label

    def render_ray(self, ray_pkg, tfrms):
        poly = []
        ray, op_delta, wvl = ray_pkg
        for i, r in enumerate(ray):
            transform_ray_seg(poly, r, tfrms[i])
        return np.array(poly)

    def update_shape(self, view):
        ray_fan = self.ray_fan
        ray_fan.update_data()
        fan = ray_fan.fan_pkg[0]

        seq_model = self.opt_model.seq_model
        tfrms = seq_model.gbl_tfrms

        ray_list = []
        for ray_item in fan:
            ray_pkg = ray_item[2]
            ray_list.append(ray_pkg)

        ray_color = lo_rgb['ray'] if ray_fan.color is None else ray_fan.color
        kwargs = {
            'linewidth': lo_lw['line'],
            'color': ray_color,
            'hilite_linewidth': lo_lw['hilite'],
            'hilite': lo_rgb['ray'],
            }

        for i, ray_pkg in enumerate(ray_list):
            global_ray = self.render_ray(ray_pkg, tfrms)
            ray_poly = view.create_polyline(global_ray, **kwargs)
            self.handles[i] = GUIHandle(ray_poly, bbox_from_poly(global_ray))

        return self.handles


class SingleRay():
    """ class for a Ray from a single field point """

    def __init__(self, opt_model, ray, start_offset, label='single ray'):
        self.opt_model = opt_model

        rayerr_filter = ray.rt_kwargs['rayerr_filter']
        ray.rayerr_filter = ('full' if rayerr_filter is None
                                 else rayerr_filter)
        self.ray = ray
        self.start_offset = start_offset
        self.label = label
        self.handles = {}
        self.actions = {}
        self.handle_actions = {}

    def listobj_str(self):
        o_str = f"{type(self).__name__}: {self.get_label()}\n"
        return o_str

    def get_label(self):
        return self.label

    def render_ray(self, ray_pkg, tfrms):
        poly = []
        ray, op_delta, wvl = ray_pkg
        for i, r in enumerate(ray):
            transform_ray_seg(poly, r, tfrms[i])
        return np.array(poly)

    def update_shape(self, view):
        ray = self.ray
        ray.update_data()

        seq_model = self.opt_model['seq_model']
        tfrms = seq_model.gbl_tfrms
        global_ray = self.render_ray(ray.ray_pkg, tfrms)

        ray_color = lo_rgb['ray'] if ray.color is None else ray.color
        kwargs = {
            'linewidth': lo_lw['line'],
            'color': ray_color,
            'hilite_linewidth': lo_lw['hilite'],
            'hilite': lo_rgb['ray'],
            }

        ray_poly = view.create_polyline(global_ray, **kwargs)
        self.handles['shape'] = GUIHandle(ray_poly, bbox_from_poly(global_ray))

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

    def listobj_str(self):
        o_str = f"{type(self).__name__}: {self.get_label()}\n"
        return o_str

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
            fig.refresh_gui(build='update')
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
        from rayoptics.elem import parttree
        self.opt_model = opt_model
        self.ray_table = None

        light_or_dark(is_dark=is_dark)

        seq_model = self.opt_model['seq_model']
        ele_model = self.opt_model['ele_model']
        part_tree = self.opt_model['part_tree']
        if len(part_tree.nodes_with_tag(tag='#element')) == 0:
            parttree.sequence_to_elements(seq_model, ele_model, part_tree)

    def sync_light_or_dark(self, is_dark):
        global lo_rgb
        lo_rgb = light_or_dark(is_dark)

    def system_length(self, ele_bbox, offset_factor=0.05):
        """ returns system length and ray start offset """
        sm = self.opt_model['seq_model']
        osp = self.opt_model['optical_spec']
        fod = self.opt_model['ar']['parax_data'].fod
        ele_length = ele_bbox[1][0] - ele_bbox[0][0]
        image_thi = abs(self.opt_model['seq_model'].gaps[-1].thi)
        obj_conj = osp.conjugate_type('object')
        if obj_conj == 'finite':
            estimated_length = sm.total_track()
            return estimated_length, offset_factor*sm.overall_length()
        elif obj_conj == 'infinite':
            if fod.efl == 0:
                estimated_length = ele_length
            else:
                estimated_length = ele_length + image_thi
            return estimated_length, offset_factor*estimated_length
        elif osp.is_afocal():
            return ele_length, offset_factor*ele_length

    def renderable_pt_nodes(self, part_filter=''):
        opm = self.opt_model
        osp = opm['optical_spec']
        pt = opm['part_tree']
        
        not_tags = ''
        if osp.conjugate_type('object') == 'infinite':
            not_tags = '#object'
        elif osp.conjugate_type('image') == 'infinite':
            not_tags = '#image'
        elif (osp.conjugate_type('object') == 'infinite' and
              osp.conjugate_type('image') == 'infinite'):
            not_tags = '#object#image'

        e_nodes = pt.nodes_with_tag(tag=part_filter,
                                    not_tag=not_tags)
        return e_nodes

    def create_element_entities(self, view, part_filter=''):
        e_nodes = self.renderable_pt_nodes(part_filter=part_filter)
        elements = [create_optical_element(self.opt_model, e_node.id)
                    for e_node in e_nodes]
        return elements

    def create_oe(self, e):
        """ opaque wrapper for create_optical_element() """
        return create_optical_element(self.opt_model, e)

    def create_ray_entities(self, view, start_offset, **kwargs):
        ray_bundles = []
        fov = self.opt_model['optical_spec']['fov']
        wvl = self.opt_model['seq_model'].central_wavelength()
        for i, fld in enumerate(fov.fields):
            fld_label = fov.index_labels[i]
            rb = RayBundle(self.opt_model, fld, fld_label, wvl, start_offset,
                           ray_table_callback=self.get_ray_table, **kwargs)
            ray_bundles.append(rb)
        return ray_bundles

    def create_ray_fan_entities(self, view, start_offset, 
                                num_rays=21, **kwargs):
        ray_fan_bundles = []
        opt_model = self.opt_model
        fov = opt_model['optical_spec']['fov']
        for i, fld in enumerate(fov.fields):
            fld_label = fov.index_labels[i]
            rayfan = RayFan(opt_model, f=fld, xyfan='y', num_rays=num_rays,
                            label=fld_label, **kwargs)
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

    def dragEnterEvent(self, canvas_fig, event):
        def glass_target_filter(artist):
            shape, handle = artist.shape
            if handle[:5] == 'shape' and 'shape' in shape.handle_actions:
                if 'glass' in shape.handle_actions['shape']:
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
                if 'glass' in shape.handle_actions['shape']:
                    action = shape.handle_actions['shape']['glass']
                    action(fig, event)
                    dropped_it = True
        canvas_fig.figure.artist_filter = None
        return dropped_it
