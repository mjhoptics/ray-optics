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

from rayoptics.optical import transform
from rayoptics.optical.trace import (trace_boundary_rays_at_field,
                                     boundary_ray_dict, RaySeg)
from rayoptics.optical.elements import (create_thinlens, create_mirror,
                                        create_lens, AirGap)
import rayoptics.optical.model_constants as mc
from rayoptics.util.rgb2mpl import rgb2mpl
import rayoptics.gui.appcmds as cmds


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


def shift_start_of_ray_bundle(start_bundle, ray_bundle, rot, t):
    """ modify ray_bundle so that rays begin "start_offset" from 1st surface

    Args:
        ray_bundle: list of rays in a bundle, i.e. all for one field.
                    ray_bundle[0] is assumed to be the chief ray
        start_offset: z distance rays should start wrt first surface.
                      positive if to left of first surface
        rot: transformation rotation
        t: transformation translation
    """

    for ri, r in enumerate(ray_bundle):
        b4_pt = rot.dot(r.ray[1].p) - t
        b4_dir = rot.dot(r.ray[0].d)
        if ri == 0:
            # For the chief ray, use the input offset.
            dst = -b4_pt[2]/b4_dir[2]
        else:
            pt0 = start_bundle[0].p
            dir0 = start_bundle[0].d
            # Calculate distance along ray to plane perpendicular to
            #  the chief ray.
            dst = -(b4_pt - pt0).dot(dir0)/b4_dir.dot(dir0)
        pt = b4_pt + dst*b4_dir
        ray0 = RaySeg(pt, b4_dir, dst, r.ray[0].nrml)
        start_bundle[ri] = ray0


def create_optical_element(opt_model, e):
    return OpticalElement(opt_model, e)


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
                                        fill_color=self.render_color(),
                                        zorder=2.5)
            elif graphics_handle.polytype == 'polyline':
                priority = 2.
                hc = 'red'
                if key == 'ct':
                    priority = 3.
                    hc = 'blue'
                p = view.create_polyline(poly_gbl, hilite=hc, zorder=priority)
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

            fig.refresh_gui()

        def on_release_shape(fig, handle, event, info):
            # print('release pt:', event.x, event.y)
            handle_actions = self.handle_actions[handle]
            add_event_data(self, event, handle)
            for key, action_obj in handle_actions.items():
                action_obj.actions['release'](fig, event)
            self.select_pt = None
            self.move_direction = None
            fig.refresh_gui()

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
        rndr_clr = [254, 197, 254, 64]  # magenta, 25%
        wvl = self.opt_model.optical_spec.spectral_region.central_wvl
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
        if abs(tfrtm0[1][2]) > self.start_offset:
            rot, t = setup_shift_of_ray_bundle(seq_model, self.start_offset)
            tfrms[0] = (rot, t)
            shift_start_of_ray_bundle(start_bundle, self.rayset.values(),
                                      rot, t)

        try:
            poly, bbox = self.render_shape(self.rayset,
                                           start_bundle, tfrms)
            cr = self.render_ray(self.rayset['00'].ray,
                                 start_bundle[0], tfrms)
            upr = self.render_ray(self.rayset['+Y'].ray,
                                  start_bundle[3], tfrms)
            lwr = self.render_ray(self.rayset['-Y'].ray,
                                  start_bundle[4], tfrms)
        finally:
            tfrms[0] = tfrtm0

        p = view.create_polygon(poly, fill_color=rndr_clr)
        self.handles['shape'] = GUIHandle(p, bbox)

        cr_poly = view.create_polyline(cr)
        self.handles['00'] = GUIHandle(cr_poly, bbox_from_poly(cr))

        upr_poly = view.create_polyline(upr)
        self.handles['+Y'] = GUIHandle(upr_poly, bbox_from_poly(upr))

        lwr_poly = view.create_polyline(lwr)
        self.handles['-Y'] = GUIHandle(lwr_poly, bbox_from_poly(lwr))

        return self.handles

    def edit_ray_bundle_actions(self):
        def on_select_ray(fig, handle, event, info):
            if handle != 'shape':
                ray_table = self.ray_table_callback()
                ray_table.root = self.rayset[handle].ray
                fig.refresh_gui()

        actions = {}
        actions['press'] = on_select_ray
        return actions


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

        priority = 3
        rp = view.create_polyline(ray_poly, color=self.color, linewidth=1,
                                  hilite=self.color, zorder=priority)
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
            fig.refresh_gui()

        def on_release_point(fig, handle, event, info):
            add_event_data(self, event, handle, info)
            self.apply_data(self.vertex, event.lcl_pt)
            fig.refresh_gui()
            self.vertex = None
            self.tfrm = None

        actions = {}
        actions['press'] = on_select_point
        actions['drag'] = on_edit_point
        actions['release'] = on_release_point
        return actions


class LensLayout():
    """ manager for live layout graphics entities """

    def __init__(self, opt_model, **kwargs):
        self.opt_model = opt_model
        self.ray_table = None

        ele_model = self.opt_model.ele_model
        if len(ele_model.elements) == 0:
            ele_model.elements_from_sequence(self.opt_model.seq_model)

    def system_length(self, ele_bbox, offset_factor=0.05):
        """ returns system length and ray start offset """
        specsheet = self.opt_model.specsheet
        osp = self.opt_model.optical_spec
        ele_length = ele_bbox[1][0] - ele_bbox[0][0]
        if specsheet.imager_defined():
            if specsheet.conjugate_type == 'finite':
                return specsheet.imager.tt, (2/3)*specsheet.imager.sp
            elif specsheet.conjugate_type == 'infinite':
                img_dist = abs(osp.parax_data[2].img_dist)
                return ele_length+img_dist, offset_factor*img_dist
                # return specsheet.imager.sp, offset_factor*specsheet.imager.sp
        else:
            img_dist = abs(osp.parax_data[2].img_dist)
            return ele_length+img_dist, offset_factor*img_dist

    def create_element_model(self, view):
        ele_model = self.opt_model.ele_model
        elements = [create_optical_element(self.opt_model, e)
                    for e in ele_model.elements]
        return elements

    def create_ray_model(self, view, start_offset):
        ray_bundles = []
        fov = self.opt_model.optical_spec.field_of_view
        wvl = self.opt_model.seq_model.central_wavelength()
        for i, fld in enumerate(fov.fields):
            fld_label = fov.index_labels[i]
            rb = RayBundle(self.opt_model, fld, fld_label, wvl, start_offset,
                           ray_table_callback=self.get_ray_table)
            ray_bundles.append(rb)
        return ray_bundles

    def create_paraxial_layout(self, view):
        parax_model = self.opt_model.parax_model

        ax_color = rgb2mpl([138, 43, 226])  # blueviolet
        ax_poly = ParaxialRay(self.opt_model, parax_model.ax, color=ax_color,
                              label='axial ray')
        pr_color = rgb2mpl([198, 113, 113])  # sgi salmon
        pr_poly = ParaxialRay(self.opt_model, parax_model.pr, color=pr_color,
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
                if handle == 'ct' and isinstance(shape.e, AirGap):
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
            fig.refresh_gui()

#        def on_edit_cmd(fig, shape, event, info):
#            add_event_data(self, event)
#            self.apply_fct(idx, event.lcl_pt)
#            fig.refresh_gui()
#        actions['drag'] = on_edit_cmd

        def on_release_cmd(fig, shape, event, info):
            nonlocal idx
#            add_event_data(self, event)
#            self.apply_fct(idx, event.lcl_pt)
#            fig.refresh_gui()
            idx = None

        actions = {}
        actions['press'] = on_select_cmd
        actions['release'] = on_release_cmd
        return actions


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
    add_elements(opt_model, idx, lcl_pt, create_thinlens, **kwargs)


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

    seq, ele = create_lens(th=t_ct)
    opt_model.insert_ifc_gp_ele(seq, ele, idx=idx, t=tk_new, insert=True)


def add_mirror(opt_model, idx, lcl_pt, **kwargs):
    add_reflector(opt_model, idx, lcl_pt, create_mirror, **kwargs)


def add_conic(opt_model, idx, lcl_pt, **kwargs):
    add_reflector(opt_model, idx, lcl_pt, create_mirror, **kwargs)


def add_doublet(opt_model, idx, lcl_pt, **kwargs):
    add_elements(opt_model, idx, lcl_pt, create_lens, **kwargs)
