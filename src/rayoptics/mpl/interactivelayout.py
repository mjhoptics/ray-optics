#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
""" Interactive layout figure with paraxial editing

.. Created on Thu Mar 14 10:20:33 2019

.. codeauthor: Michael J. Hayford
"""

import numpy as np

import rayoptics
from rayoptics.mpl.interactivefigure import InteractiveFigure
from rayoptics.gui.util import bbox_from_poly, scale_bounds


class InteractiveLayout(InteractiveFigure):
    """ Editable version of optical system layout, aka Live Layout

    Attributes:
        opt_model: parent optical model
        refresh_gui: function to be called on refresh_gui event
        offset_factor: how much to draw rays before first surface
        do_draw_rays: if True, draw edge rays
        do_paraxial_layout: if True, draw editable paraxial axial and chief ray
        entity_factory_list: list of drawable entity factories. Allows new
                             drawables to be added to the layout.
    """

    def __init__(self, opt_model, refresh_gui=None,
                 offset_factor=0.05,
                 do_draw_parts=True,
                 part_filter='#element#dummyifc#space#airgap',
                 do_draw_rays=True,
                 do_draw_beams=True,
                 do_draw_edge_rays=True,
                 do_draw_ray_fans=False,
                 num_rays_in_fan=11,
                 clip_rays=False,
                 do_paraxial_layout=False,
                 entity_factory_list=None,
                 **kwargs):
        self.refresh_gui = refresh_gui
        is_dark = kwargs['is_dark'] if 'is_dark' in kwargs else False
        self.layout = rayoptics.elem.layout.LensLayout(opt_model, 
                                                       is_dark=is_dark)

        self.do_draw_parts = do_draw_parts
        self.part_filter = part_filter

        if do_draw_rays:
            self.do_draw_beams = do_draw_beams
            self.do_draw_edge_rays = do_draw_edge_rays
        else:
            self.do_draw_beams = False
            self.do_draw_edge_rays = False
            
        self.do_draw_ray_fans = do_draw_ray_fans
        self.num_rays_in_fan = num_rays_in_fan
        self.do_paraxial_layout = do_paraxial_layout
        self.clip_rays = clip_rays
        self.offset_factor = offset_factor
        
        if entity_factory_list is None:
            self.entity_factory_list = []
        else:
            self.entity_factory_list = entity_factory_list

        super().__init__(**kwargs)

    def sync_light_or_dark(self, is_dark, **kwargs):
        self.layout.sync_light_or_dark(is_dark)
        super().sync_light_or_dark(is_dark, **kwargs)

    def update_data(self, **kwargs):
        self.artists = []
        concat_bbox = []
        layout = self.layout
        build = kwargs.get('build', 'rebuild')

        if self.do_draw_parts:
            if build == 'rebuild':
                self.ele_shapes = []
            e_nodes = layout.renderable_pt_nodes(part_filter=self.part_filter)
            e_node_set = {e.id for e in e_nodes}
            ele_shapes_set = {ele.e for ele in self.ele_shapes}
            to_remove = list(ele_shapes_set - e_node_set)
            self.ele_shapes = [ele for ele in self.ele_shapes 
                               if ele.e not in to_remove]

            to_add = list(e_node_set - ele_shapes_set)
            for e in to_add:
                ele = layout.create_oe(e)
                self.ele_shapes.append(ele)
            self.ele_bbox = self.update_patches(self.ele_shapes)
            concat_bbox.append(self.ele_bbox)

        # if self.do_draw_beams 
        #  or self.do_draw_edge_rays 
        #  or self.do_draw_ray_fans:
        self.sl_so = layout.system_length(self.ele_bbox, 
                                          offset_factor=self.offset_factor)
        system_length, start_offset = self.sl_so

        if self.do_draw_beams or self.do_draw_edge_rays:
            if build == 'rebuild':
                self.ray_shapes = layout.create_ray_entities(
                    self, start_offset, clip_rays=self.clip_rays)
            self.ray_bbox = self.update_patches(self.ray_shapes)

        if self.do_draw_ray_fans:
            if build == 'rebuild':
                self.rayfan_shapes = layout.create_ray_fan_entities(
                    self, start_offset, num_rays=self.num_rays_in_fan, 
                    clip_rays=self.clip_rays)
            self.rayfan_bbox = self.update_patches(self.rayfan_shapes)

        if self.do_paraxial_layout:
            if build == 'rebuild':
                self.parax_shapes = layout.create_paraxial_ray_entities(self)
            self.parax_bbox = self.update_patches(self.parax_shapes)
        
        for ef in self.entity_factory_list:
            # each entity factory is a tuple of callable, args, and kwargs
            # the callable signature includes the current figure.
            if build == 'rebuild':
                ef_fct, ef_args, ef_kwargs = ef
                self.ef_shapes = ef_fct(self, *ef_args, **ef_kwargs)
            ef_bbox = self.update_patches(self.ef_shapes)

        sys_bbox = np.concatenate(concat_bbox)
        self.sys_bbox = scale_bounds(bbox_from_poly(sys_bbox),
                                     self.offset_factor)

        return self

    def action_complete(self):
        super().action_complete()
        self.do_action = self.do_shape_action

    def fit_axis_limits(self):
        return self.sys_bbox
