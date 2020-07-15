#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
""" Interactive layout figure with paraxial editing

.. Created on Thu Mar 14 10:20:33 2019

.. codeauthor: Michael J. Hayford
"""

import numpy as np

from rayoptics.gui.util import bbox_from_poly

from rayoptics.mpl.interactivefigure import InteractiveFigure

from rayoptics.elem.layout import LensLayout


class InteractiveLayout(InteractiveFigure):
    """ Editable version of optical system layout, aka Live Layout

    Attributes:
        opt_model: parent optical model
        refresh_gui: function to be called on refresh_gui event
        offset_factor: how much to draw rays before first surface
        do_draw_rays: if True, draw edge rays
        do_paraxial_layout: if True, draw editable paraxial axial and chief ray
    """

    def __init__(self, opt_model, refresh_gui=None,
                 offset_factor=0.05,
                 do_draw_rays=True,
                 do_paraxial_layout=False,
                 **kwargs):
        self.refresh_gui = refresh_gui
        is_dark = kwargs['is_dark'] if 'is_dark' in kwargs else False
        self.layout = LensLayout(opt_model, is_dark=is_dark)
        self.do_draw_rays = do_draw_rays
        self.do_paraxial_layout = do_paraxial_layout
        self.offset_factor = offset_factor

        super().__init__(**kwargs)

    def sync_light_or_dark(self, is_dark, **kwargs):
        self.layout.sync_light_or_dark(is_dark)
        super().sync_light_or_dark(is_dark, **kwargs)

    def update_data(self, **kwargs):
        self.artists = []
        concat_bbox = []
        layout = self.layout

        self.ele_shapes = layout.create_element_model(self)
        self.ele_bbox = self.update_patches(self.ele_shapes)
        concat_bbox.append(self.ele_bbox)

        if self.do_draw_rays:
            sl_so = layout.system_length(self.ele_bbox,
                                         offset_factor=self.offset_factor)
            system_length, start_offset = sl_so
            self.ray_shapes = layout.create_ray_model(self, start_offset)
            self.ray_bbox = self.update_patches(self.ray_shapes)
            concat_bbox.append(self.ray_bbox)

        if self.do_paraxial_layout:
            self.parax_shapes = layout.create_paraxial_layout(self)
            self.parax_bbox = self.update_patches(self.parax_shapes)
            concat_bbox.append(self.parax_bbox)

        sys_bbox = np.concatenate(concat_bbox)
        self.sys_bbox = bbox_from_poly(sys_bbox)

        return self

    def action_complete(self):
        super().action_complete()
        self.do_action = self.do_shape_action

    def fit_axis_limits(self):
        return self.sys_bbox
