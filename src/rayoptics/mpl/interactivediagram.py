#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
"""

.. Created on Thu Oct 10 22:02:44 2019

.. codeauthor: Michael J. Hayford
"""

from rayoptics.parax.diagram import Diagram
from rayoptics.mpl.interactivefigure import InteractiveFigure, SelectInfo


class InteractiveDiagram(InteractiveFigure):
    """ Editable version of |ybar| and |nubar| diagrams

    Attributes:
        opt_model: parent optical model
        refresh_gui: function to be called on refresh_gui event
        dgm_type: diagram type, 'ht' or 'slp'
        do_barrel_constraint, bool: display the barrel diamond if True
        barrel_constraint, float: the radius for the barrel constraint
        enable_slide: Display the "bend" or "gap" constaint lines
        bend_or_gap: "bend" | "gap"
        parax_model: if None, 'ifcs' else parax_mode for layer parax_model_key
        parax_model_key: "ifcs" | "eles" | "asm" | "sys"
    """

    def __init__(self, opt_model, dgm_type, refresh_gui=None,
                 do_barrel_constraint=False, barrel_constraint=1.0,
                 enable_slide=False, bend_or_gap='bend',
                 parax_model=None, parax_model_key='ifcs', **kwargs):
        self.refresh_gui = refresh_gui
        if parax_model is None:
            self.parax_model = opt_model.parax_model
            self.parax_model_key = 'ifcs'
        else:
            self.parax_model = parax_model
            self.parax_model_key = kwargs.get('parax_model_key', 'root')

        is_dark = kwargs['is_dark'] if 'is_dark' in kwargs else False
        self.diagram = Diagram(
            opt_model, self.parax_model, self.parax_model_key, dgm_type,
            do_barrel_constraint=do_barrel_constraint,
            barrel_constraint=barrel_constraint,
            bend_or_gap=bend_or_gap, is_dark=is_dark
            )
        self.setup_dgm_type(dgm_type)
        self.enable_slide = enable_slide

        super().__init__(**kwargs)

    def setup_dgm_type(self, dgm_type):
        if dgm_type == 'ht':
            self.x_label = r'$\overline{y}$'
            self.y_label = 'y'
            self.header = r'$y-\overline{y}$ Diagram'
        elif dgm_type == 'slp':
            self.x_label = r'$\overline{\omega}$'
            self.y_label = r'$\omega$'
            self.header = r'$\omega-\overline{\omega}$ Diagram'

    def sync_light_or_dark(self, is_dark, **kwargs):
        self.diagram.sync_light_or_dark(is_dark)
        super().sync_light_or_dark(is_dark, **kwargs)

    def update_data(self, **kwargs):
        self.artists = []
        self.sys_bbox = self.diagram.update_data(self, **kwargs)
        return self

    def get_artist_for_handle(self, handle) -> SelectInfo:
        def persistent_id(shape):
            try:
                pid = shape.listobj_str()
            except AttributeError as ae:
                # print(f"persistent_id exception: {shape}")
                pid = shape[0].listobj_str()
            return pid

        (shape, hdl), info = handle
        shape_id = persistent_id(shape)
        for a in self.artists:
            if shape_id == persistent_id(a.shape):
                return SelectInfo(a, info)
        return None
    
    def action_complete(self):
        super().action_complete()
        args = tuple()
        kwargs = {'figure': self,
                  }
        self.diagram.register_commands(*args, **kwargs)

    def fit_axis_limits(self):
        return self.diagram.fit_axis_limits()
