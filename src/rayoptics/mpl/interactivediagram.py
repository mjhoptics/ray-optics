#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
"""

.. Created on Thu Oct 10 22:02:44 2019

.. codeauthor: Michael J. Hayford
"""

from rayoptics.gui.diagram import Diagram
from rayoptics.mpl.interactivefigure import InteractiveFigure


class InteractiveDiagram(InteractiveFigure):
    """ Editable version of optical system layout, aka Live Layout

    Attributes:
        opt_model: parent optical model
        refresh_gui: function to be called on refresh_gui event
        dgm_type: diagram type, 'ht' or 'slp'
    """

    def __init__(self, opt_model, dgm_type, refresh_gui=None,
                 do_barrel_constraint=False, barrel_constraint=1.0,
                 enable_slide=False, bend_or_gap='bend', **kwargs):
        self.refresh_gui = refresh_gui
        self.parax_model = opt_model.parax_model
        is_dark = kwargs['is_dark'] if 'is_dark' in kwargs else False
        self.diagram = Diagram(opt_model, dgm_type,
                               do_barrel_constraint=do_barrel_constraint,
                               barrel_constraint=barrel_constraint,
                               bend_or_gap=bend_or_gap, is_dark=is_dark)
        self.setup_dgm_type(dgm_type)
        self.enable_slide = enable_slide

        self.build = 'rebuild'

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
        self.build = 'rebuild'
        return self

    def action_complete(self):
        super().action_complete()
        self.diagram.register_commands((), figure=self)

    def fit_axis_limits(self):
        return self.diagram.fit_axis_limits()
