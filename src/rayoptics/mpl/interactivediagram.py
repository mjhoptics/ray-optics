#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
"""

.. Created on Thu Oct 10 22:02:44 2019

.. codeauthor: Michael J. Hayford
"""

import numpy as np

from rayoptics.gui.diagram import Diagram, DiagramNode, DiagramEdge
from rayoptics.gui.util import bbox_from_poly

from rayoptics.mpl.interactivefigure import InteractiveFigure

from rayoptics.optical.elements import (create_thinlens, create_mirror,
                                        create_lens)


def create_parax_design_commands(fig):
    cmds = []
    dgm = fig.diagram
    # initialize dgm with a Select command
    dgm.register_commands((), figure=fig)
    # Select an existing point
    cmds.append(('Select', (dgm.register_commands, (), {})))
    # Add thin lens
    cmds.append(('Add Thin Lens',
                 (dgm.register_commands, (),
                  {'node_init': create_thinlens,
                   'factory': create_thinlens,
                   'interact_mode': 'transmit'})))
    # Add lens
    cmds.append(('Add Lens', (dgm.register_commands, (),
                              {'node_init': create_thinlens,
                               'factory': create_lens,
                               'interact_mode': 'transmit'})))
    # Add mirror
    cmds.append(('Add Mirror',
                 (dgm.register_commands, (),
                  {'node_init': create_mirror,
                   'factory': create_mirror,
                   'interact_mode': 'reflect'})))

    return cmds


class InteractiveDiagram(InteractiveFigure):
    """ Editable version of optical system layout, aka Live Layout

    Attributes:
        opt_model: parent optical model
        refresh_gui: function to be called on refresh_gui event
        dgm_type: diagram type, 'ht' or 'slp'
    """
    def __init__(self, opt_model, refresh_gui, dgm_type,
                 **kwargs):
        self.refresh_gui = refresh_gui
        self.parax_model = opt_model.parax_model
        self.diagram = Diagram(opt_model, dgm_type)
        self.setup_dgm_type(dgm_type)

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

    def update_data(self):
        self.artists = []
        self.sys_bbox = self.diagram.update_data(self)
        self.build == 'full_rebuild'
        return self

    def action_complete(self):
        self.diagram.register_commands((), figure=self)

    def update_axis_limits(self):
        x_min, x_max = self.fit_data_range([x[0] for x in self.diagram.shape])
        y_min, y_max = self.fit_data_range([x[1] for x in self.diagram.shape])
        self.ax.set_xlim(x_min, x_max)
        self.ax.set_ylim(y_min, y_max)

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
