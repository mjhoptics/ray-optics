#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
"""

.. Created on Thu Oct 10 22:02:44 2019

.. codeauthor: Michael J. Hayford
"""

from rayoptics.gui.layout import bbox_from_poly

from rayoptics.mpl.interactivefigure import InteractiveFigure

from rayoptics.optical.paraxialdesign import Diagram
from rayoptics.optical.model_constants import ht, slp
from rayoptics.optical.elements import (create_thinlens, create_mirror,
                                        create_lens)


def create_parax_design_commands(fig):
    cmds = []
    dgm = fig.diagram
    # Select an existing point
    cmds.append(('Select', (dgm.register_commands, (),
                            {'cmd_actions': dgm.edit_diagram_actions})))
    # Add thin lens
    cmds.append(('Add Thin Lens',
                 (dgm.register_commands, (),
                  {'cmd_actions': dgm.add_element_actions,
                   'node_init': create_thinlens,
                   'factory': create_thinlens})))
    # Add lens
    cmds.append(('Add Lens', (dgm.register_commands, (),
                              {'cmd_actions': dgm.add_element_actions,
                               'node_init': create_thinlens,
                               'factory': create_lens})))
    # Add mirror
    cmds.append(('Add Mirror',
                 (dgm.register_commands, (),
                  {'cmd_actions': dgm.add_element_actions,
                   'node_init': create_mirror,
                   'factory': create_mirror})))

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

        self.skip_build = False

        super().__init__(**kwargs)

    def setup_dgm_type(self, dgm_type):
        if dgm_type == 'ht':
            self.type_sel = ht
            self.x_label = r'$\overline{y}$'
            self.y_label = 'y'
            self.apply_data = self.parax_model.apply_ht_dgm_data
            self.header = r'$y-\overline{y}$ Diagram'
        elif dgm_type == 'slp':
            self.type_sel = slp
            self.x_label = r'$\overline{\omega}$'
            self.y_label = r'$\omega$'
            self.apply_data = self.parax_model.apply_slope_dgm_data
            self.header = r'$\omega-\overline{\omega}$ Diagram'

    def update_data(self):
        self.artists = []

        if not self.skip_build:
            self.parax_model.build_lens()
        self.skip_build = False

        dgm_bbox = self.update_patches([self.diagram])

        self.sys_bbox = bbox_from_poly(dgm_bbox)

        return self
