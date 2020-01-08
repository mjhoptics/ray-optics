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
                                        create_lens, create_from_file)


def create_parax_design_commands(fig):
    cmds = []
    dgm = fig.diagram
    # initialize dgm with a Select command
    dgm.register_commands((), figure=fig)
    # Select an existing point
    cmds.append(('Select', (dgm.register_commands, (), {})))
    # Add thin lens
    cmds.append(('Add Thin Lens',
                 (dgm.register_add_replace_element, (),
                  {'node_init': create_thinlens,
                   'factory': create_thinlens,
                   'interact_mode': 'transmit'})))
    # Add lens
    cmds.append(('Add Lens', (dgm.register_add_replace_element, (),
                              {'node_init': create_thinlens,
                               'factory': create_lens,
                               'interact_mode': 'transmit'})))
    # Add mirror
    cmds.append(('Add Mirror',
                 (dgm.register_add_replace_element, (),
                  {'node_init': create_mirror,
                   'factory': create_mirror,
                   'interact_mode': 'reflect'})))
    # replace with file
    filename = '/Users/Mike/Developer/PyProjects/ray-optics/models/Sasian Triplet.roa'
    def cff(**kwargs):
        return create_from_file(filename, **kwargs)

    cmds.append(('Sasian Triplet',
                 (dgm.register_add_replace_element, (),
                  {'filename': filename,
                   'node_init': create_thinlens,
                   'factory': cff,
                   'interact_mode': 'transmit'})))

    return cmds


class InteractiveDiagram(InteractiveFigure):
    """ Editable version of optical system layout, aka Live Layout

    Attributes:
        opt_model: parent optical model
        refresh_gui: function to be called on refresh_gui event
        dgm_type: diagram type, 'ht' or 'slp'
    """
    def __init__(self, opt_model, refresh_gui, dgm_type,
                 do_barrel_constraint=False, barrel_constraint=1.0,
                 enable_slide=False, **kwargs):
        self.refresh_gui = refresh_gui
        self.parax_model = opt_model.parax_model
        self.diagram = Diagram(opt_model, dgm_type,
                               do_barrel_constraint=do_barrel_constraint,
                               barrel_constraint=barrel_constraint)
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

    def update_data(self):
        self.artists = []
        self.sys_bbox = self.diagram.update_data(self)
        self.build = 'rebuild'
        return self

    def action_complete(self):
        super().action_complete()
        self.diagram.register_commands((), figure=self)

    def fit_axis_limits(self):
        return self.diagram.fit_axis_limits()
