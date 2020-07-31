#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2020 Michael J. Hayford
"""Read an .roa file and return a OpticalModel instance

.. Created on Sun Jul 12 22:25:37 2020

.. codeauthor: Michael J. Hayford
"""

import json_tricks


module_repl_050 = {
    'rayoptics.optical.elements': 'rayoptics.elem.elements',
    'rayoptics.optical.profiles': 'rayoptics.elem.profiles',
    'rayoptics.optical.surface': 'rayoptics.elem.surface',
    'rayoptics.optical.transform': 'rayoptics.elem.transform',
    'rayoptics.gui.layout': 'rayoptics.elem.layout',
    'rayoptics.optical.etendue': 'rayoptics.parax.etendue',
    'rayoptics.optical.firstorder': 'rayoptics.parax.firstorder',
    'rayoptics.optical.idealimager': 'rayoptics.parax.idealimager',
    'rayoptics.optical.paraxialdesign': 'rayoptics.parax.paraxialdesign',
    'rayoptics.optical.specsheet': 'rayoptics.parax.specsheet',
    'rayoptics.optical.thirdorder': 'rayoptics.parax.thirdorder',
    'rayoptics.gui.diagram': 'rayoptics.parax.diagram',
    'rayoptics.optical.analyses': 'rayoptics.raytr.analyses',
    'rayoptics.optical.opticalspec': 'rayoptics.raytr.opticalspec',
    'rayoptics.optical.raytrace': 'rayoptics.raytr.raytrace',
    'rayoptics.optical.sampler': 'rayoptics.raytr.sampler',
    'rayoptics.optical.trace': 'rayoptics.raytr.trace',
    'rayoptics.optical.traceerror': 'rayoptics.raytr.traceerror',
    'rayoptics.optical.doe': 'rayoptics.oprops.doe',
    'rayoptics.optical.gap': 'rayoptics.seq.gap',
    'rayoptics.optical.interface': 'rayoptics.seq.interface',
    'rayoptics.optical.medium': 'rayoptics.seq.medium',
    'rayoptics.optical.sequential': 'rayoptics.seq.sequential',
    'rayoptics.optical.thinlens': 'rayoptics.oprops.thinlens',
    'rayoptics.optical.twoconicmirrors': 'rayoptics.seq.twoconicmirrors',
    }


def preprocess_roa(file_name, str_replacements):
    """Replace any old module references with the v0.5 ones."""
    with open(file_name, 'r') as f:
        contents = f.read()
    for old, new in str_replacements.items():
        contents = contents.replace(old, new)
    return contents


def open_roa(file_name, mapping=None):
    """ open a ray-optics file and populate an optical model with the data

    Args:
        file_name (str): a filename with a .roa extension
        mapping: dict mapping old modules to new. If None, use module_repl_050

    Returns:
        if successful, an OpticalModel instance, otherwise, None
    """
    opm = None
    str_replacements = module_repl_050 if mapping is None else mapping
    contents = preprocess_roa(file_name, str_replacements)
    obj_dict = json_tricks.loads(contents)
    if 'optical_model' in obj_dict:
        opm = obj_dict['optical_model']
        opm.sync_to_restore()
    return opm
