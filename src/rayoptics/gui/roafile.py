#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2020 Michael J. Hayford
"""Read an .roa file and return a OpticalModel instance

.. Created on Sun Jul 12 22:25:37 2020

.. codeauthor: Michael J. Hayford
"""

import json_tricks
import logging
from packaging import version
from pathlib import Path

logger = logging.getLogger(__name__)

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
    """Read and preprocess raw roa text file, returning preprocessed text."""

    # Replace any old module references with the v0.5 ones.
    with open(file_name, 'r') as f:
        contents = f.read()
    for old, new in str_replacements.items():
        contents = contents.replace(old, new)

    return contents


def postprocess_roa(opt_model, file_path, **kwargs):
    """Post processing for raw optical_model, including sync_to_restore. """

    # Force rebuild of ele_model for pre-0.8.5 models.
    old_version = False
    rebuild_from_seq = kwargs.get('rebuild_from_seq', False)
    cutoff_version = kwargs.get('cutoff_version', "0.8.5")
    if (not hasattr(opt_model, 'ro_version') or rebuild_from_seq == True or
        version.parse(opt_model.ro_version) < version.parse(cutoff_version)):
        old_version = True
        opt_model.ele_model.elements = []
        if hasattr(opt_model, 'part_tree'):
            delattr(opt_model, 'part_tree')

    for i, g in enumerate(opt_model.seq_model.gaps):
        if hasattr(g.medium, 'convert_to_OG'):
            g.medium = g.medium.convert_to_OG()

    opt_model.sync_to_restore()
    save_updated_version = kwargs.get('save_updated_version', True)
    if old_version and save_updated_version:
        save_updated_roa(file_path, opt_model)
    return opt_model


def save_updated_roa(file_path: Path, opt_model):
    """ rename file_path to file_path_version# and save new file_path version """

    forig = file_path.stem
    cur_vers =  version.parse(opt_model.ro_version).public
    fname_archive = forig + "_v" + cur_vers.replace('.', '')
    f_archive = str(file_path).replace(forig, fname_archive)
    file_path.rename(Path(f_archive))
    msg1 = f" Archived original file as {Path(f_archive).name}"
    print(msg1)
    logger.info(msg1)

    opt_model.save_model(file_path, version="0.9.0")
    new_vers =  version.parse(opt_model.ro_version).public
    msg2 = f" Updated {file_path.name} from version {cur_vers} -> {new_vers}."
    print(msg2)
    logger.info(msg2)


def open_roa(file_name, mapping=None, **kwargs):
    """ open a ray-optics file and populate an optical model with the data

    Args:
        file_name (str): a filename with a .roa extension
        mapping: dict mapping old modules to new. If None, use module_repl_050

    Returns:
        if successful, an OpticalModel instance, otherwise, None
    """
    opt_model = None
    str_replacements = module_repl_050 if mapping is None else mapping
    contents = preprocess_roa(file_name, str_replacements)
    obj_dict = json_tricks.loads(contents)
    if 'optical_model' in obj_dict:
        opt_model = obj_dict['optical_model']
        postprocess_roa(opt_model, file_name, **kwargs)
    return opt_model
