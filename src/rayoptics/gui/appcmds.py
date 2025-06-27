#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" generic ray optics commands for creating plots and tables

.. Created on Thu Nov  8 21:21:57 2018

.. codeauthor: Michael J. Hayford
"""

import logging
import math
import pathlib
from typing import Optional

from rayoptics.gui.roafile import open_roa
from rayoptics.codev import cmdproc
from rayoptics.optical import obench
from rayoptics.zemax import zmxread

from rayoptics.optical.opticalmodel import OpticalModel
from rayoptics.elem.profiles import Spherical, Conic
import rayoptics.elem.elements as ele
from rayoptics.elem import layout
from rayoptics.parax import paraxialdesign
from rayoptics.parax.firstorder import specsheet_from_parax_data
from rayoptics.parax.idealimager import ideal_imager_setup
from rayoptics.parax.etendue import create_etendue_dict
from rayoptics.parax.specsheet import SpecSheet
from rayoptics.raytr import vigcalc
from rayoptics.raytr import trace

from rayoptics.mpl import interactivefigure


logger = logging.getLogger(__name__)


def open_model(file_url, info=False, **kwargs) -> Optional[OpticalModel] | tuple[OpticalModel, tuple[dict, dict]]:
    """ open a file or url and populate an optical model with the data

    Args:
        file_url (str): a filename or url of a supported file type

            - .roa - a rayoptics JSON encoded file
            - .seq - a CODE V (TM) sequence file
            - .zmx - a Zemax (TM) lens file
            - a URL from the www.photonstophotos.net OpticalBench database
        info (bool): if true, return an info tuple with import statistics
        kwargs (dict): keyword args passed to the reader functions

    Returns:
        if successful, an OpticalModel instance, otherwise, None
    """
    file_url_pth = pathlib.Path(file_url)
    file_extension = file_url_pth.suffix.lower()
    opm: OpticalModel|None = None
    if file_extension == '.roa':
        # if we have a rayoptics file, we just read it
        opm = open_roa(file_url_pth, **kwargs)
    else:
        # if we're importing another program's file, collect import info
        if 'www.photonstophotos.net' in str(file_url_pth):
            opm, import_info = obench.read_obench_url(file_url, **kwargs)
        elif file_extension == '.seq':
            opm, import_info = cmdproc.read_lens(file_url_pth, **kwargs)
        elif file_extension == '.zmx':
            opm, import_info = zmxread.read_lens_file(file_url_pth, **kwargs)
        
        # At this point we have a complete opt_model; 
        #  all of the above call opm.update_model()
        if info:
            return opm, import_info
    return opm


def create_empty_model(**kwargs):
    """ factory function returns an instance of OpticalModel """
    opt_model = OpticalModel(**kwargs)
    return opt_model


def create_new_model():
    return create_new_optical_system()


def create_new_optical_system(efl=10.0, epd=1, fov=1.0):
    imager_inputs = {'s': -math.inf, 'f': efl}
    inf_conj_imgr = ideal_imager_setup(**imager_inputs)

    ei = create_etendue_dict()
    ei['field']['object']['angle'] = fov
    ei['aperture']['object']['epd'] = epd
    ssi = SpecSheet('infinite',
                    imager=inf_conj_imgr,
                    imager_inputs=imager_inputs,
                    etendue_inputs=ei)

    opt_model = create_new_optical_model_from_specsheet(ssi)

    sr = opt_model.optical_spec.spectral_region
    sr.set_from_list([('F', 1), ('d', 1), ('C', 1)])
    sr.reference_wvl = 1

    return opt_model


def create_new_optical_model_from_specsheet(specsheet):
    """ create an OpticalModel with a basic thinlens model, given specsheet """
    opt_model = OpticalModel(specsheet=specsheet)

    # enter a basic thinlens model for the given specsheet
    imager = specsheet.imager
    if specsheet.conjugate_type == 'finite':
        opt_model.seq_model.gaps[0].thi = -imager.s
    else:
        opt_model.seq_model.gaps[0].thi = 1.0e10

    opt_model.add_thinlens(power=1/imager.f, indx=1.5, idx=0, t=imager.sp)

    opt_model.update_model()

    return opt_model


def update_specsheet(iid, opt_model):
    specsheet = opt_model.specsheet
    specsheet_from_parax_data(opt_model, specsheet)
    iid.specsheet_dict[specsheet.conjugate_type] = specsheet
    iid.update_values()


def create_yybar_model():
    opt_model = OpticalModel()

    # put in minimum calculation defaults
    opt_model.seq_model.gaps[0].thi = 1.0
    opt_model.optical_spec.field_of_view.type = 'OBJ_HT'
    opt_model.optical_spec.field_of_view.set_from_list([0., 1.])

    opt_model.update_model()

    return opt_model


def create_live_layout_commands(fig):
    lo = fig.layout
    cmds = []
    # Add thin lens
    cmds.append(('Add Thin Lens', (lo.register_commands, (),
                 {'apply_fct': layout.add_thinlens})))
    # Add lens
    cmds.append(('Add Lens', (lo.register_commands, (),
                 {'apply_fct': layout.add_lens})))
    # Add doublet
    cmds.append(('Add Cemented Doublet', (lo.register_commands, (),
                 {'apply_fct': layout.add_doublet})))
    # Add mirror
    cmds.append(('Add Mirror', (lo.register_commands, (),
                 {'apply_fct': layout.add_mirror,
                  'profile': Spherical})))
    cmds.append(('Add Conic Mirror', (lo.register_commands, (),
                 {'apply_fct': layout.add_conic,
                  'profile': Conic})))
    cmds.append(('Add Parabola', (lo.register_commands, (),
                 {'apply_fct': layout.add_conic,
                  'profile': Conic,
                  'cc': -1.0})))

    return cmds


def create_parax_design_commands(fig):
    cmds = []
    dgm = fig.diagram
    # initialize dgm with a Select command
    args = tuple()
    kwargs = {'figure': fig,
              }
    dgm.register_commands(*args, **kwargs)

    # draw polyline to define an initial model
    args = (interactivefigure.snap_to_grid_fct(0.05), True)
    kwargs = {'factory': ele.create_thinlens,
              'opt_model': dgm.opt_model,
              'gui_fct': interactivefigure.enter_polyline,
              }
    cmds.append(('Sketch Diagram', (paraxialdesign.nodes_to_new_model, 
                                          args, kwargs)))

    # Select an existing point
    cmds.append(('Select', (dgm.register_commands, (), 
                            {'figure': fig,
                             })))

    # Add thin lens
    cmds.append(('Add Thin Lens',
                 (dgm.register_add_replace_element, (),
                  {'node_init': ele.create_thinlens,
                   'factory': ele.create_thinlens,
                   'interact_mode': 'transmit'})))
    # Add lens
    kwargs = {'node_init': ele.create_lens_from_dgm,
              'factory': ele.create_lens_from_dgm,
              'interact_mode': 'transmit',
              }
    cmds.append(('Add Lens', (dgm.register_add_replace_element, (), kwargs)))

    # Add doublet
    kwargs = {'node_init': ele.create_cemented_doublet,
              'factory': ele.create_cemented_doublet,
              'interact_mode': 'transmit'}
    cmds.append(('Add Cemented Doublet', (dgm.register_add_replace_element, 
                                          (), kwargs)))

    # Add mirror
    cmds.append(('Add Mirror',
                 (dgm.register_add_replace_element, (),
                  {'node_init': ele.create_mirror,
                   'factory': ele.create_mirror,
                   'interact_mode': 'reflect'})))

    # draw polyline
    # cmds.append(('Draw Polyline', (interactivefigure.enter_polyline, 
    #                                       (), {})))

    # Replace with file
    pth = pathlib.Path(__file__).resolve()
    try:
        rayoptics_pos = pth.parts.index('rayoptics')
    except ValueError:
        logger.debug("Can't find rayoptics: path is %s", pth)
    else:
        # models_dir = rayoptics/models
        models_dir = pathlib.Path(*pth.parts[:rayoptics_pos+1]) / 'models'
        filepath = models_dir / 'Sasian Triplet.roa'

        def cff(**kwargs):
            return ele.create_from_file(filepath, **kwargs)

        cmds.append(('Sasian Triplet',
                     (dgm.register_add_replace_element, (),
                      {'filename': filepath,
                       'node_init': ele.create_thinlens,
                       'factory': cff,
                       'interact_mode': 'transmit'})))
    finally:
        return cmds


def set_vignetting(opt_model, gui_parent=None, **kwargs):
    """ From existing fields and clear apertures, calculate vignetting. """
    vigcalc.set_vig(opt_model, **kwargs)
    if gui_parent is None:
        opt_model.update_model(src_model=opt_model['seq_model'])
    else:
        gui_parent.refresh_gui(src_model=opt_model['seq_model'])


def set_apertures(opt_model, gui_parent=None):
    """ From existing fields and vignetting, calculate clear apertures. """
    vigcalc.set_ape(opt_model)
    if gui_parent is None:
        opt_model.update_model(src_model=opt_model['seq_model'])
    else:
        gui_parent.refresh_gui(src_model=opt_model['seq_model'])


def set_pupil(opt_model, gui_parent=None):
    """ From existing stop size, calculate pupil spec and vignetting. """
    vigcalc.set_pupil(opt_model)
    if gui_parent is None:
        opt_model.update_model(src_model=opt_model['seq_model'])
    else:
        gui_parent.refresh_gui(src_model=opt_model['seq_model'])


def refocus(opt_model, gui_parent=None):
    """ Compute a focus shift bringing the axial marginal ray to zero. """
    focus_shift = trace.refocus(opt_model)
    opt_model['optical_spec']['focus'].focus_shift = focus_shift

    if gui_parent is None:
        opt_model.update_model(src_model=opt_model['optical_spec'])
    else:
        gui_parent.refresh_gui(src_model=opt_model['optical_spec'])

def set_paraxial_focus(opt_model, gui_parent=None):
    """ Set the image at the paraxial image point using the final gap."""
    opt_model['parax_model'].set_paraxial_focus(opt_model['seq_model'])
    if gui_parent is None:
        opt_model.update_model(src_model=opt_model['optical_spec'])
    else:
        gui_parent.refresh_gui(src_model=opt_model['optical_spec'])
