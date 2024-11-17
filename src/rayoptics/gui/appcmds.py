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

from opticalglass import glassmap as gm
from opticalglass import glassfactory as gfact
from opticalglass import opticalmedium as om

from rayoptics.codev import cmdproc
from rayoptics.optical import obench
from rayoptics.zemax import zmxread

import rayoptics.optical.opticalmodel as opticalmodel
from rayoptics.elem.profiles import Spherical, Conic
import rayoptics.elem.elements as ele
import rayoptics.elem.parttree as pt
from rayoptics.elem import layout
from rayoptics.parax import diagram
from rayoptics.parax import paraxialdesign
from rayoptics.parax.firstorder import specsheet_from_parax_data
from rayoptics.parax.idealimager import ideal_imager_setup
from rayoptics.parax.etendue import create_etendue_dict
from rayoptics.parax.specsheet import (conjugate_types, SpecSheet,
                                       create_specsheet, create_specsheets,
                                       create_specsheet_from_model)
from rayoptics.raytr import vigcalc
from rayoptics.raytr import trace

from rayoptics.gui.appmanager import ModelInfo
from rayoptics.gui.roafile import open_roa

from rayoptics.mpl import interactivefigure
from rayoptics.mpl.interactivelayout import InteractiveLayout
from rayoptics.mpl.axisarrayfigure import Fit
from rayoptics.mpl.axisarrayfigure import (RayFanFigure, SpotDiagramFigure,
                                           WavefrontFigure)

from rayoptics.mpl.analysisplots import FieldCurveFigure, ThirdOrderBarChart

import rayoptics.qtgui.plotview as plotview
from rayoptics.qtgui.idealimagerdialog import IdealImagerDialog
from rayoptics.qtgui.pytablemodel import PyTableModel
from rayoptics.qtgui.plotview import (create_plot_scale_panel,
                                      create_multi_plot_scale_panel,
                                      create_draw_rays_groupbox,
                                      create_diagram_controls_groupbox,
                                      create_diagram_edge_actions_groupbox,
                                      create_diagram_layers_groupbox,
                                      create_2d_figure_toolbar,
                                      )

logger = logging.getLogger(__name__)


def open_model(file_url, info=False, **kwargs):
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
    opm = None
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
    opt_model = opticalmodel.OpticalModel(**kwargs)
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
    opt_model = opticalmodel.OpticalModel(specsheet=specsheet)

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


def create_new_ideal_imager_dialog(**inputs):
    specsheets = {}
    conj_type = (inputs['conjugate_type'] if 'conjugate_type' in inputs
                 else 'finite')
    if 'opt_model' in inputs:
        opt_model = inputs['opt_model']
        specsheet = create_specsheet_from_model(opt_model)
        conj_type = specsheet.conjugate_type
        specsheets[conj_type] = specsheet
        for conj in conjugate_types:
            if conj != conj_type:
                specsheets[conj] = create_specsheet(conj)
    else:
        specsheets = create_specsheets()

    if 'gui_parent' in inputs:
        gui_parent = inputs['gui_parent']
        opt_model = gui_parent.app_manager.model
        iid = IdealImagerDialog(conj_type, specsheets,
                                cmd_fct=gui_parent.handle_ideal_imager_command,
                                parent=gui_parent)

        gui_parent.add_subwindow(iid, ModelInfo(gui_parent.app_manager.model,
                                                update_specsheet,
                                                (iid, opt_model)))
        iid.update_values()
        iid.show()
    else:
        iid = IdealImagerDialog(conj_type, specsheets)
        iid.update_values()
        iid.exec_()


def create_yybar_model():
    opt_model = opticalmodel.OpticalModel()

    # put in minimum calculation defaults
    opt_model.seq_model.gaps[0].thi = 1.0
    opt_model.optical_spec.field_of_view.type = 'OBJ_HT'
    opt_model.optical_spec.field_of_view.set_from_list([0., 1.])

    opt_model.update_model()

    return opt_model


def get_defaults_from_gui_parent(gui_parent):
    if gui_parent:
        refresh_gui = gui_parent.refresh_gui
        is_dark = gui_parent.is_dark
    else:
        refresh_gui = None
        is_dark = True
    return refresh_gui, is_dark


def create_live_layout_view(opt_model, gui_parent=None):
    refresh_gui, is_dark = get_defaults_from_gui_parent(gui_parent)
    fig = InteractiveLayout(opt_model, refresh_gui=refresh_gui,
                            do_draw_frame=True,
                            do_draw_axes=False,
                            do_draw_rays=True,
                            do_paraxial_layout=False,
                            is_dark=is_dark)
    # cmds = create_live_layout_commands(fig)
    cmds = None
    view_width = 880
    view_ht = 660
    title = "Optical Layout"
    panel_fcts = [create_2d_figure_toolbar,
                  create_draw_rays_groupbox,
                  ]
    plotview.create_plot_view(gui_parent, fig, title, view_width, view_ht,
                              add_panel_fcts=panel_fcts, commands=cmds,
                              drop_action=layout.GlassDropAction())


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


def create_paraxial_design_view_v2(opt_model, dgm_type, gui_parent=None):
    from rayoptics.mpl.interactivediagram import InteractiveDiagram
    refresh_gui, is_dark = get_defaults_from_gui_parent(gui_parent)
    fig = InteractiveDiagram(opt_model, dgm_type, refresh_gui=refresh_gui,
                             do_draw_frame=True, do_draw_axes=True,
                             aspect='auto', is_dark=is_dark)
    panel_fcts = [create_2d_figure_toolbar,
                  ]
    if dgm_type == 'ht':
        cmds = create_parax_design_commands(fig)

        panel_fcts.append(create_diagram_controls_groupbox)
        panel_fcts.append(create_diagram_edge_actions_groupbox)
        panel_fcts.append(create_diagram_layers_groupbox)
    elif dgm_type == 'slp':
        panel_fcts.append(create_diagram_layers_groupbox)
        cmds = None        
    else:
        cmds = None

    view_width = 880
    view_ht = 660
    title = "Paraxial Design View"
    plotview.create_plot_view(gui_parent, fig, title, view_width, view_ht,
                              add_panel_fcts=panel_fcts, commands=cmds,
                              drop_action=diagram.GlassDropAction(),
                              context_menu=diagram.context_menu_actions())


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


def create_ray_fan_view(opt_model, data_type, gui_parent=None):
    refresh_gui, is_dark = get_defaults_from_gui_parent(gui_parent)
    fig = RayFanFigure(opt_model, data_type,
                       scale_type=Fit.All_Same,
                       figsize=(5, 4), dpi=100, is_dark=is_dark)
    view_width = 600
    view_ht = 600
    if data_type == "Ray":
        title = "Ray Fan View"
    elif data_type == "OPD":
        title = "OPD Fan View"
    else:
        title = "bad data_type argument"
    panel_fcts = [create_multi_plot_scale_panel]
    plotview.create_plot_view(gui_parent, fig, title, view_width, view_ht,
                              add_panel_fcts=panel_fcts)


def create_ray_grid_view(opt_model, gui_parent=None):
    refresh_gui, is_dark = get_defaults_from_gui_parent(gui_parent)
    num_flds = len(opt_model.optical_spec.field_of_view.fields)

    fig = SpotDiagramFigure(opt_model, scale_type=Fit.All_Same,
                            dpi=100, is_dark=is_dark)
    view_box = 300
    view_width = view_box
    view_ht = num_flds * view_box
    title = "Spot Diagram"
    panel_fcts = [create_multi_plot_scale_panel]
    plotview.create_plot_view(gui_parent, fig, title, view_width, view_ht,
                              add_panel_fcts=panel_fcts)


def create_wavefront_view(opt_model, gui_parent=None):
    refresh_gui, is_dark = get_defaults_from_gui_parent(gui_parent)
    num_flds = len(opt_model.optical_spec.field_of_view.fields)
    num_wvls = len(opt_model.optical_spec.spectral_region.wavelengths)

    fig = WavefrontFigure(opt_model, scale_type=Fit.All_Same,
                          num_rays=32, dpi=100, is_dark=is_dark)
#                                 figsize=(2*num_wvls, 2*num_flds))
    view_box = 300
    view_width = num_wvls * view_box
    view_ht = num_flds * view_box
    title = "Wavefront Map"
    panel_fcts = [create_multi_plot_scale_panel]
    plotview.create_plot_view(gui_parent, fig, title, view_width, view_ht,
                              add_panel_fcts=panel_fcts)


def create_field_curves(opt_model, gui_parent=None):
    refresh_gui, is_dark = get_defaults_from_gui_parent(gui_parent)
    fig = FieldCurveFigure(opt_model, dpi=100, is_dark=is_dark)
    view_width = 600
    view_ht = 600
    title = "Field Curves"
    panel_fcts = [create_plot_scale_panel]
    plotview.create_plot_view(gui_parent, fig, title, view_width, view_ht,
                              add_panel_fcts=panel_fcts)


def create_3rd_order_bar_chart(opt_model, gui_parent=None):
    refresh_gui, is_dark = get_defaults_from_gui_parent(gui_parent)
    fig = ThirdOrderBarChart(opt_model, dpi=100, is_dark=is_dark)
    view_width = 600
    view_ht = 600
    title = "3rd Order Aberrations"
    panel_fcts = [create_plot_scale_panel]
    plotview.create_plot_view(gui_parent, fig, title, view_width, view_ht,
                              add_panel_fcts=panel_fcts)


def create_glass_map_view(opt_model, gui_parent=None):
    refresh_gui, is_dark = get_defaults_from_gui_parent(gui_parent)
    glass_names = set()
    glasses = list()
    for g in opt_model.seq_model.gaps:
        m = g.medium
        if not isinstance(m, om.Air):
            if m.name() not in glass_names:
                glass_names.add(m.name())
                glasses.append(m)
    glass_db = gm.GlassMapDB(glasses, gfact._catalog_list)
    plotview.create_glass_map_view(gui_parent, glass_db)


def update_table_view(table_view, **kwargs):
    table_model = table_view.model()
    table_model.endResetModel()


def create_lens_table_model(seq_model):
    def replace_glass(event, index):
        mime = event.mimeData()
        # comma separated list
        glass_name, catalog_name = mime.text().split(',')
        mat = gfact.create_glass(glass_name, catalog_name)
        seq_model.gaps[index].medium = mat

    colEvalStr = ['.ifcs[{}].interface_type()',
                  '.ifcs[{}].profile_cv',
                  '.gaps[{}].thi',
                  '.gaps[{}].medium.name()',
                  '.ifcs[{}].interact_mode',
                  '.ifcs[{}].surface_od()']
    rowHeaders = seq_model.surface_label_list()
    colHeaders = ['type', 'cv', 'thi', 'medium', 'mode', 'sd']
    colFormats = ['{:s}', '{:12.7g}', '{:12.5g}',
                  '{:s}', '{:s}', '{:12.5g}']
    drop_actions = [None]*len(colHeaders)
    drop_actions[3] = replace_glass
    return PyTableModel(seq_model, '', colEvalStr, rowHeaders,
                        colHeaders, colFormats, True,
                        get_num_rows=seq_model.get_num_surfaces,
                        get_row_headers=seq_model.surface_label_list,
                        drop_actions=drop_actions)


def create_element_table_model(opt_model):
    ele_model = opt_model.ele_model

    def get_row_headers():
        return [str(i) for i in range(ele_model.get_num_elements())]
    colEvalStr = ['.elements[{}].label', '.element_type({})',
                  '.elements[{}].medium_name',
                  '.elements[{}].tfrm[1][1]',
                  '.elements[{}].tfrm[1][2]',
                  '.elements[{}].reference_idx()']

    rowHeaders = get_row_headers()
    colHeaders = ['label', 'type', 'medium', 'y', 'z', 'idx']
    colFormats = ['{:s}', '{:s}', '{:s}', '{:12.5g}', '{:12.5g}',
                  '{:d}']
    return PyTableModel(ele_model, '', colEvalStr, rowHeaders,
                        colHeaders, colFormats, True,
                        get_num_rows=ele_model.get_num_elements,
                        get_row_headers=get_row_headers)


def create_ray_table_model(opt_model, ray):
    colEvalStr = ['[{}].p[0]', '[{}].p[1]', '[{}].p[2]',
                  '[{}].d[0]', '[{}].d[1]', '[{}].d[2]',
                  '[{}].dst']
    seq_model = opt_model.seq_model
    rowHeaders = seq_model.surface_label_list()
    colHeaders = ['x', 'y', 'z', 'l', 'm', 'n', 'length']
    colFormats = ['{:12.5g}', '{:12.5g}', '{:12.5g}', '{:9.6f}',
                  '{:9.6f}', '{:9.6f}', '{:12.5g}']
    return PyTableModel(ray, '', colEvalStr, rowHeaders,
                        colHeaders, colFormats, False,
                        get_num_rows=seq_model.get_num_surfaces,
                        get_row_headers=seq_model.surface_label_list)


def create_parax_table_model(opt_model):
    rootEvalStr = ".analysis_results['parax_data']"
    colEvalStr = ['[0][{}][0]', '[0][{}][1]', '[0][{}][2]',
                  '[1][{}][0]', '[1][{}][1]', '[1][{}][2]']
    seq_model = opt_model.seq_model
    rowHeaders = seq_model.surface_label_list()
    colHeaders = ['y', 'u', 'i', 'y-bar', 'u-bar', 'i-bar']
    colFormats = ['{:12.5g}', '{:9.6f}', '{:9.6f}', '{:12.5g}',
                  '{:9.6f}', '{:9.6f}']
    return PyTableModel(opt_model, rootEvalStr, colEvalStr, rowHeaders,
                        colHeaders, colFormats, False,
                        get_num_rows=seq_model.get_num_surfaces,
                        get_row_headers=seq_model.surface_label_list)


def create_parax_model_table(opt_model):
    colEvalStr = ['.ax[{}][0]', '.pr[{}][0]', '.ax[{}][1]', '.pr[{}][1]',
                  '.sys[{}][0]', '.sys[{}][1]', '.sys[{}][2]', '.sys[{}][3]']
    seq_model = opt_model.seq_model
    parax_model = opt_model.parax_model
    rowHeaders = seq_model.surface_label_list()
    colHeaders = ['y', 'y-bar', 'nu', 'nu-bar',
                  'pwr', 'tau', 'n after', 'mode']
    colFormats = ['{:12.5g}', '{:12.5g}', '{:9.6f}', '{:9.6f}',
                  '{:12.7g}', '{:12.5g}', '{:7.4f}', '{:s}']
    return PyTableModel(parax_model, '', colEvalStr, rowHeaders,
                        colHeaders, colFormats, True,
                        get_num_rows=parax_model.get_num_nodes,
                        get_row_headers=seq_model.surface_label_list)
