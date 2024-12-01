#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Qt Gui commands for creating plots and tables

.. Created on Sun Nov 24 21:21:57 2024

.. codeauthor: Michael J. Hayford
"""

import logging

from opticalglass import glassmap as gm
from opticalglass import glassfactory as gfact
from opticalglass import opticalmedium as om

from rayoptics.elem import layout
from rayoptics.parax import diagram

from rayoptics.parax.specsheet import (conjugate_types,
                                       create_specsheet, 
                                       create_specsheets,
                                       create_specsheet_from_model)

from rayoptics.gui import appcmds
from rayoptics.gui.appmanager import ModelInfo

from rayoptics.mpl.interactivelayout import InteractiveLayout
from rayoptics.mpl.axisarrayfigure import Fit
from rayoptics.mpl.axisarrayfigure import (RayFanFigure, 
                                           SpotDiagramFigure,
                                           WavefrontFigure)
from rayoptics.mpl.analysisplots import (FieldCurveFigure, 
                                         ThirdOrderBarChart)

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
                                                appcmds.update_specsheet,
                                                (iid, opt_model)))
        iid.update_values()
        iid.show()
    else:
        iid = IdealImagerDialog(conj_type, specsheets)
        iid.update_values()
        iid.exec_()


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
    # cmds = appcmds.create_live_layout_commands(fig)
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


def create_paraxial_design_view_v2(opt_model, dgm_type, gui_parent=None):
    from rayoptics.mpl.interactivediagram import InteractiveDiagram
    refresh_gui, is_dark = get_defaults_from_gui_parent(gui_parent)
    fig = InteractiveDiagram(opt_model, dgm_type, refresh_gui=refresh_gui,
                             do_draw_frame=True, do_draw_axes=True,
                             aspect='auto', is_dark=is_dark)
    panel_fcts = [create_2d_figure_toolbar,
                  ]
    if dgm_type == 'ht':
        cmds = appcmds.create_parax_design_commands(fig)

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
