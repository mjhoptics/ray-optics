#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" generic ray optics commands for creating plots and tables

.. Created on Thu Nov  8 21:21:57 2018

.. codeauthor: Michael J. Hayford
"""

from rayoptics.optical.opticalmodel import OpticalModel, open_model

from rayoptics.mpl.lenslayoutfigure import LensLayoutFigure
from rayoptics.mpl.axisarrayfigure import Fit
from rayoptics.mpl.axisarrayfigure import (RayFanFigure, SpotDiagramFigure,
                                           WavefrontFigure)

from rayoptics.mpl.analysisplots import FieldCurveFigure
from rayoptics.mpl.paraxdgnfigure import (ParaxialDesignFigure,
                                          create_parax_design_commands)
import rayoptics.qtgui.plotview as plotview
from rayoptics.qtgui.pytablemodel import PyTableModel


def create_new_model():
    opt_model = OpticalModel()

    # put in minimum calculation defaults
    opt_model.seq_model.gaps[0].thi = 1.0
    opt_model.optical_spec.field_of_view.type = 'OBJ_HT'
    opt_model.optical_spec.field_of_view.set_from_list([0., 1.])

    opt_model.update_model()

    return opt_model


def create_lens_layout_view(opt_model, gui_parent=None):
    fig = LensLayoutFigure(opt_model)
    view_width = 660
    view_ht = 440
    title = "Lens Layout View"
    plotview.create_plot_view(gui_parent, fig, title, view_width, view_ht,
                              add_scale_panel=False)


def create_paraxial_design_view(opt_model, dgm_type, gui_parent=None):
    fig = ParaxialDesignFigure(opt_model, gui_parent.refresh_gui, dgm_type,
                               figsize=(5, 4))
    cmds = create_parax_design_commands(fig)
    view_width = 650
    view_ht = 500
    title = "Paraxial Design View"
    plotview.create_cmdplot_view(gui_parent, fig, title, view_width, view_ht,
                                 commands=cmds)


def create_ray_fan_view(opt_model, data_type, gui_parent=None):
    fig = RayFanFigure(opt_model, data_type,
                       scale_type=Fit.All_Same,
                       figsize=(5, 4), dpi=100)
    view_width = 600
    view_ht = 600
    if data_type == "Ray":
        title = "Ray Fan View"
    elif data_type == "OPD":
        title = "OPD Fan View"
    else:
        title = "bad data_type argument"
    plotview.create_plot_view(gui_parent, fig, title, view_width, view_ht)


def create_ray_grid_view(opt_model, gui_parent=None):
    num_flds = len(opt_model.optical_spec.field_of_view.fields)

    fig = SpotDiagramFigure(opt_model, scale_type=Fit.All_Same,
                            num_rays=11, dpi=100)
    view_box = 300
    view_width = view_box
    view_ht = num_flds * view_box
    title = "Spot Diagram"
    plotview.create_plot_view(gui_parent, fig, title, view_width, view_ht)


def create_wavefront_view(opt_model, gui_parent=None):
    num_flds = len(opt_model.optical_spec.field_of_view.fields)
    num_wvls = len(opt_model.optical_spec.spectral_region.wavelengths)

    fig = WavefrontFigure(opt_model, scale_type=Fit.All_Same,
                          num_rays=32, dpi=100)
#                                 figsize=(2*num_wvls, 2*num_flds))
    view_box = 300
    view_width = num_wvls * view_box
    view_ht = num_flds * view_box
    title = "Wavefront Map"
    plotview.create_plot_view(gui_parent, fig, title, view_width, view_ht)


def create_field_curves(opt_model, gui_parent=None):
    fig = FieldCurveFigure(opt_model, dpi=100)
    view_width = 600
    view_ht = 600
    title = "Field Curves"
    plotview.create_plot_view(gui_parent, fig, title, view_width, view_ht)


def update_table_view(table_view):
    table_model = table_view.model()
    table_model.endResetModel()


def create_lens_table_model(seq_model):
    colEvalStr = ['.ifcs[{}].interface_type()',
                  '.ifcs[{}].profile_cv()',
                  '.ifcs[{}].surface_od()', '.gaps[{}].thi',
                  '.gaps[{}].medium.name()', '.ifcs[{}].refract_mode']
    rowHeaders = seq_model.surface_label_list()
    colHeaders = ['type', 'cv', 'sd', 'thi', 'medium', 'mode']
    colFormats = ['{:s}', '{:12.7g}', '{:12.5g}', '{:12.5g}',
                  '{:s}', '{:s}']
    return PyTableModel(seq_model, '', colEvalStr, rowHeaders,
                        colHeaders, colFormats, True,
                        get_num_rows=seq_model.get_num_surfaces,
                        get_row_headers=seq_model.surface_label_list)


def create_ray_table_model(opt_model, ray):
    colEvalStr = ['[{}][0][0]', '[{}][0][1]', '[{}][0][2]',
                  '[{}][1][0]', '[{}][1][1]', '[{}][1][2]',
                  '[{}][2]']
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
    rootEvalStr = ".optical_spec.parax_data"
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
