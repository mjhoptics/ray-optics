#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Ray Optics GUI Application

Relies on PyQt5

Created on Mon Feb 12 09:24:01 2018

@author: Michael J. Hayford
"""

import sys
import logging
from pathlib import Path

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import (QApplication, QAction, QMainWindow, QMdiArea,
                             QMdiSubWindow, QTextEdit, QFileDialog, QTableView,
                             QHBoxLayout, QVBoxLayout, QWidget, QGraphicsView,
                             QLineEdit, QGraphicsScene, QRadioButton,
                             QGroupBox)
from PyQt5.QtCore import pyqtSlot

import rayoptics.optical.opticalmodel as optm
import rayoptics.gui.plotcanvas as plotter
import rayoptics.gui.mpl.axisarrayfigure as aaf
import rayoptics.gui.mpl.paraxdgnfigure as pdf
from rayoptics.gui.mpl.lenslayoutfigure import LensLayoutFigure
import rayoptics.gui.pytablemodel as tbl
import rayoptics.gui.graphicsitems as gitm
from rayoptics.gui.appmanager import ModelInfo, AppManager


class MainWindow(QMainWindow):
    count = 0

    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.mdi = QMdiArea()
        self.setCentralWidget(self.mdi)

        self.app_manager = AppManager(None)
        self.mdi.subWindowActivated.connect(self.app_manager.
                                            on_window_activated)

        self.left = 100
        self.top = 50
        self.width = 1200
        self.height = 800
        self.setGeometry(self.left, self.top, self.width, self.height)

        bar = self.menuBar()

        file = bar.addMenu("File")
        file.addAction("New")
        file.addAction("Open...")
        file.addSeparator()
        file.addAction("Save")
        file.addAction("Save As...")
        file.triggered[QAction].connect(self.file_action)
        view = bar.addMenu("View")
        view.addAction("Lens Table")
        view.addAction("Lens View")
        view.addSeparator()
        view.addAction("Paraxial Height View")
        view.addAction("Paraxial Slope View")
        view.addAction("Paraxial Ray Table")
        view.addAction("Ray Table")
        view.addSeparator()
        view.addAction("Ray Fans")
        view.addAction("OPD Fans")
        view.addAction("Spot Diagram")
        view.addAction("Wavefront Map")
        view.addSeparator()
        view.triggered[QAction].connect(self.view_action)
        wnd = bar.addMenu("Window")
        wnd.addAction("Cascade")
        wnd.addAction("Tiled")
        wnd.triggered[QAction].connect(self.window_action)

        self.setWindowTitle("Ray Optics")
        self.show()

        pth = Path(__file__).resolve()
        try:
            root_pos = pth.parts.index('rayoptics')
        except ValueError:
            logging.debug("Can't find rayoptics: path is %s", pth)
        else:
            path = Path(*pth.parts[:root_pos+1])
#            self.open_file(path / "codev/test/asp46.seq")
#            self.open_file(path / "codev/test/paraboloid.seq")
#            self.open_file(path / "codev/test/paraboloid_f8.seq")
#            self.open_file(path / "codev/test/schmidt.seq")
#            self.open_file(path / "codev/test/questar35.seq")
#            self.open_file(path / "codev/test/conic_mirror.seq")
#            self.open_file(path / "codev/test/rc_f16.seq")
            self.open_file(path / "codev/test/ag_dblgauss.seq")
#            self.open_file(path / "codev/test/landscape_lens.seq")

#        try:
#            root_pos = pth.parts.index('ray-optics')
#        except ValueError:
#            logging.debug("Can't find ray-optics: path is %s", pth)
#        else:
#            path = Path(*pth.parts[:root_pos+1])
#            self.open_file(path / "test/TwoMirror.roa")
#            self.open_file(path / "test/Sasian Triplet.roa")
#            self.open_file(path / "test/singlet_f5.roa")
#            self.open_file(path / "test/Ritchey_Chretien.roa")

    def add_subwindow(self, widget, model_info):
            sub_wind = self.mdi.addSubWindow(widget)
            self.app_manager.add_window(sub_wind, model_info)
            MainWindow.count += 1
            return sub_wind

    def delete_subwindow(self, sub_wind):
            self.app_manager.delete_window(sub_wind)
            self.mdi.removeSubWindow(sub_wind)
            MainWindow.count -= 1

    def initial_window_offset(self):
        offset_x = 50
        offset_y = 25
        orig_x = (MainWindow.count - 1)*offset_x
        orig_y = (MainWindow.count - 1)*offset_y
        return orig_x, orig_y

    def file_action(self, q):
        if q.text() == "New":
            opt_model = optm.OpticalModel()
            sub = self.add_subwindow(QTextEdit(),
                                     ModelInfo(opt_model, None, None))
            sub.setWindowTitle("subwindow"+str(MainWindow.count))
            sub.show()

        if q.text() == "Open...":
            options = QFileDialog.Options()
            # options |= QFileDialog.DontUseNativeDialog
            fileName, _ = QFileDialog.getOpenFileName(self,
                          "QFileDialog.getOpenFileName()",
                          "",
                          "CODE V Files (*.seq);;Ray-Optics Files (*.roa)",
                          options=options)
            if fileName:
                logging.debug("open file: %s", fileName)
                self.open_file(fileName)

        if q.text() == "Save As...":
            options = QFileDialog.Options()
            # options |= QFileDialog.DontUseNativeDialog
            fileName, _ = QFileDialog.getSaveFileName(self,
                          "QFileDialog.getSaveFileName()",
                          "",
                          "Ray-Optics Files (*.roa);;All Files (*)",
                          options=options)
            if fileName:
                logging.debug("save file: %s", fileName)
                self.save_file(fileName)

    def open_file(self, file_name):
        self.cur_filename = file_name
        self.app_manager.model = optm.open_model(file_name)
        self.is_changed = True
        self.create_lens_table()
#        self.create_lens_layout_view()
        self.create_2D_lens_view()

    def save_file(self, file_name):
        self.app_manager.model.save_model(file_name)
        self.cur_filename = file_name
        self.is_changed = False

    def view_action(self, q):
        if q.text() == "Lens Table":
            self.create_lens_table()

        if q.text() == "Lens View":
            self.create_lens_layout_view()
#            self.create_2D_lens_view()

        if q.text() == "Ray Fans":
            self.create_ray_fan_view("Ray")

        if q.text() == "OPD Fans":
            self.create_ray_fan_view("OPD")

        if q.text() == "Spot Diagram":
            self.create_ray_grid_view()

        if q.text() == "Wavefront Map":
            self.create_wavefront_view()

        if q.text() == "Paraxial Height View":
            self.create_paraxial_design_view(pdf.ht_dgm)

        if q.text() == "Paraxial Slope View":
            self.create_paraxial_design_view(pdf.slp_dgm)

        if q.text() == "Paraxial Ray Table":
            root = self.app_manager.model
            rootEvalStr = ".seq_model.optical_spec.parax_data"
            colEvalStr = ['[0][{}][0]', '[0][{}][1]', '[0][{}][2]',
                          '[1][{}][0]', '[1][{}][1]', '[1][{}][2]']
            rowHeaders = self.app_manager.model.seq_model.surface_label_list()
            colHeaders = ['y', 'u', 'i', 'y-bar', 'u-bar', 'i-bar']
            colFormats = ['{:12.5g}', '{:9.6f}', '{:9.6f}', '{:12.5g}',
                          '{:9.6f}', '{:9.6f}']
            model = tbl.PyTableModel(root, rootEvalStr, colEvalStr, rowHeaders,
                                     colHeaders, colFormats, False)
            self.create_table_view(model, "Paraxial Ray Table")

        if q.text() == "Ray Table":
            self.create_ray_table()

    def window_action(self, q):
        if q.text() == "Cascade":
            self.mdi.cascadeSubWindows()

        if q.text() == "Tiled":
            self.mdi.tileSubWindows()

    def create_lens_table(self):
        seq_model = self.app_manager.model.seq_model
        colEvalStr = ['.ifcs[{}].interface_type()',
                      '.ifcs[{}].profile_cv()',
                      '.ifcs[{}].surface_od()', '.gaps[{}].thi',
                      '.gaps[{}].medium.name()', '.ifcs[{}].refract_mode']
        rowHeaders = seq_model.surface_label_list()
        colHeaders = ['type', 'cv', 'sd', 'thi', 'medium', 'mode']
        colFormats = ['{:s}', '{:12.7g}', '{:12.5g}', '{:12.5g}',
                      '{:s}', '{:s}']
        model = tbl.PyTableModel(seq_model, '', colEvalStr, rowHeaders,
                                 colHeaders, colFormats, True)
        self.create_table_view(model, "Surface Data Table")

    def create_ray_table(self):
        sm = self.app_manager.model.seq_model
        osp = sm.optical_spec
        pupil = [0., 0.]
        fi = 1
        wl = osp.spectral_region.reference_wvl
        fld, wvl, foc = osp.lookup_fld_wvl_focus(fi, wl)
        ray, ray_op, wvl, opd = osp.trace_with_opd(sm, pupil, fld, wvl, foc)

        colEvalStr = ['[{}][0][0]', '[{}][0][1]', '[{}][0][2]',
                      '[{}][1][0]', '[{}][1][1]', '[{}][1][2]',
                      '[{}][2]']
        rowHeaders = sm.surface_label_list()
        colHeaders = ['x', 'y', 'z', 'l', 'm', 'n', 'length']
        colFormats = ['{:12.5g}', '{:12.5g}', '{:12.5g}', '{:9.6f}',
                      '{:9.6f}', '{:9.6f}', '{:12.5g}']
        model = tbl.PyTableModel(ray, '', colEvalStr, rowHeaders,
                                 colHeaders, colFormats, False)
        self.create_table_view(model, "Ray Table")

    def create_2D_lens_view(self):
        scene2d = QGraphicsScene()
        self.create_element_model(scene2d)
        self.create_ray_model(scene2d)
        scene2d.setBackgroundBrush(QColor(237, 243, 254))  # light blue
        sceneRect2d = scene2d.sceneRect()

        # construct the top level widget
        widget = QWidget()
        # construct the top level layout
        layout = QVBoxLayout(widget)

        # set the layout on the widget
        widget.setLayout(layout)

        sub = self.add_subwindow(widget, ModelInfo(self.app_manager.model,
                                         MainWindow.update_2D_lens_view,
                                         scene2d))
        sub.setWindowTitle("2D Lens View")
        view_width = 660
        view_ht = 440
        view_ratio = view_width/view_ht
        orig_x, orig_y = self.initial_window_offset()
        sub.setGeometry(orig_x, orig_y, view_width, view_ht)

        self.gview2d = QGraphicsView(scene2d)
        scene_ratio = sceneRect2d.width()/sceneRect2d.height()
        oversize_fraction = 1.2
        if scene_ratio > view_ratio:
            view_scale = view_width/(oversize_fraction*sceneRect2d.width())
        else:
            view_scale = view_ht/(oversize_fraction*sceneRect2d.height())

        self.gview2d.scale(view_scale, view_scale)
        layout.addWidget(self.gview2d)

        sub.show()

    def update_2D_lens_view(scene2d):
        for gi in scene2d.items():
            gi.prepareGeometryChange()
            gi.update_shape()

    def create_element_model(self, gscene):
        ele_model = self.app_manager.model.ele_model
        ele_model.elements_from_sequence(self.app_manager.model.seq_model)
        for e in ele_model.elements:
            ge = gitm.OpticalElement(e)
            gscene.addItem(ge)

    def create_ray_model(self, gscene, start_surf=1):
        seq_model = self.app_manager.model.seq_model

        img_dist = abs(seq_model.optical_spec.parax_data[2].img_dist)
        start_offset = 0.05*(gscene.sceneRect().width() + img_dist)

        fov = seq_model.optical_spec.field_of_view
        for fi, f in enumerate(fov.fields):
            rb = gitm.RayBundle(seq_model, fi, start_offset)
            gscene.addItem(rb)

    def create_lens_layout_view(self):
        fig = LensLayoutFigure(self.app_manager.model)
        view_width = 660
        view_ht = 440
        title = "Lens Layout View"
        self.create_plot_view(fig, title, view_width, view_ht,
                              add_scale_panel=False)

    def update_figure_view(plotFigure):
        plotFigure.refresh()

    def create_paraxial_design_view(self, dgm_type):
        seq_model = self.app_manager.model.seq_model
        fig = pdf.ParaxialDesignFigure(seq_model, self.refresh_gui, dgm_type,
                                       figsize=(5, 4))
        view_width = 500
        view_ht = 500
        title = "Paraxial Design View"
        self.create_plot_view(fig, title, view_width, view_ht,
                              add_scale_panel=False)

    def create_ray_fan_view(self, data_type):
        seq_model = self.app_manager.model.seq_model
        fig = aaf.RayFanFigure(seq_model, data_type,
                               scale_type=aaf.Fit_All_Same,
                               figsize=(5, 4), dpi=100)
        view_width = 600
        view_ht = 600
        if data_type == "Ray":
            title = "Ray Fan View"
        elif data_type == "OPD":
            title = "OPD Fan View"
        else:
            title = "bad data_type argument"
        self.create_plot_view(fig, title, view_width, view_ht)

    def create_ray_grid_view(self):
        seq_model = self.app_manager.model.seq_model
        num_flds = len(seq_model.optical_spec.field_of_view.fields)

        fig = aaf.SpotDiagramFigure(seq_model, scale_type=aaf.Fit_All_Same,
                                    num_rays=11, dpi=100)
        view_box = 300
        view_width = view_box
        view_ht = num_flds * view_box
        title = "Spot Diagram"
        self.create_plot_view(fig, title, view_width, view_ht)

    def create_wavefront_view(self):
        seq_model = self.app_manager.model.seq_model
        num_flds = len(seq_model.optical_spec.field_of_view.fields)
        num_wvls = len(seq_model.optical_spec.spectral_region.wavelengths)

        fig = aaf.WavefrontFigure(seq_model, scale_type=aaf.Fit_All_Same,
                                  num_rays=32, dpi=100)
#                                  figsize=(2*num_wvls, 2*num_flds))
        view_box = 300
        view_width = num_wvls * view_box
        view_ht = num_flds * view_box
        title = "Wavefront Map"
        self.create_plot_view(fig, title, view_width, view_ht)

    def create_plot_view(self, fig, title, view_width, view_ht,
                         add_scale_panel=True):
        pc = plotter.PlotCanvas(self, fig)
        # construct the top level widget
        widget = QWidget()
        # construct the top level layout
        layout = QVBoxLayout(widget)

        # set the layout on the widget
        widget.setLayout(layout)

        if add_scale_panel:
            psp = self.create_plot_scale_panel(pc)
            layout.addWidget(psp)

        sub = self.add_subwindow(widget, ModelInfo(self.app_manager.model,
                                         MainWindow.update_figure_view, fig))
        sub.setWindowTitle(title)
        orig_x, orig_y = self.initial_window_offset()
        sub.setGeometry(orig_x, orig_y, view_width, view_ht)

        layout.addWidget(pc)

        sub.show()

    def create_plot_scale_panel(self, pc):
        groupBox = QGroupBox("Plot Scale", self)

        user_scale_wdgt = QLineEdit()
        user_scale_wdgt.setReadOnly(True)
        pf = pc.figure
        cntxt = pf, user_scale_wdgt
        user_scale_wdgt.editingFinished.connect(lambda: self.
                                                on_plot_scale_changed(cntxt))
        fit_all_btn = QRadioButton("Fit All")
        fit_all_btn.setChecked(pf.scale_type == aaf.Fit_All)
        fit_all_btn.toggled.connect(lambda: self.
                                    on_plot_scale_toggled(cntxt, aaf.Fit_All))
        fit_all_same_btn = QRadioButton("Fit All Same")
        fit_all_same_btn.setChecked(pf.scale_type == aaf.Fit_All_Same)
        fit_all_same_btn.toggled.connect(lambda: self.on_plot_scale_toggled(
                                                 cntxt, aaf.Fit_All_Same))
        user_scale_btn = QRadioButton("User Scale")
        user_scale_btn.setChecked(pf.scale_type == aaf.User_Scale)
        user_scale_btn.toggled.connect(lambda: self.on_plot_scale_toggled(
                                       cntxt, aaf.User_Scale))
        box = QHBoxLayout()
        box.addWidget(fit_all_btn)
        box.addWidget(fit_all_same_btn)
        box.addWidget(user_scale_btn)
        box.addWidget(user_scale_wdgt)

        groupBox.setLayout(box)

        return groupBox

    def on_plot_scale_toggled(self, cntxt, scale_type):
        plotFigure, scale_wdgt = cntxt
        plotFigure.scale_type = scale_type
        if scale_type == aaf.User_Scale:
            scale_wdgt.setReadOnly(False)
            scale_wdgt.setText('{:7.4f}'.format(plotFigure.user_scale_value))
        else:
            scale_wdgt.setReadOnly(True)

        plotFigure.plot()

    def on_plot_scale_changed(self, cntxt):
        plotFigure, scale_wdgt = cntxt
        eval_str = scale_wdgt.text()
        try:
            val = eval(eval_str)
            plotFigure.user_scale_value = val
            scale_wdgt.setText('{:7.4f}'.format(val))
        except IndexError:
            return ''

        plotFigure.plot()

    def create_table_view(self, table_model, table_title):
        # construct the top level widget
        widget = QWidget()
        # construct the top level layout
        layout = QVBoxLayout(widget)

        tableView = QTableView()
        tableView.setAlternatingRowColors(True)

        # Add table to box layout
        layout.addWidget(tableView)

        # set the layout on the widget
        widget.setLayout(layout)

        sub = self.add_subwindow(widget, ModelInfo(self.app_manager.model,
                                         MainWindow.update_table_view,
                                         tableView))
        sub.setWindowTitle(table_title)

        tableView.setModel(table_model)

        tableView.setMinimumWidth(tableView.horizontalHeader().length() +
                                  tableView.horizontalHeader().height())
#                                  The following line should work but returns 0
#                                  tableView.verticalHeader().width())

        view_width = tableView.width()
        view_ht = tableView.height()
        orig_x, orig_y = self.initial_window_offset()
        sub.setGeometry(orig_x, orig_y, view_width, view_ht)

        # table data updated successfully
        table_model.update.connect(self.on_data_changed)

        sub.show()

    def update_table_view(table_view):
        table_model = table_view.model()
        table_model.endResetModel()

    def refresh_gui(self):
        self.app_manager.refresh_gui()

    @pyqtSlot(object, int)
    def on_data_changed(self, rootObj, index):
        self.refresh_gui()

    @pyqtSlot(QMdiSubWindow)
    def on_subwindow_activated(self, window):
        self.app_manager.on_window_activated(window)


def main():
    app = QApplication(sys.argv)
    logging.basicConfig(filename='rayoptics.log',
                        filemode='w',
                        level=logging.INFO)
    ex = MainWindow()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
