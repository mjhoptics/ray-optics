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

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import (QApplication, QAction, QMainWindow, QMdiArea,
                             QMdiSubWindow, QTextEdit, QFileDialog, QTableView,
                             QHBoxLayout, QVBoxLayout, QWidget, QGraphicsView,
                             QLineEdit, QGraphicsScene, QRadioButton,
                             QGroupBox)
from PyQt5.QtCore import pyqtSlot

import codev.cmdproc as cvp
import optical.opticalmodel as optm
import optical.sequential as seq
import optical.elements as ele
import gui.plotcanvas as plotter
import gui.mpl.axisarrayfigure as aaf
import gui.mpl.paraxdgnfigure as pdf
import gui.pytablemodel as tbl
import gui.graphicsitems as gitm


class MainWindow(QMainWindow):
    count = 0

    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.opt_model = None
        self.mdi = QMdiArea()
        self.setCentralWidget(self.mdi)

        self.wnd_dict = {}
        self.mdi.subWindowActivated.connect(self.on_subwindow_activated)

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
        view.addAction("Ray Fans")
        view.addAction("Paraxial Design View")
        view.addAction("Paraxial Ray Table")
        view.addAction("Ray Table")
        view.addSeparator()
        view.triggered[QAction].connect(self.view_action)
        wnd = bar.addMenu("Window")
        wnd.addAction("Cascade")
        wnd.addAction("Tiled")
        wnd.triggered[QAction].connect(self.window_action)

        self.setWindowTitle("Ray Optics")
        self.show()
#        self.open_file("/Users/Mike/Developer/PyProjects/ray-optics/"
#                       "codev/test/asp46.seq")
#        self.open_file("/Users/Mike/Developer/PyProjects/ray-optics/"
#                       "codev/test/paraboloid.seq")
#        self.open_file("/Users/Mike/Developer/PyProjects/ray-optics/"
#                       "codev/test/schmidt.seq")
#        self.open_file("/Users/Mike/Developer/PyProjects/ray-optics/"
#                       "codev/test/schmidt_sph.seq")
        self.open_file("/Users/Mike/Developer/PyProjects/ray-optics/"
                       "codev/test/ag_dblgauss.seq")

    def add_subwindow(self, widget, model_info):
            sub_wind = self.mdi.addSubWindow(widget)
            self.wnd_dict[sub_wind] = model_info
            MainWindow.count += 1
            return sub_wind

    def delete_subwindow(self, sub_wind):
            self.mdi.removeSubWindow(sub_wind)
            del self.wnd_dict[sub_wind]
            MainWindow.count -= 1

    def initial_window_offset(self):
        offset_x = 50
        offset_y = 25
        orig_x = (MainWindow.count - 1)*offset_x
        orig_y = (MainWindow.count - 1)*offset_y
        return orig_x, orig_y

    def file_action(self, q):
        if q.text() == "New":
            self.opt_model = optm.OpticalModel()
            sub = self.add_subwindow(QTextEdit(),
                                     (self.opt_model, ))
            sub.setWindowTitle("subwindow"+str(MainWindow.count))
            sub.show()

        if q.text() == "Open...":
            options = QFileDialog.Options()
            # options |= QFileDialog.DontUseNativeDialog
            fileName, _ = QFileDialog.getOpenFileName(self,
                          "QFileDialog.getOpenFileName()",
                          "",
                          "CODE V Files (*.seq)",
                          options=options)
            if fileName:
                logging.debug("open file: %s", fileName)
                self.open_file(fileName)

    def open_file(self, file_name):
        self.opt_model = optm.OpticalModel()
        self.cur_filename = file_name
        self.is_changed = True
        cvp.read_lens(self.opt_model, file_name)
        self.create_lens_table()
        self.create_2D_lens_view()

    def view_action(self, q):
        if q.text() == "Lens Table":
            self.create_lens_table()

        if q.text() == "Lens View":
            self.create_2D_lens_view()

        if q.text() == "Ray Fans":
            self.create_ray_fan_view()

        if q.text() == "Paraxial Design View":
            self.create_paraxial_design_view()

        if q.text() == "Paraxial Ray Table":
            parax_data = self.opt_model.seq_model.optical_spec.parax_data
            colEvalStr = ['[0][{}][0]', '[0][{}][1]', '[0][{}][2]',
                          '[1][{}][0]', '[1][{}][1]', '[1][{}][2]']
            rowHeaders = self.opt_model.seq_model.surface_label_list()
            colHeaders = ['y', 'u', 'i', 'y-bar', 'u-bar', 'i-bar']
            colFormats = ['{:12.5g}', '{:9.6f}', '{:9.6f}', '{:12.5g}',
                          '{:9.6f}', '{:9.6f}']
            model = tbl.PyTableModel(parax_data, colEvalStr, rowHeaders,
                                     colHeaders, colFormats, False)
            self.create_table_view(model, "Paraxial Ray Table")

        if q.text() == "Ray Table":
            seq_model = self.opt_model.seq_model
            r2f1, _ = seq_model.optical_spec.trace(seq_model, [0., 1.], 0)
            colEvalStr = ['[{}][0][0]', '[{}][0][1]', '[{}][0][2]',
                          '[{}][1][0]', '[{}][1][1]', '[{}][1][2]',
                          '[{}][2]']
            rowHeaders = seq_model.surface_label_list()
            colHeaders = ['x', 'y', 'z', 'l', 'm', 'n', 'length']
            colFormats = ['{:12.5g}', '{:12.5g}', '{:12.5g}', '{:9.6f}',
                          '{:9.6f}', '{:9.6f}', '{:12.5g}']
            model = tbl.PyTableModel(r2f1, colEvalStr, rowHeaders,
                                     colHeaders, colFormats, False)
            self.create_table_view(model, "Ray Table")

    def window_action(self, q):
        if q.text() == "Cascade":
            self.mdi.cascadeSubWindows()

        if q.text() == "Tiled":
            self.mdi.tileSubWindows()

    def create_lens_table(self):
        seq_model = self.opt_model.seq_model
        colEvalStr = ['.surfs[{}].profile.type', '.surfs[{}].profile.cv',
                      '.surfs[{}].surface_od()', '.gaps[{}].thi',
                      '.gaps[{}].medium.name()', '.surfs[{}].refract_mode']
        rowHeaders = seq_model.surface_label_list()
        colHeaders = ['type', 'cv', 'sd', 'thi', 'medium', 'mode']
        colFormats = ['{:s}', '{:12.7g}', '{:12.5g}', '{:12.5g}',
                      '{:s}', '{:s}']
        model = tbl.PyTableModel(seq_model, colEvalStr, rowHeaders,
                                 colHeaders, colFormats, True)
        self.create_table_view(model, "Surface Data Table")

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

        sub = self.add_subwindow(widget, (self.opt_model,
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
        ele_model = self.opt_model.ele_model
        ele_model.elements_from_sequence(self.opt_model.seq_model)
        for e in ele_model.elements:
            ge = gitm.OpticalElement(e)
            gscene.addItem(ge)

    def create_ray_model(self, gscene, start_surf=1):
        seq_model = self.opt_model.seq_model

        img_dist = abs(seq_model.optical_spec.parax_data[2].img_dist)
        start_offset = 0.05*(gscene.sceneRect().width() + img_dist)

        fov = seq_model.optical_spec.field_of_view
        for fi, f in enumerate(fov.fields):
            rb = gitm.RayBundle(seq_model, fi, start_offset)
            gscene.addItem(rb)

    def create_paraxial_design_view(self):
        seq_model = self.opt_model.seq_model
        fig = pdf.ParaxialDesignFigure(seq_model, figsize=(5, 4))
        pc = plotter.PlotCanvas(self, fig)
        # construct the top level widget
        widget = QWidget()
        # construct the top level layout
        layout = QVBoxLayout(widget)

        # set the layout on the widget
        widget.setLayout(layout)

        sub = self.add_subwindow(widget,
                                 (self.opt_model,
                                  MainWindow.update_paraxial_design_view, fig))
        sub.setWindowTitle("Paraxial Design View")
        view_width = 500
        view_ht = 500
        orig_x, orig_y = self.initial_window_offset()
        sub.setGeometry(orig_x, orig_y, view_width, view_ht)

        layout.addWidget(pc)

        sub.show()

    def update_paraxial_design_view(plotFigure):
        plotFigure.update_data()
        plotFigure.plot()

    def create_ray_fan_view(self):
        seq_model = self.opt_model.seq_model
        fig = aaf.AxisArrayFigure(seq_model, aaf.Fit_All_Same,
                                  figsize=(5, 4), dpi=100)
        pc = plotter.PlotCanvas(self, fig)
        psp = self.create_plot_scale_panel(pc)
        # construct the top level widget
        widget = QWidget()
        # construct the top level layout
        layout = QVBoxLayout(widget)

        # set the layout on the widget
        widget.setLayout(layout)

        layout.addWidget(psp)

        sub = self.add_subwindow(widget, (self.opt_model,
                                          MainWindow.update_ray_fan_view, pc))
        sub.setWindowTitle("Ray Fan View")
        view_width = 600
        view_ht = 600
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

    def update_ray_fan_view(plotFigure):
        plotFigure.update_data()
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

        sub = self.add_subwindow(widget, (self.opt_model, ))
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

    @pyqtSlot(object, int)
    def on_data_changed(self, rootObj, index):
        self.opt_model.update_model()
        for model_info in self.wnd_dict.values():
            if model_info[0] == self.opt_model:
                num_items = len(model_info)
                if num_items == 2:
                    model_info[1]()
                elif num_items == 3:
                    model_info[1](model_info[2])

    @pyqtSlot(QMdiSubWindow)
    def on_subwindow_activated(self, window):
        if window is not None:
            try:
                model_info = self.wnd_dict[window]
            except KeyError:
                logging.debug('Window "%s" not in wnd_dict',
                              window.windowTitle())
            else:
                opt_model = model_info[0]
                logging.debug("on_subwindow_activated: %s, %s" %
                              (opt_model.system_spec.title,
                               window.windowTitle()))
                if opt_model and opt_model != self.opt_model:
                    self.opt_model = opt_model
                    logging.debug("switch opt_model to",
                                  opt_model.system_spec.title)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    logging.basicConfig(filename='rayoptics.log',
                        filemode='w',
                        level=logging.INFO)
    ex = MainWindow()
    sys.exit(app.exec())
