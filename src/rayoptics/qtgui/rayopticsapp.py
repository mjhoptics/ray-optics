#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Ray Optics GUI Application

Relies on PyQt5

.. Created on Mon Feb 12 09:24:01 2018

.. codeauthor: Michael J. Hayford
"""

import sys
import logging
from pathlib import Path

from PyQt5.QtCore import Qt as qt
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import (QApplication, QAction, QMainWindow, QMdiArea,
                             QMdiSubWindow, QFileDialog, QTableView, QWidget,
                             QVBoxLayout, QGraphicsView, QGraphicsScene)
from PyQt5.QtCore import pyqtSlot
from traitlets.config.configurable import MultipleInstanceError

from rayoptics.optical.opticalmodel import open_model
import rayoptics.gui.appcmds as cmds
from rayoptics.gui.appmanager import ModelInfo, AppManager
from rayoptics.mpl.paraxdgnfigure import Dgm
import rayoptics.qtgui.dockpanels as dock
from rayoptics.qtgui.graphicsitems import OpticalElement, RayBundle
from rayoptics.qtgui.ipyconsole import create_ipython_console
from rayoptics.optical import trace as trace


class MainWindow(QMainWindow):
    count = 0

    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.mdi = QMdiArea()
        self.setCentralWidget(self.mdi)

        self.app_manager = AppManager(None, gui_parent=self)
        self.mdi.subWindowActivated.connect(self.app_manager.
                                            on_view_activated)

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
        file.addAction("Close")
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
        view.addAction("Astigmatism Curves")
        view.addSeparator()
        view.triggered[QAction].connect(self.view_action)

        wnd = bar.addMenu("Window")
        wnd.addAction("Cascade")
        wnd.addAction("Tiled")
        wnd.addSeparator()

        dock.create_dock_windows(self)
        for pi in dock.panels.values():
            wnd.addAction(pi.menu_action)

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
#            self.open_file(path / "codev/tests/asp46.seq")
#            self.open_file(path / "codev/tests/paraboloid.seq")
#            self.open_file(path / "codev/tests/paraboloid_f8.seq")
#            self.open_file(path / "codev/tests/schmidt.seq")
#            self.open_file(path / "codev/tests/questar35.seq")
#            self.open_file(path / "codev/tests/conic_mirror.seq")
#            self.open_file(path / "codev/tests/rc_f16.seq")
            self.open_file(path / "codev/tests/ag_dblgauss.seq")
#            self.open_file(path / "codev/tests/landscape_lens.seq")
#            self.open_file(path / "optical/tests/cell_phone_camera.roa")
#            self.open_file(path / "optical/tests/singlet_f3.roa")

#        try:
#            root_pos = pth.parts.index('ray-optics')
#        except ValueError:
#            logging.debug("Can't find ray-optics: path is %s", pth)
#        else:
#            path = Path(*pth.parts[:root_pos+1])
#            self.open_file(path / "models/TwoMirror.roa")
#            self.open_file(path / "models/TwoSphericalMirror.roa")
#            self.open_file(path / "models/Sasian Triplet.roa")
#            self.open_file(path / "models/singlet_f5.roa")
#            self.open_file(path / "models/Ritchey_Chretien.roa")
        finally:
            try:
                create_ipython_console(self, 'iPython console', 600, 400)
            except MultipleInstanceError:
                logging.debug("Unable to open iPython console. "
                              "MultipleInstanceError")
            except Exception as inst:
                print(type(inst))    # the exception instance
                print(inst.args)     # arguments stored in .args
                print(inst)          # __str__ allows args to be printed directly,
                pass                 # but may be overridden in exception subclasses

    def add_subwindow(self, widget, model_info):
            sub_wind = self.mdi.addSubWindow(widget)
            self.app_manager.add_view(sub_wind, model_info)
            MainWindow.count += 1
            return sub_wind

    def delete_subwindow(self, sub_wind):
            self.app_manager.delete_view(sub_wind)
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
            opt_model = cmds.create_new_model()
            self.app_manager.model = opt_model
            self.create_lens_table()
#            self.create_2D_lens_view()
            cmds.create_paraxial_design_view(opt_model, Dgm.ht,
                                             gui_parent=self)
            self.refresh_app_ui()

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

        if q.text() == "Close":
            self.close_model()

    def open_file(self, file_name):
        self.cur_filename = file_name
        self.app_manager.model = open_model(file_name)
        self.is_changed = True
        self.create_lens_table()
#        self.create_lens_layout_view()
        self.create_2D_lens_view()
        self.refresh_app_ui()

    def save_file(self, file_name):
        self.app_manager.model.save_model(file_name)
        self.cur_filename = file_name
        self.is_changed = False

    def close_model(self):
        """ NOTE: this does not check to save a modified model """
        self.app_manager.close_model(self.delete_subwindow)

    def view_action(self, q):
        opt_model = self.app_manager.model
        seq_model = self.app_manager.model.seq_model
        if q.text() == "Lens Table":
            self.create_lens_table()

        if q.text() == "Lens View":
            cmds.create_lens_layout_view(self.app_manager.model,
                                         gui_parent=self)
#            self.create_2D_lens_view()

        if q.text() == "Ray Fans":
            cmds.create_ray_fan_view(opt_model, "Ray", gui_parent=self)

        if q.text() == "OPD Fans":
            cmds.create_ray_fan_view(opt_model, "OPD", gui_parent=self)

        if q.text() == "Spot Diagram":
            cmds.create_ray_grid_view(opt_model, gui_parent=self)

        if q.text() == "Wavefront Map":
            cmds.create_wavefront_view(opt_model, gui_parent=self)

        if q.text() == "Astigmatism Curves":
            cmds.create_field_curves(opt_model, gui_parent=self)

        if q.text() == "Paraxial Height View":
            cmds.create_paraxial_design_view(opt_model, Dgm.ht,
                                             gui_parent=self)

        if q.text() == "Paraxial Slope View":
            cmds.create_paraxial_design_view(opt_model, Dgm.slp,
                                             gui_parent=self)

        if q.text() == "Paraxial Ray Table":
            model = cmds.create_parax_table_model(opt_model)
            self.create_table_view(model, "Paraxial Ray Table")

        if q.text() == "Ray Table":
            self.create_ray_table(opt_model)

    def window_action(self, q):
        if q.text() == "Cascade":
            self.mdi.cascadeSubWindows()

        if q.text() == "Tiled":
            self.mdi.tileSubWindows()

    def create_lens_table(self):
        seq_model = self.app_manager.model.seq_model
        model = cmds.create_lens_table_model(seq_model)
        self.create_table_view(model, "Surface Data Table")

    def create_ray_table(self, opt_model):
        osp = opt_model.optical_spec
        pupil = [0., 0.]
        fi = 1
        wl = osp.spectral_region.reference_wvl
        fld, wvl, foc = osp.lookup_fld_wvl_focus(fi, wl)
        ray, ray_op, wvl, opd = trace.trace_with_opd(opt_model, pupil,
                                                     fld, wvl, foc)

        cr = trace.RayPkg(ray, ray_op, wvl)
        s, t = trace.trace_coddington_fan(opt_model, cr, foc)

        model = cmds.create_ray_table_model(opt_model, ray)
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
                                         (scene2d,)))
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
            ge = OpticalElement(e)
            gscene.addItem(ge)

    def create_ray_model(self, gscene, start_surf=1):
        opt_model = self.app_manager.model

        img_dist = abs(opt_model.optical_spec.parax_data[2].img_dist)
        start_offset = 0.05*(gscene.sceneRect().width() + img_dist)

        fov = opt_model.optical_spec.field_of_view
        for fi, f in enumerate(fov.fields):
            rb = RayBundle(opt_model, fi, start_offset)
            gscene.addItem(rb)

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
                                                   cmds.update_table_view,
                                                   (tableView,)))
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

    def refresh_gui(self):
        self.app_manager.refresh_gui()

    def refresh_app_ui(self):
        dock.update_dock_windows(self)

    @pyqtSlot(object, int)
    def on_data_changed(self, rootObj, index):
        self.refresh_gui()


def main():
    qtapp = QApplication(sys.argv)
    logging.basicConfig(filename='rayoptics.log',
                        filemode='w',
                        level=logging.INFO)
    qtwnd = MainWindow()
    return qtapp.exec()


if __name__ == '__main__':
    sys.exit(main())
