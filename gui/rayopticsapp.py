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
from PyQt5.QtGui import (QColor, QStandardItemModel, QStandardItem)
from PyQt5.QtWidgets import (QApplication, QAction, QMainWindow, QMdiArea,
                             QMdiSubWindow, QTextEdit, QFileDialog, QTableView,
                             QVBoxLayout, QWidget, QGraphicsView,
                             QGraphicsScene)
from PyQt5.QtCore import pyqtSlot

import codev.cmdproc as cvp
import optical.opticalmodel as optm
import optical.sequential as seq
import optical.elements as ele
import gui.plotcanvas as plotter
import gui.pytablemodel as tbl
import gui.graphicsitems as gitm


class MainWindow(QMainWindow):
    CURVATURE, THICKNESS, MATERIAL, SEMIDIAM = range(4)
    count = 0

    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.opt_model = optm.OpticalModel()
        self.mdi = QMdiArea()
        self.setCentralWidget(self.mdi)

        self.left = 100
        self.top = 50
        self.offset_x = 100
        self.offset_y = 25
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

    def file_action(self, q):
        if q.text() == "New":
            MainWindow.count = MainWindow.count+1
            sub = QMdiSubWindow()
            sub.setWidget(QTextEdit())
            sub.setWindowTitle("subwindow"+str(MainWindow.count))
            self.mdi.addSubWindow(sub)
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
        self.cur_filename = file_name
        self.is_changed = True
        cvp.read_lens(self.opt_model, file_name)
        self.create_lens_table()
        self.create_2D_lens_view()

    def view_action(self, q):
        if q.text() == "Lens Table":
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
            self.create_table_view(model, "Lens Table")

        if q.text() == "Lens View":
            self.create_2D_lens_view()

        if q.text() == "Ray Fans":
            self.create_ray_fan_view()

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
        print("window triggered")

        if q.text() == "Cascade":
            self.mdi.cascadeSubWindows()

        if q.text() == "Tiled":
            self.mdi.tileSubWindows()

    def create_lens_table(self):
        # construct the top level widget
        widget = QWidget()
        # construct the top level layout
        layout = QVBoxLayout(widget)

        tableView = QTableView()
        tableView.setAlternatingRowColors(True)
        # table selection change
#        self.tableView.doubleClicked.connect(self.on_click)

        # Add table to box layout
        layout.addWidget(tableView)

        # set the layout on the widget
        widget.setLayout(layout)

        sub = self.mdi.addSubWindow(widget)
        sub.setWindowTitle("Surface Data Table")

        model = self.createSurfaceModel(self)
        tableView.setModel(model)
        seq_model = self.opt_model.seq_model
        for s in range(len(seq_model.surfs)):
            self.addSurface(model, seq_model.list_surface_and_gap(s))

        tableView.setMinimumWidth(tableView.horizontalHeader().length() +
                                  tableView.horizontalHeader().height())
#                                  The following line should work but returns 0
#                                  tableView.verticalHeader().width())
        MainWindow.count += 1
        sub.show()

    def create_2D_lens_view(self):
        self.scene2d = QGraphicsScene()
        self.create_element_model(self.scene2d)
        self.create_ray_model(self.scene2d)
        self.scene2d.setBackgroundBrush(QColor(237, 243, 254))  # light blue
        sceneRect2d = self.scene2d.sceneRect()
#        print("Scene rect1:", sceneRect2d.width()/sceneRect2d.height(),
#              sceneRect2d.x(), sceneRect2d.y(),
#              sceneRect2d.width(), sceneRect2d.height())

        # construct the top level widget
        widget = QWidget()
        # construct the top level layout
        layout = QVBoxLayout(widget)

        # set the layout on the widget
        widget.setLayout(layout)

        sub = self.mdi.addSubWindow(widget)
        sub.setWindowTitle("2D Lens View")
        view_width = 600
        view_ht = 400
        view_ratio = view_width/view_ht
        orig_x = MainWindow.count*self.offset_x
        orig_y = MainWindow.count*self.offset_y
        sub.setGeometry(orig_x, orig_y, view_width, view_ht)

        self.gview2d = QGraphicsView(self.scene2d)
#        self.gview2d.setGeometry(100, 50, view_width, view_ht)
        scene_ratio = sceneRect2d.width()/sceneRect2d.height()
        oversize_fraction = 1.2
        if scene_ratio > view_ratio:
            view_scale = view_width/(oversize_fraction*sceneRect2d.width())
        else:
            view_scale = view_ht/(oversize_fraction*sceneRect2d.height())

#        print(view_ratio, scene_ratio, view_scale)
        frame_before = self.gview2d.frameGeometry()
#        print("Frame before:", frame_before.x(), frame_before.y(),
#              frame_before.width(), frame_before.height())
        self.gview2d.scale(view_scale, view_scale)
        layout.addWidget(self.gview2d)

        MainWindow.count += 1
        sub.show()

    def update_2D_lens_view(self):
        for gi in self.scene2d.items():
            gi.prepareGeometryChange()
            gi.update_shape()

    def create_ray_fan_view(self):
        # construct the top level widget
        widget = QWidget()
        # construct the top level layout
        layout = QVBoxLayout(widget)

        # set the layout on the widget
        widget.setLayout(layout)

        sub = self.mdi.addSubWindow(widget)
        sub.setWindowTitle("Ray Fan View")
        view_width = 600
        view_ht = 600
        orig_x = MainWindow.count*self.offset_x
        orig_y = MainWindow.count*self.offset_y
        sub.setGeometry(orig_x, orig_y, view_width, view_ht)

        seq_model = self.opt_model.seq_model
        pc = plotter.PlotCanvas(self, seq_model, width=5, height=4)
        layout.addWidget(pc)

        MainWindow.count += 1
        sub.show()

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

        sub = self.mdi.addSubWindow(widget)
        sub.setWindowTitle(table_title)

        tableView.setModel(table_model)

        tableView.setMinimumWidth(tableView.horizontalHeader().length() +
                                  tableView.horizontalHeader().height())
#                                  The following line should work but returns 0
#                                  tableView.verticalHeader().width())

        # table data updated successfully
        table_model.update.connect(self.on_data_changed)

        MainWindow.count += 1
        sub.show()

    @pyqtSlot(object, int)
    def on_data_changed(self, rootObj, index):
#        print("on_data_changed - index:", index)
        self.opt_model.update_model()
        self.update_2D_lens_view()

    def createSurfaceModel(self, parent):
        model = QStandardItemModel(0, 4, parent)
        model.setHeaderData(self.CURVATURE, Qt.Horizontal, "Curvature")
        model.setHeaderData(self.THICKNESS, Qt.Horizontal, "Thickness")
        model.setHeaderData(self.MATERIAL, Qt.Horizontal, "Material")
        model.setHeaderData(self.SEMIDIAM, Qt.Horizontal, "Semi-Diameter")
        return model

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

    def addSurface(self, model, surf_gap):
        itemList = []
        for i in range(4):
            item = QStandardItem()
            item.setData(surf_gap[i], Qt.DisplayRole)
            itemList.append(item)
        model.appendRow(itemList)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    logging.basicConfig(filename='rayoptics.log',
                        filemode='w',
                        level=logging.INFO)
    ex = MainWindow()
    sys.exit(app.exec())
