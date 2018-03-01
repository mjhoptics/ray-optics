#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 09:24:01 2018

@author: Mike
"""

import sys
import logging

from PyQt5.QtCore import Qt, QPointF, QRectF
from PyQt5.QtGui import (QBrush, QPen, QColor, QStandardItemModel,
                         QStandardItem, QPolygonF)
from PyQt5.QtWidgets import (QApplication, QAction, QMainWindow, QMdiArea,
                             QMdiSubWindow, QTextEdit, QFileDialog, QTableView,
                             QVBoxLayout, QWidget, QGraphicsView,
                             QGraphicsScene, QGraphicsPolygonItem,
                             QGraphicsRectItem)

import codev.cmdproc as cvp
import optical.sequential as seq
import optical.elements as ele
import gui.rgbtable as rgbt
import numpy as np


class MainWindow(QMainWindow):
    CURVATURE, THICKNESS, MATERIAL, SEMIDIAM = range(4)
    count = 0

    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.seq_model = seq.SequentialModel()
        self.mdi = QMdiArea()
        self.setCentralWidget(self.mdi)

        bar = self.menuBar()

        file = bar.addMenu("File")
        file.addAction("New")
        file.addAction("Open...")
        file.addSeparator()
        file.addAction("Save")
        file.addAction("Save As...")
        file.triggered[QAction].connect(self.file_action)
        view = bar.addMenu("View")
        view.addAction("Table")
        view.addAction("Lens View")
        view.triggered[QAction].connect(self.view_action)
        wnd = bar.addMenu("Window")
        wnd.addAction("Cascade")
        wnd.addAction("Tiled")
        wnd.triggered[QAction].connect(self.window_action)
        self.setWindowTitle("Ray Optics")
        self.show()
        self.open_file("/Users/Mike/Developer/PyProjects/ray-optics/codev/test/ag_dblgauss.seq")

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
        cvp.read_lens(self.seq_model, file_name)
        self.create_lens_table()
        self.create_2D_lens_view()

    def view_action(self, q):
        print("view triggered")

        if q.text() == "Table":
            MainWindow.count = MainWindow.count+1
            self.create_lens_table()

        if q.text() == "Lens View":
            MainWindow.count = MainWindow.count+1
            self.create_2D_lens_view()

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
        for s in range(len(self.seq_model.surfs)):
            self.addSurface(model, self.seq_model.list_surface_and_gap(s))

        tableView.setMinimumWidth(tableView.horizontalHeader().length() +
                                  tableView.horizontalHeader().height())
#                                  The following line should work but returns 0
#                                  tableView.verticalHeader().width())
        sub.show()

    def create_2D_lens_view(self):
        self.scene2d = QGraphicsScene()
        self.create_element_model(self.scene2d)
        self.create_ray_model(self.scene2d)
        self.scene2d.setBackgroundBrush(QColor(237, 243, 254))  # light blue
        sceneRect2d = self.scene2d.sceneRect()
        print("Scene rect1:", sceneRect2d.width()/sceneRect2d.height(),
              sceneRect2d.x(), sceneRect2d.y(),
              sceneRect2d.width(), sceneRect2d.height())

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
        sub.setGeometry(100, 50, view_width, view_ht)

        self.gview2d = QGraphicsView(self.scene2d)
#        self.gview2d.setGeometry(100, 50, view_width, view_ht)
        scene_ratio = sceneRect2d.width()/sceneRect2d.height()
        oversize_fraction = 1.2
        if scene_ratio > view_ratio:
            view_scale = view_width/(oversize_fraction*sceneRect2d.width())
        else:
            view_scale = view_ht/(oversize_fraction*sceneRect2d.height())

        print(view_ratio, scene_ratio, view_scale)
        frame_before = self.gview2d.frameGeometry()
        print("Frame before:", frame_before.x(), frame_before.y(),
              frame_before.width(), frame_before.height())
        self.gview2d.scale(view_scale, view_scale)
        layout.addWidget(self.gview2d)

        sub.show()

    def createSurfaceModel(self, parent):
        model = QStandardItemModel(0, 4, parent)
        model.setHeaderData(self.CURVATURE, Qt.Horizontal, "Curvature")
        model.setHeaderData(self.THICKNESS, Qt.Horizontal, "Thickness")
        model.setHeaderData(self.MATERIAL, Qt.Horizontal, "Material")
        model.setHeaderData(self.SEMIDIAM, Qt.Horizontal, "Semi-Diameter")
        return model

    def create_element_model(self, gscene):
        clut = rgbt.RGBTable(filename='gui/red_blue64.csv',
                             data_range=[10.0, 100.])
        ele_model = ele.ElementModel()
        ele_model.elements_from_sequence(self.seq_model)
        pen = QPen()
        pen.setCosmetic(True)
        for e in ele_model.elements:
            poly = e.shape()
            polygon = QPolygonF()
            for p in poly:
                polygon.append(QPointF(p[0], p[1]))
            gpoly = QGraphicsPolygonItem()
            gpoly.setPolygon(polygon)
            # set element color based on V-number
            gc = float(e.medium.glass_code())
            vnbr = round(100.0*(gc - int(gc)), 3)
            ergb = clut.get_color(vnbr)
            gpoly.setBrush(QColor(*ergb))
            gpoly.setPen(pen)

            t = e.tfrm[1]
            gpoly.setPos(QPointF(t[2], -t[1]))
            gscene.addItem(gpoly)

    def create_ray_model(self, gscene, start_surf=1):
        tfrms = self.seq_model.compute_global_coords(start_surf)
        rayset = self.seq_model.trace_boundary_rays()

        start_offset = 0.1*gscene.sceneRect().width()
        if abs(tfrms[0][1][2]) > start_offset:
            tfrms[0] = self.seq_model.shift_start_of_rayset(rayset,
                                                            start_offset)

        pen = QPen()
        pen.setCosmetic(True)
        for rays in rayset:
            poly1 = []
            for i, r in enumerate(rays[3][0][0:]):
                rot, trns = tfrms[i]
                p = rot.dot(r[0]) + trns
#                print(i, r[0], rot, trns, p)
                poly1.append(QPointF(p[2], -p[1]))

            poly2 = []
            for i, r in enumerate(rays[4][0][0:]):
                rot, trns = tfrms[i]
                p = rot.dot(r[0]) + trns
#                print(i, r[0], rot, trns, p)
                poly2.append(QPointF(p[2], -p[1]))

            poly2.reverse()
            poly1.extend(poly2)
            polygon = QPolygonF()
            for p in poly1:
                polygon.append(p)
            gpoly = QGraphicsPolygonItem()
            gpoly.setPolygon(polygon)
            gpoly.setBrush(QColor(254, 197, 254, 64))  # magenta, 25%
            gpoly.setPen(pen)
            gscene.addItem(gpoly)

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
