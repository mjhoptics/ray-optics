#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 12:50:03 2018

@author: Mike
"""

import sys

from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QHBoxLayout,
                             QVBoxLayout, QSizePolicy, QGroupBox, QCheckBox,
                             QTableWidget, QTableWidgetItem)
from PyQt5.QtCore import Qt as qt

from matplotlib.backends.backend_qt5agg \
     import (FigureCanvasQTAgg as FigureCanvas,
             NavigationToolbar2QT as NavigationToolbar)

from matplotlib.figure import Figure
# import matplotlib.pyplot as plt

# import numpy as np

import glass.hoya as h
import glass.ohara as o
import glass.schott as s

Hoya, Ohara, Schott = range(3)
pickTableHeader = ["Catalog", "Glass", "Nd", "Vd"]


class GlassMapViewer(QMainWindow):
    def __init__(self):
        super().__init__()
        self.title = 'Glass Map Viewer'
        self.left = 50
        self.top = 150
        self.width = 1100
        self.height = 650
        self.dataSets = GlassMapModel()
        self.displayDataSets = [True, True, True]
        self.initUI()

    def initUI(self):
        self._main = QWidget()
        self.setCentralWidget(self._main)
        layout = QHBoxLayout(self._main)

        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        self.gmt = self.createTable()
        self.gm = PlotCanvas(self, width=5, height=4)
        self.addToolBar(qt.BottomToolBarArea,
                        NavigationToolbar(self.gm, self))
        layout.addWidget(self.gm)

        rightBar = QVBoxLayout()
        layout.addLayout(rightBar)

        catalogGroup = self.createCatalogGroupBox()
        rightBar.addWidget(catalogGroup)

        rightBar.addWidget(self.gmt)

    def createCatalogGroupBox(self):
        groupBox = QGroupBox("Glass Catalogs", self)

        hoya_checkBox = QCheckBox("&Hoya")
        hoya_checkBox.setChecked(True)
        hoya_checkBox.stateChanged.connect(self.hoya_check)
        ohara_checkBox = QCheckBox("&Ohara")
        ohara_checkBox.setChecked(True)
        ohara_checkBox.stateChanged.connect(self.ohara_check)
        schott_checkBox = QCheckBox("&Schott")
        schott_checkBox.setChecked(True)
        schott_checkBox.stateChanged.connect(self.schott_check)

        vbox = QVBoxLayout()
        vbox.addWidget(hoya_checkBox)
        vbox.addWidget(ohara_checkBox)
        vbox.addWidget(schott_checkBox)

        groupBox.setLayout(vbox)

        return groupBox

    def createTable(self):
        table = PickTable()
        return table

    def hoya_check(self, state):
        checked = state == qt.Checked
        self.gm.displayDataSets[Hoya] = checked
        self.gm.updateVisibility(Hoya, checked)

    def ohara_check(self, state):
        checked = state == qt.Checked
        self.gm.displayDataSets[Ohara] = checked
        self.gm.updateVisibility(Ohara, checked)

    def schott_check(self, state):
        checked = state == qt.Checked
        self.gm.displayDataSets[Schott] = checked
        self.gm.updateVisibility(Schott, checked)


class GlassMapModel():

    def __init__(self):
        self.dataSetList = []
        self.dataSetList.append((h.HoyaCatalog(), 'Hoya'))
        self.dataSetList.append((o.OharaCatalog(), 'Ohara'))
        self.dataSetList.append((s.SchottCatalog(), 'Schott'))
        self.hoyaCat = h.HoyaCatalog()

    def get_data_at(self, i):
        return self.dataSetList[i][0].glass_map_data()

    def get_data_set_label_at(self, i):
        return self.dataSetList[i][1]


class PickTable(QTableWidget):
    def __init__(self):
        super().__init__(16, 4)
        self.rowFill = 16
        self.setHorizontalHeaderLabels(pickTableHeader)
        for i, w in enumerate([52, 96, 61, 48]):
            self.setColumnWidth(i, w)
        self.setAlternatingRowColors(True)
        self.setMinimumWidth(225)
        self.setMaximumWidth(275)

    def setRowCount(self, count):
        if count > self.rowFill:
            super().setRowCount(count)

    def resizeEvent(self, event):
        super(QTableWidget, self).resizeEvent(event)
        rowsToFill = int(self.height()/self.rowHeight(0))
        if rowsToFill > self.rowFill:
            self.rowFill = rowsToFill
            if rowsToFill > self.rowCount():
                super().setRowCount(rowsToFill)


class PlotCanvas(FigureCanvas):
    dsc = [(113/255, 113/255, 198/255),  # sgi slateblue
           (102/255, 205/255, 0),  # chartreuse 3
           (255/255, 114/255, 86/255)]  # coral 1

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        # fig = Figure(figsize=(width, height), dpi=dpi)
        fig = Figure()
        self.axes = fig.add_subplot(111)
        self.rawData = []
        self.data = parent.dataSets
        self.displayDataSets = parent.displayDataSets
        self.pickTable = parent.gmt
        self.needsClear = True
        self.pickRow = 0

        for i, display in enumerate(self.displayDataSets):
            if display:
                x, y, lbl = self.data.get_data_at(i)
                dsLabel = self.data.get_data_set_label_at(i)
                self.rawData.append([dsLabel, (x, y, lbl)])

        super().__init__(fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)

        FigureCanvas.updateGeometry(self)
        self.plot()

    def plot(self):
        self.axes.set_title('Glass Map')
        for i, display in enumerate(self.displayDataSets):
            if display:
                self.axes.plot(self.rawData[i][1][0], self.rawData[i][1][1],
                               linestyle='None', marker='o', markersize=5,
                               color=self.dsc[i], alpha=0.5, gid=i,
                               label=self.rawData[i][0], picker=5)

        self.figure.canvas.mpl_connect('pick_event', self.on_pick)
        self.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.axes.invert_xaxis()
        self.axes.grid()
        self.axes.legend()
        self.draw()

    def clearPickTable(self):
        self.pickTable.clearContents()
        self.pickRow = 0
        self.needsClear = False

    def on_press(self, event):
        if self.needsClear:
            # If needsClear is still set, there have been no pick events so
            #  this is a click in an empty region of the plot.
            #  Clear the pick table
            self.clearPickTable()
        else:
            # on_press event happens after on_pick events. Set needsClear for
            #  next on_pick, i.e. a new selection, to handle
            self.needsClear = True

    def on_pick(self, event):
        if self.needsClear:
            self.clearPickTable()
        line = event.artist
        id = line.get_gid()
        if self.displayDataSets[id]:
            ind = event.ind
            dsLabel = self.rawData[id][0]
            x, y, lbl = self.rawData[id][1]
            self.pickTable.setRowCount(self.pickRow+len(ind))
            for k in ind:
                glass = (dsLabel, lbl[k], y[k], x[k])
                for j in range(4):
                    item = QTableWidgetItem(str(glass[j]))
                    self.pickTable.setItem(self.pickRow, j, item)
                self.pickRow += 1

    def updateVisibility(self, indx, state):
        self.axes.lines[indx].set_visible(state)
        self.draw()


if __name__ == '__main__':
    qapp = QApplication(sys.argv)
    app = GlassMapViewer()
    app.show()
    qapp.exec_()
