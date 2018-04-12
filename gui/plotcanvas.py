#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 09:36:19 2018

@author: Mike
"""

from PyQt5.QtWidgets import QSizePolicy

from matplotlib.backends.backend_qt5agg \
     import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

import numpy as np

Fit_All, Fit_All_Same, User_Scale = range(3)


class PlotCanvas(FigureCanvas):

    def __init__(self, parent, seq_model, width=5, height=4, dpi=100):
        self.seq_model = seq_model
        self.max_value_all = 0.0
        self.user_scale_value = 0.1
        self.scale_type = Fit_All

        self.fig = Figure(figsize=(width, height), dpi=dpi)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

        self.num_rows = len(self.seq_model.optical_spec.field_of_view.fields)
        self.num_cols = 2

        self.update_data()
        self.plot()

    def eval_axis_data(self, i, j):
        """ returns xData, yData, and the max Y value """
        return self.seq_model.trace_fan(i, j)

    def init_axis(self, ax):
        ax.grid(True)
        ax.set_xlim(-1., 1.)
        ax.axvline(0, c='black', lw=1)
        ax.axhline(0, c='black', lw=1)

    def construct_plot_array(self, m, n):
        arr = []
        k = 1
        for i in reversed(range(m)):
            row = []
            for j in reversed(range(n)):
                ax = self.fig.add_subplot(m, n, k)
                self.init_axis(ax)
#                title = "["+str(i)+"],["+str(j)+"]"
#                print("title, id(ax):", title, id(ax))
#                ax.set_title(title)
                row.append(ax)
                k += 1
            arr.append(row)
        return arr

    def update_data(self):
        print("update_data")
        self.axis_data_array = []
        for i in reversed(range(self.num_rows)):
            row = []
            for j in reversed(range(self.num_cols)):
                ax = self.eval_axis_data(i, j)
                row.append(ax)
            self.axis_data_array.append(row)

    def plot(self):
        self.fig.clf()

        m = self.num_rows - 1
        n = self.num_cols - 1
        self.ax_arr = self.construct_plot_array(self.num_rows, self.num_cols)

        self.max_value_all = 0.0
        for i in reversed(range(self.num_rows)):
            for j in reversed(range(self.num_cols)):
                x_data, y_data, max_value = self.axis_data_array[m-i][n-j]
#                x_data, y_data, max_value = self.eval_axis_data(i, j)
                ax = self.ax_arr[m-i][n-j]
                ax.plot(x_data, y_data)
#                print("i, j, ax, title:", i, j, id(ax), ax.get_title())
#                print("i, j, max", i, j, max_value, ax.get_title())
#                print("i, j, xdata", i, j, x_data)
#                print("i, j, ydata", i, j, y_data)

                if max_value > self.max_value_all:
                    self.max_value_all = max_value

        if self.scale_type == Fit_All:
            print("Fit_All", self.max_value_all)
#            [[ax.set_ylim(-mv, mv) for ax in r] for r in self.ax_arr]
        if self.scale_type == Fit_All_Same:
            mv = self.max_value_all
            print("Fit_All_Same", mv)
            [[ax.set_ylim(-mv, mv) for ax in r] for r in self.ax_arr]
        if self.scale_type == User_Scale and self.user_scale_value is not None:
            us = self.user_scale_value
            print("User_Scale", us)
            [[ax.set_ylim(-us, us) for ax in r] for r in self.ax_arr]

        self.draw()
