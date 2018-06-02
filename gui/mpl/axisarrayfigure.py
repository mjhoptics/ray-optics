#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
"""
Created on Fri Apr  13 10:43:19 2018

@author: Michael J. Hayford
"""

from matplotlib.figure import Figure
from scipy.interpolate import spline

import numpy as np

Fit_All, Fit_All_Same, User_Scale = range(3)


def clip_to_range(rgb_list, lower, upper):
    rgbc_list = []
    for rgb in rgb_list:
        rgbc = []
        for i, clr in enumerate(rgb):
            rgbc.append(clr)
            if clr < lower:
                rgbc[i] = lower
            if upper < clr:
                rgbc[i] = upper
        rgbc_list.append(rgbc)
    return rgbc


class AxisArrayFigure(Figure):

    def __init__(self, seq_model,
                 scale_type=Fit_All,
                 user_scale_value=0.1, **kwargs):
        self.seq_model = seq_model
        self.user_scale_value = user_scale_value
        self.scale_type = scale_type
        self.max_value_all = 0.0

        Figure.__init__(self, **kwargs)

        self.num_rows = len(self.seq_model.optical_spec.field_of_view.fields)
        self.num_cols = 2
        self.update_data()

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
                ax = self.add_subplot(m, n, k)
                self.init_axis(ax)
#                title = "["+str(i)+"],["+str(j)+"]"
#                print("title, id(ax):", title, id(ax))
#                ax.set_title(title)
                row.append(ax)
                k += 1
            arr.append(row)
        return arr

    def update_data(self):
        self.axis_data_array = []
        for i in reversed(range(self.num_rows)):
            row = []
            for j in reversed(range(self.num_cols)):
                x_smooth = []
                y_smooth = []
                x_data, y_data, max_value, rc = self.eval_axis_data(i, j)
#                rc = clip_to_range(rc, 0.0, 1.0)
                for k in range(len(x_data)):
                    x_smooth.append(np.linspace(x_data[k].min(),
                                                x_data[k].max(), 100))
                    y_smooth.append(spline(x_data[k], y_data[k], x_smooth[k]))
                row.append((x_smooth, y_smooth, max_value, rc))
            self.axis_data_array.append(row)

    def plot(self):
        self.clf()

        m = self.num_rows - 1
        n = self.num_cols - 1
        self.ax_arr = self.construct_plot_array(self.num_rows, self.num_cols)

        self.max_value_all = 0.0
        for i in reversed(range(self.num_rows)):
            for j in reversed(range(self.num_cols)):
                x_data, y_data, max_value, rc = self.axis_data_array[m-i][n-j]
                ax = self.ax_arr[m-i][n-j]
                for k in range(len(x_data)):
                    ax.plot(x_data[k], y_data[k], c=rc[k])

                if max_value > self.max_value_all:
                    self.max_value_all = max_value

        if self.scale_type == Fit_All:
            pass
#            print("Fit_All", self.max_value_all)
#            [[ax.set_ylim(-mv, mv) for ax in r] for r in self.ax_arr]
        if self.scale_type == Fit_All_Same:
            mv = self.max_value_all
#            print("Fit_All_Same", mv)
            [[ax.set_ylim(-mv, mv) for ax in r] for r in self.ax_arr]
        if self.scale_type == User_Scale and self.user_scale_value is not None:
            us = self.user_scale_value
#            print("User_Scale", us)
            [[ax.set_ylim(-us, us) for ax in r] for r in self.ax_arr]

        self.canvas.draw()

        return self
