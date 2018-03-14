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


class PlotCanvas(FigureCanvas):

    def __init__(self, parent, seq_model, width=5, height=4, dpi=100):
        self.seq_model = seq_model

        self.fig = Figure(figsize=(width, height), dpi=dpi)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.plot()

    def plot(self):
        num_plot = 1
        num_flds = len(self.seq_model.optical_spec.field_of_view.fields)
        for fi in reversed(range(num_flds)):
            for xy in reversed(range(2)):
                fans_x, fans_y = self.seq_model.trace_fan(fi, xy)
                ax = self.fig.add_subplot(num_flds, 2, num_plot)
                ax.plot(fans_x, fans_y)
                ax.grid(True)
                ax.set_xlim(-1., 1.)
                num_plot += 1

#        self.axes.set_title('Ray Fan Plot')
        self.draw()
