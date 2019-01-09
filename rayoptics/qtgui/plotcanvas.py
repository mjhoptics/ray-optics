#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
"""

.. Created on Fri Mar  2 09:36:19 2018

.. codeauthor: Michael J. Hayford
"""

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QSizePolicy

from matplotlib.backends.backend_qt5agg \
     import FigureCanvasQTAgg as FigureCanvas


class PlotCanvas(FigureCanvas):

    def __init__(self, parent, figure):
        FigureCanvas.__init__(self, figure)
        self.setParent(parent)

        # Next 2 lines are needed so that key press events are correctly
        #  passed with mouse events
        # https://github.com/matplotlib/matplotlib/issues/707/
        self.setFocusPolicy(Qt.ClickFocus)
        self.setFocus()

        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

        self.figure.plot()
