#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
""" mpl implementations of common optical analyses

.. codeauthor: Michael J. Hayford
"""

import numpy as np
from matplotlib.figure import Figure

from rayoptics.mpl.axisarrayfigure import Fit

from rayoptics.optical.opticalspec import Field
from rayoptics.optical.trace import trace_astigmatism


class FieldCurveFigure(Figure):
    """ Plot of astigmatism curves """
    def __init__(self, opt_model, eval_fct=trace_astigmatism, **kwargs):
        self.opt_model = opt_model
        self.scale_type = Fit.All
        self.eval_fct = eval_fct

        Figure.__init__(self, **kwargs)

        self.update_data()

    def refresh(self):
        self.update_data()
        self.plot()

    def update_data(self):
        self.s_data = []
        self.t_data = []
        self.field_data = []

        osp = self.opt_model.optical_spec
        _, wvl, foc = osp.lookup_fld_wvl_focus(0)
        fld = Field()
        max_field = osp.field_of_view.max_field()[0]
        for f in np.linspace(0., max_field, num=11):
            fld.y = f
            s_foc, t_foc = self.eval_fct(self.opt_model, fld, wvl, foc)
            self.s_data.append(s_foc)
            self.t_data.append(t_foc)
            self.field_data.append(f)
        return self

    def plot(self):
        self.clf()
        self.ax = self.add_subplot(1, 1, 1)

        self.ax.set_title("Astigmatic Field Plot", pad=10.0, fontsize=18)

        self.ax.plot(self.s_data, self.field_data, label='sagittal')
        self.ax.plot(self.t_data, self.field_data, label='tangential')
        self.ax.set_xlabel('focus')
        self.ax.set_ylabel('field height')

        self.ax.legend()

        self.canvas.draw()

        return self


# experimental - something usable from qt and jupyter
class AnalysisPlot():
    """ abstract api for matplotlib axes customized for specific analyses """
    def __init__(self, ax, opt_model):
        self.opt_model = opt_model
        self.ax = ax

    def refresh(self):
        self.update_data()
        self.plot()

    def update_data(self):
        pass

    def plot(self):
        pass


class AstigmatismCurvePlot(AnalysisPlot):
    def __init__(self, ax, opt_model,
                 eval_fct=trace_astigmatism, **kwargs):
        super(AnalysisPlot, self).__init__(ax, opt_model)
        self.scale_type = Fit.All
        self.eval_fct = eval_fct

        self.update_data()

    def update_data(self):
        self.s_data = []
        self.t_data = []
        self.field_data = []

        osp = self.opt_model.optical_spec
        _, wvl, foc = osp.lookup_fld_wvl_focus(0)
        fld = Field()
        max_field = osp.field_of_view.max_field()[0]
        for f in np.linspace(0., max_field, num=11):
            fld.y = f
            s_foc, t_foc = self.eval_fct(self.opt_model, fld, wvl, foc)
            self.s_data.append(s_foc)
            self.t_data.append(t_foc)
            self.field_data.append(f)

    def plot(self):
        self.ax.cla()

        self.ax.set_title("Astigmatic Field Plot", pad=10.0, fontsize=18)

        self.ax.plot(self.s_data, self.field_data, label='sagittal')
        self.ax.plot(self.t_data, self.field_data, label='tangential')
        self.ax.set_xlabel('focus')
        self.ax.set_ylabel('field height')

        self.ax.legend()
