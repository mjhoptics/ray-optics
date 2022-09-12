#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
""" mpl implementations of common optical analyses

.. codeauthor: Michael J. Hayford
"""

import abc

import numpy as np

from rayoptics.mpl.axisarrayfigure import Fit
from rayoptics.mpl.styledfigure import StyledFigure

from rayoptics.raytr.opticalspec import Field
from rayoptics.raytr.trace import trace_astigmatism, trace_astigmatism_curve
from rayoptics.parax.thirdorder import compute_third_order


class FieldCurveFigure(StyledFigure):
    """ Plot of astigmatism curves """

    def __init__(self, opt_model,
                 user_scale_value=0.1, num_points=21,
                 **kwargs):
        self.opt_model = opt_model
        self.scale_type = Fit.All
        self.user_scale_value = user_scale_value
        self.num_points = num_points

        super().__init__(**kwargs)

        self.update_data()

    def refresh(self, **kwargs):
        self.update_data(**kwargs)
        self.plot()
        return self

    def update_data(self, **kwargs):
        results = trace_astigmatism_curve(self.opt_model, self.num_points)
        self.field_data, self.s_data, self.t_data = results
        return self

    def plot(self):
        self.clf()
        self.ax = self.add_subplot(1, 1, 1)

        self.ax.set_title("Astigmatic Field Plot", pad=10.0, fontsize=18)

        self.ax.plot(self.s_data, self.field_data,
                     label='sagittal',
                     c=self._rgb['violet'])
        self.ax.plot(self.t_data, self.field_data,
                     label='tangential',
                     c=self._rgb['magenta'])
        self.ax.set_xlabel('focus')
        self.ax.set_ylabel('field height')

        if self.scale_type == Fit.All:
            pass
        if self.scale_type == Fit.User_Scale and self.user_scale_value is not None:
            us = self.user_scale_value
            self.ax.set_xlim(-us, us)

        self.ax.legend()

        self.canvas.draw()

        return self


class ThirdOrderBarChart(StyledFigure):
    def __init__(self, opt_model,
                 user_scale_value=0.1,
                 **kwargs):
        super().__init__(**kwargs)
        self.opt_model = opt_model
        self.scale_type = Fit.All
        self.user_scale_value = user_scale_value

        self.update_data()

    def refresh(self, **kwargs):
        self.update_data(**kwargs)
        self.plot()
        return self

    def update_data(self, **kwargs):
        self.to_pkg = compute_third_order(self.opt_model)

    def plot(self):
        self.clf()
        self.ax = self.add_subplot(1, 1, 1)

        self.ax.set_xlabel('Surface')
        self.ax.set_ylabel('third order aberration')
        self.ax.set_title('Surface by surface third order aberrations')
        self.to_pkg.plot.bar(ax=self.ax, rot=0)
        self.ax.grid(True)

        if self.scale_type == Fit.All:
            pass
        if self.scale_type == Fit.User_Scale and self.user_scale_value is not None:
            us = self.user_scale_value
            self.ax.set_ylim(-us, us)

        self.tight_layout()

        self.canvas.draw()

        return self


# experimental - something usable from qt and jupyter
class AnalysisPlot(abc.ABC):
    """ abstract api for matplotlib axes customized for specific analyses """

    def __init__(self, opt_model):
        self.opt_model = opt_model

    def refresh(self, **kwargs):
        """ called by the app manager to refresh the plot """
        self.update_data(**kwargs)
        self.plot()
        return self

    @abc.abstractmethod
    def update_data(self):
        """ function to update the backend data needed for the plot """

    @abc.abstractmethod
    def plot(self):
        """ function that executes the plotting commands """


class AstigmatismCurvePlot(AnalysisPlot):
    def __init__(self, opt_model, eval_fct=trace_astigmatism, **kwargs):
        super().__init__(opt_model)
        self.scale_type = Fit.All
        self.eval_fct = eval_fct

        self.update_data()

    def update_data(self, **kwargs):
        self.s_data = []
        self.t_data = []
        self.field_data = []

        osp = self.opt_model.optical_spec
        _, wvl, foc = osp.lookup_fld_wvl_focus(0)
        fld = Field()
        max_field = osp.field_of_view.max_field()[0]
        for f in np.linspace(0., max_field, num=21):
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

#        fig.canvas.draw()
