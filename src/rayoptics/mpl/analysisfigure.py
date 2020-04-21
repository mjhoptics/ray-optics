#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2020 Michael J. Hayford
"""Single panel, multiplot MPL figure, with support for line and surface plots

.. Created on Tue Mar 17 16:21:27 2020

.. codeauthor: Michael J. Hayford
"""

import numpy as np
from matplotlib.colors import PowerNorm

from rayoptics.mpl.styledfigure import StyledFigure

from rayoptics.optical import analyses


class AnalysisFigure(StyledFigure):
    """Containing Figure for single panel plots, supports update_data."""

    def __init__(self, data_objs, **kwargs):

        self.data_objs = data_objs
        self.subplots = []
        super().__init__(**kwargs)

    def refresh(self, build='rebuild'):
        self.update_data(build=build)
        self.plot()

    def update_data(self, build='rebuild'):
        for data_obj in self.data_objs:
            data_obj.update_data(build=build)

        return self

    def plot(self):
        self.clf()
        for subplot in self.subplots:
            subplot.plot()

        self.canvas.draw_idle()


class RayFanPlot():
    """Single axis line plot, supporting data display from multiple RayFans.

    Attributes:
        fig: the parent figure for this panel
        gs: gridspec for this panel
        fan_list: list of (fan, data_type, kwargs)

            - fan: a RayFan instance
            - data_type: 'x', 'y', 'opd' to be extracted from fan
            - kwargs: passed to axis.plot() call

        scale_type: if 'fit', set scale to encompass largest data value
        user_scale_value: max scale to apply if scale_type is 'user'
        title: title, if desired, of this plot panel
        yaxis_ticks_position: 'left' or 'right', default is 'left'
        kwargs: passed to plot call

    """

    def __init__(self, fig, gs, fan_list,
                 user_scale_value=0.1, scale_type='fit',
                 yaxis_ticks_position='left',
                 **kwargs):
        self.fig = fig
        self.fig.subplots.append(self)

        self.gs = gs
        self.fan_list = fan_list

        if 'title' in kwargs:
            self.title = kwargs.pop('title', None)
        self.user_scale_value = user_scale_value
        self.scale_type = scale_type
        self.yaxis_ticks_position = yaxis_ticks_position

        self.update_data()

    def init_axis(self, ax):
        ax.grid(True)
        ax.set_xlim(-1., 1.)
        ax.axvline(0, c=self.fig._rgb['foreground'], lw=1)
        ax.axhline(0, c=self.fig._rgb['foreground'], lw=1)
        ax.tick_params(labelbottom=False)
        ax.yaxis.set_ticks_position(self.yaxis_ticks_position)
        if self.title is not None:
            ax.set_title(self.title, fontsize='small')
        return ax

    def refresh(self, build='rebuild'):
        self.update_data(build=build)
        self.plot()

    def update_data(self, build='rebuild'):
        return self

    def plot(self):
        ax = self.fig.add_subplot(self.gs)
        self.init_axis(ax)

        for fan_pkg in self.fan_list:
            fan, data_type, kws = fan_pkg
            if data_type == 'dx':
                dt = 0
            elif data_type == 'dy':
                dt = 1
            elif data_type == 'opd':
                dt = 2
            fan_plot = analyses.select_plot_data(fan.fan, fan.xyfan, dt)
            if 'num_points' in kws:
                num_pts = kws.pop('num_points')
                fan_plot = analyses.smooth_plot_data(*fan_plot,
                                                     num_points=num_pts)

            ax.plot(*fan_plot, **kws)  # x_data, y_data = fan_plot

        if self.scale_type == 'fit':
            y_data = fan_plot[1]
            max_value = max(np.nanmax(y_data), -np.nanmin(y_data))
            ax.set_ylim(-max_value, max_value)
        elif self.user_scale_value is not None:
            us = self.user_scale_value
            ax.set_ylim(-us, us)

        return self


class RayGeoPSF():
    """Single axis spot diagram or 2d histogram.

    Attributes:
        fig: the parent figure for this panel
        gs: gridspec for this panel
        ray_list: a RayList instance
        dsp_type: display type, either 'spot' or 'hist2d'
        scale_type: if 'fit', set scale to encompass largest data value
        user_scale_value: max scale to apply if scale_type is 'user'
        title: title, if desired, of this plot panel
        yaxis_ticks_position: 'left' or 'right', default is 'left'
        kwargs: passed to plot call
    """

    def __init__(self, fig, gs, ray_list,
                 user_scale_value=0.1, scale_type='fit',
                 yaxis_ticks_position='left', dsp_typ='hist2d',
                 **kwargs):
        self.fig = fig
        self.fig.subplots.append(self)

        self.gs = gs
        self.ray_list = ray_list
        self.dsp_typ = dsp_typ

        if 'title' in kwargs:
            self.title = kwargs.pop('title', None)

        if 'norm' in kwargs:
            self.norm = kwargs.pop('norm', None)
        else:
            gamma = kwargs.pop('gamma', 0.5)
            vmax = kwargs.get('vmax') if 'vmax' in kwargs else None
            self.norm = PowerNorm(gamma, vmin=0., vmax=vmax)

        self.plot_kwargs = kwargs

        self.user_scale_value = user_scale_value
        self.scale_type = scale_type
        self.yaxis_ticks_position = yaxis_ticks_position

        self.update_data()

    def init_axis(self, ax):
        if self.dsp_typ == 'spot':
            ax.grid(True)
            linewidth = 0.5
            ax.axvline(0, c=self.fig._rgb['foreground'], lw=linewidth)
            ax.axhline(0, c=self.fig._rgb['foreground'], lw=linewidth)

        ax.yaxis.set_ticks_position(self.yaxis_ticks_position)
        if self.title is not None:
            ax.set_title(self.title, fontsize='small')
        return ax

    def refresh(self, build='rebuild'):
        self.update_data(build=build)
        self.plot()

    def update_data(self, build='rebuild'):
        self.ray_list.update_data(build=build)
        return self

    def plot(self):
        ax = self.fig.add_subplot(self.gs)
        self.init_axis(ax)

        if self.scale_type == 'fit':
            x_data = self.ray_list.ray_abr[0]
            y_data = self.ray_list.ray_abr[1]
            max_value = max(max(np.nanmax(x_data), -np.nanmin(x_data)),
                            max(np.nanmax(y_data), -np.nanmin(y_data)))
        elif self.user_scale_value is not None:
            max_value = self.user_scale_value

        ax.set_xlim(-max_value, max_value)
        ax.set_ylim(-max_value, max_value)

        if self.dsp_typ == 'spot':
            ax.scatter(*self.ray_list.ray_abr, **self.plot_kwargs)
        elif self.dsp_typ == 'hist2d':
            edges = np.linspace(-max_value, max_value, num=100)
            h, xedges, yedges, qmesh = ax.hist2d(*self.ray_list.ray_abr,
                                                 bins=edges,
                                                 norm=self.norm,
                                                 **self.plot_kwargs)
            ax.set_facecolor(qmesh.cmap(0))

        ax.set_aspect('equal')

        return self


class Wavefront():
    """Single axis wavefront map.

    Attributes:
        fig: the parent figure for this panel
        gs: gridspec for this panel
        ray_grid: a RayGrid instance
        do_contours: if True, display contour plot, else plot data grid as an
        image
        scale_type: if 'fit', set scale to encompass largest data value
        user_scale_value: max scale to apply if scale_type is 'user'
        title: title, if desired, of this plot panel
        yaxis_ticks_position: 'left' or 'right', default is 'left'
        cmap: color map for plot, defaults to 'RdBu_r'
        kwargs: passed to plot call
    """

    def __init__(self, fig, gs, ray_grid,
                 do_contours=False, user_scale_value=None,
                 **kwargs):
        self.fig = fig
        self.fig.subplots.append(self)

        self.gs = gs
        self.ray_grid = ray_grid

        if 'title' in kwargs:
            self.title = kwargs.pop('title', None)
        kwargs['cmap'] = kwargs.get('cmap', "RdBu_r")
        self.plot_kwargs = kwargs

        self.do_contours = do_contours
        self.user_scale = user_scale_value

        self.update_data()

    def init_axis(self, ax):
        ax.set_xlim(-1., 1.)
        ax.set_ylim(-1., 1.)
        ax.tick_params(labelbottom=False, labelleft=False)
        if self.title is not None:
            ax.set_title(self.title, fontsize='small')
        return ax

    def refresh(self, build='rebuild'):
        self.update_data(build=build)
        self.plot()

    def update_data(self, build='rebuild'):
        self.ray_grid.update_data(build=build)
        return self

    def plot(self):
        ax = self.fig.add_subplot(self.gs)
        self.init_axis(ax)

        grid = self.ray_grid.grid
        max_value = max(np.nanmax(grid[2]), -np.nanmin(grid[2]))
        max_value = self.user_scale if self.user_scale else max_value

        if self.do_contours:
            hmap = ax.contourf(grid[0],
                               grid[1],
                               grid[2],
                               vmin=-max_value,
                               vmax=max_value,
                               **self.plot_kwargs)
        else:
            # transpose and origin=lower needed to match display orientation
            # of contours
            hmap = ax.imshow(np.transpose(grid[2]),
                             origin='lower',
                             vmin=-max_value,
                             vmax=max_value,
                             extent=[-1., 1., -1., 1.],
                             **self.plot_kwargs
                             )
        self.fig.colorbar(hmap, ax=ax, use_gridspec=True)
        ax.set_aspect('equal')

        return self
