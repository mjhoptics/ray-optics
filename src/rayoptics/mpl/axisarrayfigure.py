#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
"""

.. Created on Fri Apr  13 10:43:19 2018

.. codeauthor: Michael J. Hayford
"""

from enum import Enum, auto

from scipy.interpolate import interp1d

import numpy as np

from rayoptics.mpl.styledfigure import StyledFigure

from rayoptics.raytr.waveabr import wave_abr_full_calc
import rayoptics.optical.model_constants as mc


class Fit(Enum):
    All = auto()
    All_Same = auto()
    User_Scale = auto()


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


class AxisArrayFigure(StyledFigure):

    def __init__(self, opt_model,
                 num_rays=21,
                 scale_type=Fit.All,
                 user_scale_value=0.1,
                 num_rows=1, num_cols=1,
                 eval_fct=None, **kwargs):
        self.opt_model = opt_model
        self.num_rays = num_rays
        self.user_scale_value = user_scale_value
        self.scale_type = scale_type

        super().__init__(**kwargs)

        self.eval_fct = eval_fct
        self.num_rows = num_rows
        self.num_cols = num_cols
        self.update_data(**kwargs)

    def init_axis(self, ax):
        ax.grid(True)
        ax.set_xlim(-1., 1.)
        ax.axvline(0, c=self._rgb['foreground'], lw=1)
        ax.axhline(0, c=self._rgb['foreground'], lw=1)

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

    def wvl_to_sys_units(self, wvl):
        return self.opt_model.nm_to_sys_units(wvl)

    def refresh(self, **kwargs):
        self.update_data(**kwargs)
        self.plot()
        return self

    def update_data(self, **kwargs):
        pass

    def plot(self):
        pass


class RayFanFigure(AxisArrayFigure):

    def __init__(self, opt_model, data_type, override_style=True,
                 do_smoothing=True, **kwargs):
        self.max_value_all = 0.0
        self.override_style = override_style
        self.do_smoothing = do_smoothing
        seq_model = opt_model.seq_model
        osp = opt_model.optical_spec
        central_wvl = osp.spectral_region.central_wvl
        convert_to_waves = 1/opt_model.nm_to_sys_units(central_wvl)

        def ray_abr(p, xy, ray_pkg, fld, wvl, foc):
            if ray_pkg[mc.ray] is not None:
                image_pt = fld.ref_sphere[0]
                ray = ray_pkg[mc.ray]
                dist = foc / ray[-1][mc.d][2]
                defocused_pt = ray[-1][mc.p] + dist*ray[-1][mc.d]
                t_abr = defocused_pt - image_pt
                return t_abr[xy]
            else:
                return None

        def opd(p, xy, ray_pkg, fld, wvl, foc):
            if ray_pkg[mc.ray] is not None:
                fod = opt_model['analysis_results']['parax_data'].fod
                opd = wave_abr_full_calc(fod, fld, wvl, foc, ray_pkg,
                                         fld.chief_ray, fld.ref_sphere)
                return convert_to_waves*opd
            else:
                return None

        def eval_abr_fan(i, j):
            fld, wvl, foc = osp.lookup_fld_wvl_focus(i)
            return seq_model.trace_fan(ray_abr, i, j, num_rays=self.num_rays)

        def eval_opd_fan(i, j):
            fld, wvl, foc = osp.lookup_fld_wvl_focus(i)
            return seq_model.trace_fan(opd, i, j, num_rays=self.num_rays)

        if data_type == 'Ray':
            eval_fan = eval_abr_fan
        elif data_type == 'OPD':
            eval_fan = eval_opd_fan
        num_flds = len(osp.field_of_view.fields)
        super().__init__(opt_model, eval_fct=eval_fan,
                         num_rows=num_flds, num_cols=2, **kwargs)

    def update_data(self, build='rebuild', **kwargs):
        do_smoothing = kwargs.get('do_smoothing', self.do_smoothing)
        self.axis_data_array = []
        for i in reversed(range(self.num_rows)):
            row = []
            for j in reversed(range(self.num_cols)):
                x_smooth = []
                y_smooth = []
                x_data, y_data, max_value, rc = self.eval_fct(i, j)
#                x_data, y_data, max_value, rc = self.eval_axis_data(i, j)
#                rc = clip_to_range(rc, 0.0, 1.0)
                for k in range(len(x_data)):
                    if do_smoothing:
                        interpolator = interp1d(x_data[k], y_data[k],
                                                kind='cubic',
                                                assume_sorted=True)
                        x_samp = np.linspace(x_data[k].min(),
                                             x_data[k].max(), 100)
                        y_fit = interpolator(x_samp)
                        x_smooth.append(x_samp)
                        y_smooth.append(y_fit)
                    else:
                        x_smooth.append(x_data[k])
                        y_smooth.append(y_data[k])
                        
                row.append((x_smooth, y_smooth, max_value, rc))
            self.axis_data_array.append(row)
        return self

    def plot(self):
        if hasattr(self, 'ax_arr'):
            self.clf()

        m = self.num_rows - 1
        n = self.num_cols - 1
        self.ax_arr = self.construct_plot_array(self.num_rows, self.num_cols)

        self.max_rho_all = 0.0
        self.max_value_all = 0.0
        for i in reversed(range(self.num_rows)):
            for j in reversed(range(self.num_cols)):
                x_data, y_data, max_value, rc = self.axis_data_array[m-i][n-j]
                ax = self.ax_arr[m-i][n-j]
                for k in range(len(x_data)):
                    if self.override_style:
                        ax.plot(x_data[k], y_data[k], c=rc[k])
                    else:
                        ax.plot(x_data[k], y_data[k])

                max_rho_val, max_y_val = max_value
                if max_rho_val > self.max_rho_all:
                    self.max_rho_all = max_rho_val
                if max_y_val > self.max_value_all:
                    self.max_value_all = max_y_val

        rv = self.max_rho_all
        [[ax.set_xlim(-rv, rv) for ax in r] for r in self.ax_arr]
        if self.scale_type == Fit.All:
            pass
#            print("Fit.All", self.max_value_all)
#            [[ax.set_ylim(-mv, mv) for ax in r] for r in self.ax_arr]
        if self.scale_type == Fit.All_Same:
            mv = self.max_value_all
#            print("Fit.All_Same", mv)
            [[ax.set_ylim(-mv, mv) for ax in r] for r in self.ax_arr]
        if self.scale_type == Fit.User_Scale and self.user_scale_value is not None:
            us = self.user_scale_value
#            print("User_Scale", us)
            [[ax.set_ylim(-us, us) for ax in r] for r in self.ax_arr]

        self.canvas.draw()

        return self


class SpotDiagramFigure(AxisArrayFigure):

    def __init__(self, opt_model, **kwargs):
        self.max_value_all = 0.0
        seq_model = opt_model.seq_model
        osp = opt_model.optical_spec

        def spot(p, wi, ray_pkg, fld, wvl, foc):
            if ray_pkg is not None:
                image_pt = fld.ref_sphere[0]
                ray = ray_pkg[mc.ray]
                dist = foc / ray[-1][mc.d][2]
                defocused_pt = ray[-1][mc.p] + dist*ray[-1][mc.d]
                t_abr = defocused_pt - image_pt
                return np.array([t_abr[0], t_abr[1]])
            else:
                return None

        def eval_grid(i, j):
            fld, wvl, foc = osp.lookup_fld_wvl_focus(i, wl=j)
            return seq_model.trace_grid(spot, i, num_rays=self.num_rays,
                                        form='list', append_if_none=False)

        num_flds = len(osp.field_of_view.fields)
        super().__init__(opt_model, eval_fct=eval_grid,
                         num_rows=num_flds, num_cols=1, **kwargs)

    def init_axis(self, ax):
        ax.grid(True)
        ax.axvline(0, c=self._rgb['foreground'], lw=1)
        ax.axhline(0, c=self._rgb['foreground'], lw=1)

    def update_data(self, build='rebuild', **kwargs):
        self.axis_data_array = []
        for i in reversed(range(self.num_rows)):
            row = []
            for j in reversed(range(self.num_cols)):
                grids, rc = self.eval_fct(i, j)
                max_val = max([max(np.max(g), -np.min(g)) for g in grids])
                row.append((grids, max_val, rc))
            self.axis_data_array.append(row)
        return self

    def plot(self):
        if hasattr(self, 'ax_arr'):
            self.clf()

        m = self.num_rows - 1
        n = self.num_cols - 1
        self.ax_arr = self.construct_plot_array(self.num_rows, self.num_cols)

        self.max_value_all = 0.0
        for i in reversed(range(self.num_rows)):
            for j in reversed(range(self.num_cols)):
                grids, max_value, rc = self.axis_data_array[m-i][n-j]
                ax = self.ax_arr[m-i][n-j]
                for k in range(len(rc)):
                    ax.plot(np.transpose(grids[k])[0],
                            np.transpose(grids[k])[1],
                            c=rc[k], linestyle='None',
                            marker='o', markersize=2)
                    ax.set_aspect('equal')

                if max_value > self.max_value_all:
                    self.max_value_all = max_value

        if self.scale_type == Fit.All:
            pass
#            print("Fit.All", self.max_value_all)
#            [[ax.set_ylim(-mv, mv) for ax in r] for r in self.ax_arr]
        if self.scale_type == Fit.All_Same:
            mv = self.max_value_all
#            print("Fit.All_Same", mv)
            [[ax.set_xlim(-mv, mv) for ax in r] for r in self.ax_arr]
            [[ax.set_ylim(-mv, mv) for ax in r] for r in self.ax_arr]
        if self.scale_type == Fit.User_Scale and self.user_scale_value is not None:
            us = self.user_scale_value
#            print("User_Scale", us)
            [[ax.set_xlim(-us, us) for ax in r] for r in self.ax_arr]
            [[ax.set_ylim(-us, us) for ax in r] for r in self.ax_arr]

        self.tight_layout()

        self.canvas.draw()

        return self


class WavefrontFigure(AxisArrayFigure):

    def __init__(self, opt_model, **kwargs):
        self.max_value_all = 0.0
        seq_model = opt_model.seq_model
        osp = opt_model.optical_spec
        central_wvl = osp.spectral_region.central_wvl
        convert_to_waves = 1/opt_model.nm_to_sys_units(central_wvl)

        def wave(p, wi, ray_pkg, fld, wvl, foc):
            x = p[0]
            y = p[1]
            if ray_pkg is not None:
                fod = opt_model['analysis_results']['parax_data'].fod
                opd = wave_abr_full_calc(fod, fld, wvl, foc, ray_pkg,
                                         fld.chief_ray, fld.ref_sphere)
                opd = convert_to_waves*opd
            else:
                opd = 0.0
            return np.array([x, y, opd])

        def eval_grid(i, j):
            fld, wvl, foc = osp.lookup_fld_wvl_focus(i, wl=j)
            return seq_model.trace_grid(wave, i, wl=j, num_rays=self.num_rays,
                                        form='grid', append_if_none=True)

        num_flds = len(osp.field_of_view.fields)
        num_wvls = len(osp.spectral_region.wavelengths)
        super().__init__(opt_model, eval_fct=eval_grid,
                         num_rows=num_flds, num_cols=num_wvls, **kwargs)

    def init_axis(self, ax):
#        ax.grid(True)
        ax.set_xlim(-1., 1.)
        ax.set_ylim(-1., 1.)
#        ax.axvline(0, c='black', lw=1)
#        ax.axhline(0, c='black', lw=1)

    def update_data(self, build='rebuild', **kwargs):
        self.axis_data_array = []
        self.max_value_all = 0.0
        for i in reversed(range(self.num_rows)):
            row = []
            for j in reversed(range(self.num_cols)):
                grids, rc = self.eval_fct(i, j)
                g = grids[0]
                g = np.rollaxis(g, 2)
                max_value = max(np.max(g[2]), -np.min(g[2]))
                row.append((g, max_value, rc))

                if max_value > self.max_value_all:
                    self.max_value_all = max_value

            self.axis_data_array.append(row)
        return self

    def plot(self):
        if hasattr(self, 'ax_arr'):
            self.clf()

        m = self.num_rows - 1
        n = self.num_cols - 1
        self.ax_arr = self.construct_plot_array(self.num_rows, self.num_cols)

        for i in reversed(range(self.num_rows)):
            for j in reversed(range(self.num_cols)):
                grid, max_value, rc = self.axis_data_array[m-i][n-j]
                if self.scale_type == Fit.All_Same:
                    max_value = self.max_value_all
                elif (self.scale_type == Fit.User_Scale and
                      self.user_scale_value is not None):
                    max_value = self.user_scale_value
                ax = self.ax_arr[m-i][n-j]
                hmap = ax.contourf(grid[0],
                                   grid[1],
                                   grid[2],
                                   cmap="RdBu_r",
                                   vmin=-max_value,
                                   vmax=max_value)
                # hmap = ax.imshow(np.transpose(grid[2]),
                #                  cmap="RdBu_r",
                #                  vmin=-max_value,
                #                  vmax=max_value,
                #                  extent=[-1., 1., -1., 1.],
                #                  origin='lower')
                self.colorbar(hmap, ax=ax)

                ax.set_aspect('equal')

        self.tight_layout()

        self.canvas.draw()

        return self
