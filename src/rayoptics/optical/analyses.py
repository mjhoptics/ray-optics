#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2020 Michael J. Hayford
"""Aberration calculations for (fld, wvl, foc), including focus and image shift

    This module refactors some existing ray trace and aberration calculations
    in other modules to be expressed for a single field point and wavelength.
    The ability to apply focus and image shifts to an already acquired data set
    is provided for use interactively and in other performance critical areas.

.. Created on Sat Feb 22 22:01:56 2020

.. codeauthor: Michael J. Hayford
"""
from math import sqrt
import numpy as np
from scipy.interpolate import interp1d
from rayoptics.util.misc_math import normalize

import rayoptics.optical.model_constants as mc

from rayoptics.optical.raytrace import eic_distance
from rayoptics.optical.trace import trace_base, aim_chief_ray, trace_chief_ray


def get_chief_ray_pkg(opt_model, fld, wvl, foc):
    if fld.chief_ray is None:
        aim_chief_ray(opt_model, fld, wvl=wvl)
        chief_ray_pkg = trace_chief_ray(opt_model, fld, wvl, foc)
        fld.chief_ray = chief_ray_pkg
    elif fld.chief_ray[0][2] != wvl:
        chief_ray_pkg = trace_chief_ray(opt_model, fld, wvl, foc)
        fld.chief_ray = chief_ray_pkg
    else:
        chief_ray_pkg = fld.chief_ray
    return chief_ray_pkg


def setup_exit_pupil_coords(opt_model, fld, wvl, foc,
                            chief_ray_pkg, image_pt_2d=None):

    osp = opt_model.optical_spec
    fod = osp.parax_data.fod

    cr, cr_exp_seg = chief_ray_pkg
    # cr_exp_pt: E upper bar prime: pupil center for pencils from Q
    # cr_exp_pt, cr_b4_dir, cr_dst
    # cr_exp_pt = cr_exp_seg[mc.p]

    if image_pt_2d is None:
        # get distance along cr corresponding to a z shift of the defocus
        dist = foc / cr.ray[-1][mc.d][2]
        image_pt = cr.ray[-1][mc.p] - dist*cr.ray[-1][mc.d]
    else:
        image_pt = np.array([image_pt_2d[0], image_pt_2d[1], foc])

    img_pt = np.array(image_pt)
    img_pt[2] += fod.img_dist

    # R' radius of reference sphere for O'
    ref_sphere_vec = img_pt - cr_exp_seg[mc.p]
    ref_sphere_radius = np.linalg.norm(ref_sphere_vec)
    ref_dir = normalize(ref_sphere_vec)

    ref_sphere = (image_pt, ref_dir, ref_sphere_radius)

    return ref_sphere


def wave_abr_full_calc(fod, fld, wvl, foc, ray_pkg, chief_ray_pkg, ref_sphere):
    """Given a ray, a chief ray and an image pt, evaluate the OPD."""
    image_pt, ref_dir, ref_sphere_radius = ref_sphere
    cr, cr_exp_seg = chief_ray_pkg
    chief_ray, chief_ray_op, wvl = cr
    cr_exp_pt, cr_exp_dir, cr_exp_dist = cr_exp_seg

    ray, ray_op, wvl = ray_pkg

    k = -2  # last interface in sequence

    # eq 3.12
    e1 = eic_distance((ray[1][mc.p], ray[0][mc.d]),
                      (chief_ray[1][mc.p], chief_ray[0][mc.d]))
    # eq 3.13
    ekp = eic_distance((ray[k][mc.p], ray[k][mc.d]),
                       (chief_ray[k][mc.p], chief_ray[k][mc.d]))

    dst = ekp - cr_exp_dist

    eic_exp_pt = ray[k][mc.p] - dst*ray[k][mc.d]
    p_coord = eic_exp_pt - cr_exp_pt
    F = ref_dir.dot(ray[k][mc.d]) - ray[k][mc.d].dot(p_coord)/ref_sphere_radius
    J = p_coord.dot(p_coord)/ref_sphere_radius - 2.0*ref_dir.dot(p_coord)
    ep = J/(F + sqrt(F**2 + J/ref_sphere_radius))

    n_obj = abs(fod.n_obj)
    n_img = abs(fod.n_img)
    opd = -n_obj*e1 - ray_op + n_img*ekp + chief_ray_op - n_img*ep
    return opd


def wave_abr_pre_calc(fod, fld, wvl, foc, ray_pkg, chief_ray_pkg):
    """Pre-calculate the part of the OPD calc independent of focus."""
    cr, cr_exp_seg = chief_ray_pkg
    chief_ray, chief_ray_op, wvl = cr
    cr_exp_pt, cr_exp_dir, cr_exp_dist = cr_exp_seg

    ray, ray_op, wvl = ray_pkg

    k = -2  # last interface in sequence

    # eq 3.12
    e1 = eic_distance((ray[1][mc.p], ray[0][mc.d]),
                      (chief_ray[1][mc.p], chief_ray[0][mc.d]))
    # eq 3.13
    ekp = eic_distance((ray[k][mc.p], ray[k][mc.d]),
                       (chief_ray[k][mc.p], chief_ray[k][mc.d]))

    pre_opd = -abs(fod.n_obj)*e1 - ray_op + abs(fod.n_img)*ekp + chief_ray_op

    dst = ekp - cr_exp_dist
    eic_exp_pt = ray[k][mc.p] - dst*ray[k][mc.d]
    p_coord = eic_exp_pt - cr_exp_pt

    return pre_opd, p_coord


def wave_abr_calc(fod, fld, wvl, foc, ray_pkg, pre_opd_pkg, ref_sphere):
    """Given pre-calculated info and a ref. sphere, return the ray's OPD."""
    image_pt, ref_dir, ref_sphere_radius = ref_sphere
    pre_opd, p_coord = pre_opd_pkg
    ray, ray_op, wvl = ray_pkg

    k = -2  # last interface in sequence

    F = ref_dir.dot(ray[k][mc.d]) - ray[k][mc.d].dot(p_coord)/ref_sphere_radius
    J = p_coord.dot(p_coord)/ref_sphere_radius - 2.0*ref_dir.dot(p_coord)
    ep = J/(F + sqrt(F**2 + J/ref_sphere_radius))

    opd = pre_opd - abs(fod.n_img)*ep
    return opd


class RayFan():

    def __init__(self, opt_model, f=0, wl=None, foc=None, image_pt_2d=None,
                 num_rays=21, xyfan='y', **kwargs):
        self.opt_model = opt_model
        osp = opt_model.optical_spec
        self.fld = osp.field_of_view.fields[f] if isinstance(f, int) else f
        self.wvl = osp.spectral_region.central_wvl if wl is None else wl

        self.foc = osp.defocus.focus_shift if foc is None else foc
        self.image_pt_2d = image_pt_2d

        self.num_rays = num_rays

        if xyfan == 'x':
            self.xyfan = 0
        elif xyfan == 'y':
            self.xyfan = 1
        else:
            self.xyfan = int(xyfan)

        self.update_data()

    def update_data(self, build='rebuild'):
        """Set the fan attribute to a list of (pupil coords), dx, dy, opd."""
        if build == 'rebuild':
            self.fan_pkg = trace_fan(
                self.opt_model, self.fld, self.wvl, self.foc, self.xyfan,
                image_pt_2d=self.image_pt_2d, num_rays=self.num_rays)

        self.fan = focus_fan(self.opt_model, self.fan_pkg,
                             self.fld, self.wvl, self.foc,
                             image_pt_2d=self.image_pt_2d)
        return self


def select_plot_data(fan, xyfan, data_type):
    f_x = []
    f_y = []
    for p, val in fan:
        f_x.append(p[xyfan])
        f_y.append(val[data_type])
    f_x = np.array(f_x)
    f_y = np.array(f_y)
    return f_x, f_y


def smooth_plot_data(f_x, f_y, num_points=100):
    interpolator = interp1d(f_x, f_y,
                            kind='cubic', assume_sorted=True)
    x_sample = np.linspace(f_x.min(), f_x.max(), num_points)
    y_fit = interpolator(x_sample)

    return x_sample, y_fit


def trace_ray_fan(opt_model, fan_rng, fld, wvl, foc, **kwargs):
    """ xy determines whether x (=0) or y (=1) fan """
    start = np.array(fan_rng[0])
    stop = fan_rng[1]
    num = fan_rng[2]
    step = (stop - start)/(num - 1)
    fan = []
    for r in range(num):
        pupil = np.array(start)
        ray_pkg = trace_base(opt_model, pupil, fld, wvl, **kwargs)
        fan.append([pupil[0], pupil[1], ray_pkg])
        start += step
    return fan


def eval_fan(opt_model, fld, wvl, foc, xy,
             image_pt_2d=None, num_rays=21):
    """Trace a fan of rays and evaluate dx, dy, & OPD across the fan."""
    fod = opt_model.optical_spec.parax_data.fod
    cr_pkg = get_chief_ray_pkg(opt_model, fld, wvl, foc)
    ref_sphere = setup_exit_pupil_coords(opt_model, fld, wvl, foc, cr_pkg,
                                         image_pt_2d=image_pt_2d)
    fld.chief_ray = cr_pkg
    fld.ref_sphere = ref_sphere

    fan_start = np.array([0., 0.])
    fan_stop = np.array([0., 0.])
    fan_start[xy] = -1.0
    fan_stop[xy] = 1.0
    fan_def = [fan_start, fan_stop, num_rays]

    fan = trace_ray_fan(opt_model, fan_def, fld, wvl, foc)

    convert_to_opd = 1/opt_model.nm_to_sys_units(wvl)

    def rfc(fi):
        pupil_x, pupil_y, ray_pkg = fi
        if ray_pkg is not None:
            image_pt = ref_sphere[0]
            ray = ray_pkg[mc.ray]
            dist = foc / ray[-1][mc.d][2]
            defocused_pt = ray[-1][mc.p] - dist*ray[-1][mc.d]
            t_abr = defocused_pt - image_pt

            opdelta = wave_abr_full_calc(fod, fld, wvl, foc, ray_pkg,
                                         cr_pkg, ref_sphere)
            opd = convert_to_opd*opdelta
            return (pupil_x, pupil_y), (t_abr[0], t_abr[1], opd)
        else:
            return pupil_x, pupil_y, np.NaN
    fan_data = [rfc(i) for i in fan]

    return fan_data


def trace_fan(opt_model, fld, wvl, foc, xy,
              image_pt_2d=None, num_rays=21):
    """Trace a grid of rays and evaluate the OPD across the wavefront."""
    fod = opt_model.optical_spec.parax_data.fod
    cr_pkg = get_chief_ray_pkg(opt_model, fld, wvl, foc)
    ref_sphere = setup_exit_pupil_coords(opt_model, fld, wvl, foc, cr_pkg,
                                         image_pt_2d=image_pt_2d)
    fld.chief_ray = cr_pkg
    fld.ref_sphere = ref_sphere

    fan_start = np.array([0., 0.])
    fan_stop = np.array([0., 0.])
    fan_start[xy] = -1.0
    fan_stop[xy] = 1.0
    fan_def = [fan_start, fan_stop, num_rays]

    fan = trace_ray_fan(opt_model, fan_def, fld, wvl, foc)

    def wpc(fi):
        pupil_x, pupil_y, ray_pkg = fi
        if ray_pkg is not None:
            pre_opd_pkg = wave_abr_pre_calc(fod, fld, wvl, foc,
                                            ray_pkg, cr_pkg)
            return pre_opd_pkg
        else:
            return None

    upd_fan = [wpc(i) for i in fan]

    return fan, upd_fan


def focus_fan(opt_model, fan_pkg, fld, wvl, foc, image_pt_2d=None):
    """Trace a grid of rays and evaluate the OPD across the wavefront."""
    fod = opt_model.optical_spec.parax_data.fod
    fan, upd_fan = fan_pkg
    cr_pkg = get_chief_ray_pkg(opt_model, fld, wvl, foc)
    ref_sphere = setup_exit_pupil_coords(opt_model, fld, wvl, foc, cr_pkg,
                                         image_pt_2d=image_pt_2d)
    convert_to_opd = 1/opt_model.nm_to_sys_units(wvl)

    def rfc(fi, fiu):
        pupil_x, pupil_y, ray_pkg = fi
        if ray_pkg is not None:
            image_pt = ref_sphere[0]
            ray = ray_pkg[mc.ray]
            dist = foc / ray[-1][mc.d][2]
            defocused_pt = ray[-1][mc.p] - dist*ray[-1][mc.d]
            t_abr = defocused_pt - image_pt

            opdelta = wave_abr_calc(fod, fld, wvl, foc,
                                    ray_pkg, fiu, ref_sphere)
            opd = convert_to_opd*opdelta
            return (pupil_x, pupil_y), (t_abr[0], t_abr[1], opd)
        else:
            return pupil_x, pupil_y, np.NaN
    fan_data = [rfc(fi, fiu) for fi, fiu in zip(fan, upd_fan)]
    return fan_data


def trace_ray_grid(opt_model, grid_rng, fld, wvl, foc, **kwargs):
    """Trace a grid of rays at fld and wvl and return ray_pkgs in 2d list."""
    start = np.array(grid_rng[0])
    stop = grid_rng[1]
    num = grid_rng[2]
    step = np.array((stop - start)/(num - 1))
    grid = []
    for i in range(num):
        grid_row = []

        for j in range(num):
            pupil = np.array(start)
            if (pupil[0]**2 + pupil[1]**2) < 1.0:
                ray_pkg = trace_base(opt_model, pupil, fld, wvl, **kwargs)
                grid_row.append([pupil[0], pupil[1], ray_pkg])
            else:  # ray outside pupil
                grid_row.append([pupil[0], pupil[1], None])

            start[1] += step[1]

        grid.append(grid_row)
        start[0] += step[0]
        start[1] = grid_rng[0][1]

    return grid


def eval_wavefront(opt_model, fld, wvl, foc,
                   image_pt_2d=None, num_rays=21):
    """Trace a grid of rays and evaluate the OPD across the wavefront."""
    fod = opt_model.optical_spec.parax_data.fod
    cr_pkg = get_chief_ray_pkg(opt_model, fld, wvl, foc)
    ref_sphere = setup_exit_pupil_coords(opt_model, fld, wvl, foc, cr_pkg,
                                         image_pt_2d=image_pt_2d)
    fld.chief_ray = cr_pkg
    fld.ref_sphere = ref_sphere

    grid_start = np.array([-1., -1.])
    grid_stop = np.array([1., 1.])
    grid_def = [grid_start, grid_stop, num_rays]

    grid = trace_ray_grid(opt_model, grid_def, fld, wvl, foc)

    convert_to_opd = 1/opt_model.nm_to_sys_units(wvl)

    def rfc(gij):
        pupil_x, pupil_y, ray_pkg = gij
        if ray_pkg is not None:
            opdelta = wave_abr_full_calc(fod, fld, wvl, foc, ray_pkg,
                                         cr_pkg, ref_sphere)
            opd = convert_to_opd*opdelta
            return pupil_x, pupil_y, opd
        else:
            return pupil_x, pupil_y, np.NaN
    opd_grid = [[rfc(j) for j in i] for i in grid]

    return np.array(opd_grid)


def trace_wavefront(opt_model, fld, wvl, foc,
                    image_pt_2d=None, num_rays=21):
    """Trace a grid of rays and pre-calculate data needed for rapid refocus."""
    fod = opt_model.optical_spec.parax_data.fod
    cr_pkg = get_chief_ray_pkg(opt_model, fld, wvl, foc)
    ref_sphere = setup_exit_pupil_coords(opt_model, fld, wvl, foc, cr_pkg,
                                         image_pt_2d=image_pt_2d)
    fld.chief_ray = cr_pkg
    fld.ref_sphere = ref_sphere

    grid_start = np.array([-1., -1.])
    grid_stop = np.array([1., 1.])
    grid_def = [grid_start, grid_stop, num_rays]

    grid = trace_ray_grid(opt_model, grid_def, fld, wvl, foc)

    def wpc(gij):
        pupil_x, pupil_y, ray_pkg = gij
        if ray_pkg is not None:
            pre_opd_pkg = wave_abr_pre_calc(fod, fld, wvl, foc,
                                            ray_pkg, cr_pkg)
            return pre_opd_pkg
        else:
            return None
    upd_grid = [[wpc(j) for j in i] for i in grid]

    return grid, upd_grid


def focus_wavefront(opt_model, grid_pkg, fld, wvl, foc, image_pt_2d=None):
    """Given pre-calculated info and a ref. sphere, return the ray's OPD."""
    fod = opt_model.optical_spec.parax_data.fod
    grid, upd_grid = grid_pkg
    cr_pkg = get_chief_ray_pkg(opt_model, fld, wvl, foc)
    ref_sphere = setup_exit_pupil_coords(opt_model, fld, wvl, foc, cr_pkg,
                                         image_pt_2d=image_pt_2d)
    convert_to_opd = 1/opt_model.nm_to_sys_units(wvl)

    def rfc(gij, uij):
        pupil_x, pupil_y, ray_pkg = gij
        if ray_pkg is not None:
            opdelta = wave_abr_calc(fod, fld, wvl, foc,
                                    ray_pkg, uij, ref_sphere)
            opd = convert_to_opd*opdelta
            return pupil_x, pupil_y, opd
        else:
            return pupil_x, pupil_y, np.NaN
    refocused_grid = [[rfc(jg, ju) for jg, ju in zip(ig, iu)]
                      for ig, iu in zip(grid, upd_grid)]

    return np.array(refocused_grid)


class RayGrid():

    def __init__(self, opt_model, f=0, wl=None, foc=None, image_pt_2d=None,
                 num_rays=21, **kwargs):
        self.opt_model = opt_model
        osp = opt_model.optical_spec
        self.fld = osp.field_of_view.fields[f] if isinstance(f, int) else f
        self.wvl = osp.spectral_region.central_wvl if wl is None else wl

        self.foc = osp.defocus.focus_shift if foc is None else foc
        self.image_pt_2d = image_pt_2d

        self.num_rays = num_rays

        self.update_data()

    def update_data(self, build='rebuild'):
        if build == 'rebuild':
            self.grid_pkg = trace_wavefront(
                self.opt_model, self.fld, self.wvl, self.foc,
                image_pt_2d=self.image_pt_2d, num_rays=self.num_rays)

        opd = focus_wavefront(self.opt_model, self.grid_pkg,
                              self.fld, self.wvl, self.foc,
                              image_pt_2d=self.image_pt_2d)

        self.grid = np.rollaxis(opd, 2)

        return self
