#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2020 Michael J. Hayford
"""Aberration calculations for (fld, wvl, foc), including focus and image shift

    This module refactors some existing ray trace and aberration calculations
    in other modules to be expressed for a single field point and wavelength.
    The ability to apply focus and image shifts to an already acquired data set
    is provided for use interactively and in other performance critical areas.

    The following classes are implemented in this module:

        - :class:`~.Ray`: trace a single ray
        - :class:`~.RayFan`: trace a fan of rays in either the x or y meridian
        - :class:`~.RayList`: trace a list of rays from an object point
        - :class:`~.RayGrid`: trace a rectilinear grid of rays

    All but the `Ray` class are supported by a group of functions to trace the
    rays, accumulate the data (trace_*), and refocus (focus_*) the data. A
    all-in-one function (eval_*) to trace and apply focus is supplied also.
    These are used in the update_data methods of the classes to generate the
    ray data.

    This module also has functions to calculate chief ray and reference sphere
    information as well as functions for calculating the monochromatic PSF of
    the model.

.. Created on Sat Feb 22 22:01:56 2020

.. codeauthor: Michael J. Hayford
"""
import numpy as np
from numpy.fft import fftshift, fft2

from scipy.interpolate import interp1d

import rayoptics.optical.model_constants as mc

from rayoptics.raytr import sampler
from rayoptics.raytr import trace
from rayoptics.raytr import traceerror as terr
from rayoptics.raytr import waveabr


# --- Single ray
class Ray():
    """A ray at the given field and wavelength.

    Attributes:
        opt_model: :class:`~.OpticalModel` instance
        p: relative 2d pupil coordinates
        f: index into :class:`~.FieldSpec` or a :class:`~.Field` instance
        wl: wavelength (nm) to trace the ray, or central wavelength if None
        foc: focus shift to apply to the results
        image_pt_2d: base image point. if None, the chief ray is used
        image_delta: image offset to apply to image_pt_2d
        srf_save:

            'single': save the ray data for surface srf_indx
            'all': save all of the surface by surface ray data

        srf_indx: for single surface retention, the surface index to save
    """

    def __init__(self, opt_model, p, f=0, wl=None, foc=None, image_pt_2d=None,
                 image_delta=None, srf_indx=-1, srf_save='single', 
                 output_filter=None, rayerr_filter=None, color=None, 
                 clip_rays=False):
        self.opt_model = opt_model
        osp = opt_model.optical_spec
        self.pupil = p
        self.fld = osp['fov'].fields[f] if isinstance(f, int) else f
        self.wvl = osp['wvls'].central_wvl if wl is None else wl

        self.foc = osp['focus'].focus_shift if foc is None else foc
        self.image_pt_2d = image_pt_2d
        self.image_delta = image_delta

        self.output_filter = output_filter
        self.rayerr_filter = rayerr_filter
        self.clip_rays = clip_rays

        self.color = color

        self.srf_save = srf_save
        self.srf_indx = srf_indx

        self.update_data()

    def update_data(self, **kwargs):
        """Trace the ray and calculate transverse aberrations. """
        ref_sphere, cr_pkg = trace.setup_pupil_coords(
            self.opt_model, self.fld, self.wvl, self.foc, 
            image_pt=self.image_pt_2d, image_delta=self.image_delta
            )
        build = kwargs.pop('build', 'rebuild')
        if build == 'rebuild':
            ray_result = trace.trace_safe(
                self.opt_model, self.pupil, self.fld, self.wvl, 
                self.output_filter, self.rayerr_filter, 
                use_named_tuples=True, check_apertures=self.clip_rays, 
                **kwargs)
            ray_pkg, ray_err = ray_result
            self.ray_seg = ray_pkg.ray[self.srf_indx]

            if self.srf_save == 'all':
                self.ray_pkg = ray_pkg


        ray_seg = self.ray_seg
        dist = self.foc / ray_seg[mc.d][2]
        defocused_pt = ray_seg[mc.p] + dist*ray_seg[mc.d]
        reference_image_pt = ref_sphere[0]
        self.t_abr = defocused_pt[:2] - reference_image_pt[:2]

        return self


# --- Fan of rays
class RayFan():
    """A fan of rays across the pupil at the given field and wavelength.

    Attributes:
        opt_model: :class:`~.OpticalModel` instance
        f: index into :class:`~.FieldSpec` or a :class:`~.Field` instance
        wl: wavelength (nm) to trace the fan, or central wavelength if None
        foc: focus shift to apply to the results
        image_pt_2d: base image point. if None, the chief ray is used
        image_delta: image offset to apply to image_pt_2d
        num_rays: number of samples along the fan
        xyfan: 'x' or 'y', specifies the axis the fan is sampled on
    """

    def __init__(self, opt_model, f=0, wl=None, foc=None, image_pt_2d=None,
                 image_delta=None, num_rays=21, xyfan='y', output_filter=None,
                 rayerr_filter=None, color=None, clip_rays=False, 
                 **kwargs):
        self.opt_model = opt_model
        osp = opt_model.optical_spec
        self.fld = osp.field_of_view.fields[f] if isinstance(f, int) else f
        self.wvl = osp.spectral_region.central_wvl if wl is None else wl

        self.foc = osp.defocus.focus_shift if foc is None else foc
        self.image_pt_2d = image_pt_2d
        self.image_delta = image_delta

        self.num_rays = num_rays

        if xyfan == 'x':
            self.xyfan = 0
        elif xyfan == 'y':
            self.xyfan = 1
        else:
            self.xyfan = int(xyfan)

        self.color = color

        self.rt_kwargs = kwargs
        self.rt_kwargs['output_filter'] = output_filter
        self.rt_kwargs['rayerr_filter'] = rayerr_filter
        self.rt_kwargs['check_apertures'] = clip_rays

        self.update_data()

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['opt_model']
        del attrs['fan_pkg']
        return attrs

    def update_data(self, **kwargs):
        """Set the fan attribute to a list of (pupil coords), dx, dy, opd."""
        build = kwargs.get('build', 'rebuild')
        if build == 'rebuild':
            self.fan_pkg = trace_fan(
                self.opt_model, self.fld, self.wvl, self.foc, self.xyfan,
                image_pt_2d=self.image_pt_2d, image_delta=self.image_delta, 
                num_rays=self.num_rays,
                **self.rt_kwargs)

        self.fan = focus_fan(self.opt_model, self.fan_pkg,
                             self.fld, self.wvl, self.foc,
                             image_pt_2d=self.image_pt_2d,
                             image_delta=self.image_delta,
                             **self.rt_kwargs)
        return self


def select_plot_data(fan, xyfan, data_type):
    """Given a fan of data, select the sample points and the resulting data."""
    f_x = []
    f_y = []
    for p, val in fan:
        f_x.append(p[xyfan])
        f_y.append(val[data_type])
    f_x = np.array(f_x)
    f_y = np.array(f_y)
    return f_x, f_y


def smooth_plot_data(f_x, f_y, num_points=100):
    """Interpolate fan data points and return a smoothed version."""
    interpolator = interp1d(f_x, f_y,
                            kind='cubic', assume_sorted=True)
    x_sample = np.linspace(f_x.min(), f_x.max(), num_points)
    y_fit = interpolator(x_sample)

    return x_sample, y_fit


def trace_ray_fan(opt_model, fan_rng, fld, wvl, foc,
                  output_filter=None, rayerr_filter=None, **kwargs):
    """Trace a fan of rays, according to fan_rng. """
    start = np.array(fan_rng[0])
    stop = fan_rng[1]
    num = fan_rng[2]
    step = (stop - start)/(num - 1)
    fan = []
    for r in range(num):
        pupil = np.array(start)
        ray_result = trace.trace_safe(opt_model, pupil, fld, wvl, 
                                      output_filter, rayerr_filter, 
                                      use_named_tuples=True, **kwargs)
        ray_pkg, ray_err = ray_result

        if ray_pkg is not None:
            fan.append([pupil[0], pupil[1], ray_pkg])
        start += step
    return fan


def eval_fan(opt_model, fld, wvl, foc, xy,
             image_pt_2d=None, image_delta=None, num_rays=21,
             output_filter=None, rayerr_filter=None, **kwargs):
    """Trace a fan of rays and evaluate dx, dy, & OPD across the fan."""
    fod = opt_model['analysis_results']['parax_data'].fod
    ref_sphere, cr_pkg = trace.setup_pupil_coords(opt_model, fld, wvl, foc, 
                                                  image_pt=image_pt_2d,
                                                  image_delta=image_delta)
    fld.chief_ray = cr_pkg
    fld.ref_sphere = ref_sphere

    fan_start = np.array([0., 0.])
    fan_stop = np.array([0., 0.])
    fan_start[xy] = -1.0
    fan_stop[xy] = 1.0
    fan_def = [fan_start, fan_stop, num_rays]

    fan = trace_ray_fan(opt_model, fan_def, fld, wvl, foc,
                        output_filter=output_filter,
                        rayerr_filter=rayerr_filter, **kwargs)

    central_wvl = opt_model.optical_spec.spectral_region.central_wvl
    convert_to_opd = 1/opt_model.nm_to_sys_units(central_wvl)

    def rfc(fi):
        pupil_x, pupil_y, ray_pkg = fi
        if ray_pkg is not None:
            image_pt = ref_sphere[0]
            ray = ray_pkg[mc.ray]
            dist = foc / ray[-1][mc.d][2]
            defocused_pt = ray[-1][mc.p] + dist*ray[-1][mc.d]
            t_abr = defocused_pt - image_pt

            opdelta = waveabr.wave_abr_full_calc(fod, fld, wvl, foc, ray_pkg,
                                                 cr_pkg, ref_sphere)
            opd = convert_to_opd*opdelta
            return (pupil_x, pupil_y), (t_abr[0], t_abr[1], opd)
        else:
            return pupil_x, pupil_y, np.NaN
    fan_data = [rfc(i) for i in fan]

    return fan_data


def trace_fan(opt_model, fld, wvl, foc, xy,
              image_pt_2d=None, image_delta=None, num_rays=21,
              output_filter=None, rayerr_filter=None, **kwargs):
    """Trace a fan of rays and precalculate data for rapid refocus later."""
    fod = opt_model['analysis_results']['parax_data'].fod
    ref_sphere, cr_pkg = trace.setup_pupil_coords(opt_model, fld, wvl, foc, 
                                                  image_pt=image_pt_2d,
                                                  image_delta=image_delta)
    fld.chief_ray = cr_pkg
    fld.ref_sphere = ref_sphere

    """ xy determines whether x (=0) or y (=1) fan """
    fan_start = np.array([0., 0.])
    fan_stop = np.array([0., 0.])
    fan_start[xy] = -1.0
    fan_stop[xy] = 1.0
    fan_def = [fan_start, fan_stop, num_rays]

    fan = trace_ray_fan(opt_model, fan_def, fld, wvl, foc,
                        output_filter=output_filter,
                        rayerr_filter=rayerr_filter, **kwargs)

    def wpc(fi):
        pupil_x, pupil_y, ray_pkg = fi
        if ray_pkg is not None and not isinstance(ray_pkg, terr.TraceError):
            pre_opd_pkg = waveabr.wave_abr_pre_calc(fod, fld, wvl, foc,
                                                    ray_pkg, cr_pkg, 
                                                    ref_sphere)
            return pre_opd_pkg
        else:
            return None

    upd_fan = [wpc(i) for i in fan]

    return fan, upd_fan


def focus_fan(opt_model, fan_pkg, fld, wvl, foc, 
              image_pt_2d=None, image_delta=None, **kwargs):
    """Refocus the fan of rays and return the tranverse abr. and OPD."""
    fod = opt_model['analysis_results']['parax_data'].fod
    fan, upd_fan = fan_pkg
    ref_sphere, cr_pkg = trace.setup_pupil_coords(opt_model, fld, wvl, foc, 
                                                  image_pt=image_pt_2d,
                                                  image_delta=image_delta)
    central_wvl = opt_model.optical_spec.spectral_region.central_wvl
    convert_to_opd = 1/opt_model.nm_to_sys_units(central_wvl)

    def rfc(fi, fiu):
        pupil_x, pupil_y, ray_pkg = fi
        if ray_pkg is not None and not isinstance(ray_pkg, terr.TraceError):
            image_pt = ref_sphere[0]
            ray = ray_pkg[mc.ray]
            dist = foc / ray[-1][mc.d][2]
            defocused_pt = ray[-1][mc.p] + dist*ray[-1][mc.d]
            t_abr = defocused_pt - image_pt

            opdelta = waveabr.wave_abr_calc(fod, fld, wvl, foc,
                                            ray_pkg, cr_pkg, fiu, ref_sphere)
            opd = convert_to_opd*opdelta
            return (pupil_x, pupil_y), (t_abr[0], t_abr[1], opd)
        else:
            return pupil_x, pupil_y, np.NaN
    fan_data = [rfc(fi, fiu) for fi, fiu in zip(fan, upd_fan)]
    return fan_data


# --- List of rays
class RayList():
    """Container class for a list of rays produced from a list or generator

    Attributes:
        opt_model: :class:`~.OpticalModel` instance
        pupil_gen: (fct, args, kwargs), where:

            - fct: a function returning a generator returning a 2d coordinate
            - args: passed to fct
            - kwargs: passed to fct

        pupil_coords: list of 2d coordinates. If None, filled in by calling
                      pupil_gen.
        num_rays: number of samples side of grid. Used only if pupil_coords and
                  pupil_gen are None.
        f: index into :class:`~.FieldSpec` or a :class:`~.Field` instance
        wl: wavelength (nm) to trace the fan, or central wavelength if None
        foc: focus shift to apply to the results
        image_pt_2d: base image point. if None, the chief ray is used
        image_delta: image offset to apply to image_pt_2d
        apply_vignetting: whether to apply vignetting factors to pupil coords
    """

    def __init__(self, opt_model,
                 pupil_gen=None, pupil_coords=None, num_rays=21,
                 f=0, wl=None, foc=None, image_pt_2d=None, image_delta=None, 
                 output_filter=None, rayerr_filter=None, clip_rays=False, 
                 apply_vignetting=True, **kwargs):
        self.opt_model = opt_model
        osp = opt_model.optical_spec
        if pupil_coords is not None and pupil_gen is None:
            self.pupil_coords = pupil_coords
            self.pupil_gen = None
        else:
            if pupil_gen is not None:
                self.pupil_gen = pupil_gen
            else:
                grid_start = np.array([-1., -1.])
                grid_stop = np.array([1., 1.])
                grid_def = [grid_start, grid_stop, num_rays]
                self.pupil_gen = (sampler.csd_grid_ray_generator,
                                  (grid_def,), {})
            fct, args, kwa = self.pupil_gen
            self.pupil_coords = fct(*args, **kwa)

        self.fld = osp.field_of_view.fields[f] if isinstance(f, int) else f
        self.wvl = osp.spectral_region.central_wvl if wl is None else wl

        self.foc = osp.defocus.focus_shift if foc is None else foc
        self.image_pt_2d = image_pt_2d
        self.image_delta = image_delta

        self.rt_kwargs = kwargs
        self.rt_kwargs['apply_vignetting'] = apply_vignetting
        self.rt_kwargs['output_filter'] = output_filter
        self.rt_kwargs['rayerr_filter'] = rayerr_filter
        self.rt_kwargs['check_apertures'] = clip_rays

        self.update_data()

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['opt_model']
        del attrs['pupil_gen']
        del attrs['pupil_coords']
        del attrs['ray_list']
        return attrs

    def update_data(self, **kwargs):
        build = kwargs.get('build', 'rebuild')
        if build == 'rebuild':
            if self.pupil_gen:
                fct, args, kwa = self.pupil_gen
                self.pupil_coords = fct(*args, **kwa)

            self.ray_list = trace_pupil_coords(
                self.opt_model, self.pupil_coords,
                self.fld, self.wvl, self.foc,
                image_pt_2d=self.image_pt_2d, 
                image_delta=self.image_delta, 
                **self.rt_kwargs)

        ray_list_data = focus_pupil_coords(
            self.opt_model, self.ray_list,
            self.fld, self.wvl, self.foc,
            image_pt_2d=self.image_pt_2d,
            image_delta=self.image_delta,
            **self.rt_kwargs)

        self.ray_abr = np.rollaxis(ray_list_data, 1)

        return self


def trace_ray_list(opt_model, pupil_coords, fld, wvl, foc,
                   append_if_none=False,
                   output_filter=None, rayerr_filter=None,
                   **kwargs):
    """Trace a list of rays at fld and wvl and return ray_pkgs in a list."""

    ray_list = []
    for pupil in pupil_coords:
        ray_result = trace.trace_safe(opt_model, pupil, fld, wvl, 
                                      output_filter, rayerr_filter, 
                                      **kwargs)
        ray_pkg, ray_err = ray_result
        if ray_pkg is not None:
            ray_list.append([pupil[0], pupil[1], ray_pkg])
        else:  # ray outside pupil or failed
            if append_if_none:
                ray_list.append([pupil[0], pupil[1], None])

    return ray_list


def trace_list_of_rays(opt_model, rays,
                       output_filter=None, rayerr_filter=None,
                       **kwargs):
    """Trace a list of rays (pt, dir, wvl) and return ray_pkgs in a list.

    Args:
        opt_model: :class:`~.OpticalModel` instance
        rays: list of (pt0, dir0, wvl)

            - pt0: starting point in coords of first interface
            - dir0: starting direction cosines in coords of first interface
            - wvl: wavelength in nm

        output_filter: None, "last", or a callable. See below
        **kwargs: kwyword args passed to the trace function

    The output_filter keyword argument controls what ray data is returned to
    the caller.

        - if None, returns the entire traced ray
        - if "last", returns the ray data from the last interface
        - if a callable, it must take a ray_pkg as an argument and return the
          desired data or None

    Returns:
        A list with an entry for each ray in rays
    """
    ray_list = []
    for ray in rays:
        pt0, dir0, wvl = ray
        try:
            ray_pkg = trace.trace(opt_model.seq_model, pt0, dir0, wvl, 
                                  **kwargs)
        except terr.TraceError as rayerr:
            if rayerr_filter is None:
                pass
            elif rayerr_filter == 'full':
                ray_list.append((ray, rayerr))
            elif rayerr_filter == 'summary':
                rayerr.ray_pkg = None
                ray_list.append((ray, rayerr))
            else:
                pass
        else:
            if output_filter is None:
                ray_list.append(ray_pkg)
            elif output_filter == 'last':
                ray, op_delta, wvl = ray_pkg
                ray_list.append((ray[-1], op_delta, wvl))
            else:
                ray_list.append(output_filter(ray_pkg))

    return ray_list


def eval_pupil_coords(opt_model, fld, wvl, foc, image_pt_2d=None, 
                      image_delta=None, num_rays=21, **kwargs):
    """Trace a list of rays and return the transverse abr."""
    ref_sphere, cr_pkg = trace.setup_pupil_coords(opt_model, fld, wvl, foc, 
                                                  image_pt=image_pt_2d,
                                                  image_delta=image_delta)
    fld.chief_ray = cr_pkg
    fld.ref_sphere = ref_sphere

    grid_start = np.array([-1., -1.])
    grid_stop = np.array([1., 1.])
    grid_def = [grid_start, grid_stop, num_rays]

    kwargs['check_apertures'] = kwargs.get('check_apertures', True)
    ray_list = trace_ray_list(opt_model, sampler.grid_ray_generator(grid_def),
                              fld, wvl, foc, **kwargs)

    def rfc(ri):
        pupil_x, pupil_y, ray_pkg = ri
        if ray_pkg is not None:
            image_pt = ref_sphere[0]
            ray = ray_pkg[mc.ray]
            dist = foc / ray[-1][mc.d][2]
            defocused_pt = ray[-1][mc.p] + dist*ray[-1][mc.d]
            t_abr = defocused_pt - image_pt
            return t_abr[0], t_abr[1]
        else:
            return np.NaN
    ray_list_data = [rfc(ri) for ri in ray_list]
    return np.array(ray_list_data)


def trace_pupil_coords(opt_model, pupil_coords, fld, wvl, foc,
                       image_pt_2d=None, image_delta=None, **kwargs):
    """Trace a list of rays and return data needed for rapid refocus."""
    ref_sphere, cr_pkg = trace.setup_pupil_coords(opt_model, fld, wvl, foc, 
                                                  image_pt=image_pt_2d,
                                                  image_delta=image_delta)
    fld.chief_ray = cr_pkg
    fld.ref_sphere = ref_sphere

    kwargs['check_apertures'] = kwargs.get('check_apertures', True)
    ray_list = trace_ray_list(opt_model, pupil_coords,
                              fld, wvl, foc, **kwargs)

    return ray_list


def focus_pupil_coords(opt_model, ray_list, fld, wvl, foc, 
                       image_pt_2d=None, image_delta=None, **kwargs):
    """Given pre-traced rays and a ref. sphere, return the transverse abr."""
    ref_sphere, cr_pkg = trace.setup_pupil_coords(opt_model, fld, wvl, foc, 
                                                  image_pt=image_pt_2d,
                                                  image_delta=image_delta)

    def rfc(ri):
        pupil_x, pupil_y, ray_pkg = ri
        if ray_pkg is not None:
            image_pt = ref_sphere[0]
            ray = ray_pkg[mc.ray]
            dist = foc / ray[-1][mc.d][2]
            defocused_pt = ray[-1][mc.p] + dist*ray[-1][mc.d]
            t_abr = defocused_pt - image_pt
            return t_abr[0], t_abr[1]
        else:
            return np.NaN
    ray_list_data = [rfc(ri) for ri in ray_list]
    return np.array(ray_list_data)


# --- Square grid of rays
class RayGrid():
    """Container for a square grid of rays.

    Attributes:
        opt_model: :class:`~.OpticalModel` instance
        f: index into :class:`~.FieldSpec` or a :class:`~.Field` instance
        wl: wavelength (nm) to trace the fan, or central wavelength if None
        foc: focus shift to apply to the results
        image_pt_2d: base image point. if None, the chief ray is used
        image_delta: image offset to apply to image_pt_2d
        num_rays: number of samples along the side of the grid
    """

    def __init__(self, opt_model, f=0, wl=None, foc=None, image_pt_2d=None,
                 image_delta=None, output_filter=None, rayerr_filter=None, 
                 num_rays=21, clip_rays=True, value_if_none=np.NaN, 
                 **kwargs):
        self.opt_model = opt_model
        osp = opt_model.optical_spec
        self.fld = osp.field_of_view.fields[f] if isinstance(f, int) else f
        self.wvl = osp.spectral_region.central_wvl if wl is None else wl

        self.foc = osp.defocus.focus_shift if foc is None else foc
        self.image_pt_2d = image_pt_2d
        self.image_delta = image_delta

        self.num_rays = num_rays
        self.value_if_none = value_if_none

        self.rt_kwargs = kwargs
        self.rt_kwargs['output_filter'] = output_filter
        self.rt_kwargs['rayerr_filter'] = rayerr_filter
        self.rt_kwargs['check_apertures'] = clip_rays

        self.update_data()

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['opt_model']
        del attrs['grid_pkg']
        return attrs

    def update_data(self, **kwargs):
        build = kwargs.get('build', 'rebuild')
        if build == 'rebuild':
            self.grid_pkg = trace_wavefront(
                self.opt_model, self.fld, self.wvl, self.foc,
                image_pt_2d=self.image_pt_2d, image_delta=self.image_delta, 
                num_rays=self.num_rays,
                **self.rt_kwargs)

        opd = focus_wavefront(self.opt_model, self.grid_pkg,
                              self.fld, self.wvl, self.foc,
                              image_pt_2d=self.image_pt_2d,
                              image_delta=self.image_delta, 
                              value_if_none=self.value_if_none,
                              **self.rt_kwargs)

        self.grid = np.rollaxis(opd, 2)

        return self


def trace_ray_grid(opt_model, grid_rng, fld, wvl, foc, append_if_none=True,
                   output_filter=None, rayerr_filter=None, **kwargs):
    """Trace a grid of rays at fld and wvl and return ray_pkgs in 2d list."""
    start = np.array(grid_rng[0])
    stop = grid_rng[1]
    num = grid_rng[2]
    step = np.array((stop - start)/(num - 1))
    grid = []
    kwargs['apply_vignetting'] = kwargs.get('apply_vignetting', False)
    for i in range(num):
        grid_row = []

        for j in range(num):
            pupil = np.array(start)
            ray_result = trace.trace_safe(opt_model, pupil, fld, wvl, 
                                          output_filter, rayerr_filter, 
                                          **kwargs)
            ray_pkg, ray_err = ray_result
            if ray_pkg is not None:
                    grid_row.append([pupil[0], pupil[1], ray_pkg])
            else:  # ray outside pupil or failed
                if append_if_none:
                    grid_row.append([pupil[0], pupil[1], None])

            start[1] += step[1]

        grid.append(grid_row)
        start[0] += step[0]
        start[1] = grid_rng[0][1]

    return grid


def eval_wavefront(opt_model, fld, wvl, foc, image_pt_2d=None, 
                   image_delta=None, num_rays=21, value_if_none=np.NaN, 
                   **kwargs):
    """Trace a grid of rays and evaluate the OPD across the wavefront."""
    fod = opt_model['analysis_results']['parax_data'].fod
    ref_sphere, cr_pkg = trace.setup_pupil_coords(opt_model, fld, wvl, foc, 
                                                  image_pt=image_pt_2d,
                                                  image_delta=image_delta)
    fld.chief_ray = cr_pkg
    fld.ref_sphere = ref_sphere

    vig_bbox = fld.vignetting_bbox(opt_model['osp']['pupil'])
    vig_grid_def = [vig_bbox[0], vig_bbox[1], num_rays]

    kwargs['check_apertures'] = kwargs.get('check_apertures', True)
    grid = trace_ray_grid(opt_model, vig_grid_def, 
                          fld, wvl, foc, **kwargs)

    central_wvl = opt_model.optical_spec.spectral_region.central_wvl
    convert_to_opd = 1/opt_model.nm_to_sys_units(central_wvl)

    def rfc(gij):
        pupil_x, pupil_y, ray_pkg = gij
        if ray_pkg is not None:
            opdelta = waveabr.wave_abr_full_calc(fod, fld, wvl, foc, ray_pkg,
                                                 cr_pkg, ref_sphere)
            opd = convert_to_opd*opdelta
            return pupil_x, pupil_y, opd
        else:
            return pupil_x, pupil_y, value_if_none
    opd_grid = [[rfc(j) for j in i] for i in grid]

    return np.array(opd_grid)


def trace_wavefront(opt_model, fld, wvl, foc,
                    image_pt_2d=None, image_delta=None, num_rays=21,
                    **kwargs):
    """Trace a grid of rays and pre-calculate data needed for rapid refocus."""
    fod = opt_model['analysis_results']['parax_data'].fod
    ref_sphere, cr_pkg = trace.setup_pupil_coords(opt_model, fld, wvl, foc, 
                                                  image_pt=image_pt_2d,
                                                  image_delta=image_delta)
    fld.chief_ray = cr_pkg
    fld.ref_sphere = ref_sphere

    vig_bbox = fld.vignetting_bbox(opt_model['osp']['pupil'])
    vig_grid_def = [vig_bbox[0], vig_bbox[1], num_rays]

    kwargs['check_apertures'] = kwargs.get('check_apertures', True)
    grid = trace_ray_grid(opt_model, vig_grid_def, 
                          fld, wvl, foc, **kwargs)

    def wpc(gij):
        pupil_x, pupil_y, ray_pkg = gij
        if ray_pkg is not None:
            pre_opd_pkg = waveabr.wave_abr_pre_calc(fod, fld, wvl, foc,
                                                    ray_pkg, cr_pkg, 
                                                    ref_sphere)
            return pre_opd_pkg
        else:
            return None
    upd_grid = [[wpc(j) for j in i] for i in grid]

    return grid, upd_grid


def focus_wavefront(opt_model, grid_pkg, fld, wvl, foc, image_pt_2d=None,
                    image_delta=None, value_if_none=np.NaN, **kwargs):
    """Given pre-traced rays and a ref. sphere, return the ray's OPD."""
    fod = opt_model['analysis_results']['parax_data'].fod
    grid, upd_grid = grid_pkg
    ref_sphere, cr_pkg = trace.setup_pupil_coords(opt_model, fld, wvl, foc, 
                                                  image_pt=image_pt_2d,
                                                  image_delta=image_delta)
    central_wvl = opt_model.optical_spec.spectral_region.central_wvl
    convert_to_opd = 1/opt_model.nm_to_sys_units(central_wvl)

    def rfc(gij, uij):
        pupil_x, pupil_y, ray_pkg = gij
        if ray_pkg is not None:
            opdelta = waveabr.wave_abr_calc(fod, fld, wvl, foc,
                                            ray_pkg, cr_pkg, uij, ref_sphere)
            opd = convert_to_opd*opdelta
            return pupil_x, pupil_y, opd
        else:
            return pupil_x, pupil_y, value_if_none
    refocused_grid = [[rfc(jg, ju) for jg, ju in zip(ig, iu)]
                      for ig, iu in zip(grid, upd_grid)]

    return np.array(refocused_grid)


# --- PSF calculation
def psf_sampling(n=None, n_pupil=None, n_airy=None):
    """Given 2 of 3 parameters, calculate the third.

    Args:
        n: The total width of the sampling grid
        n_pupil: The sampling across the pupil
        n_airy: The sampling across the central peak of the Airy disk

    Returns: (n, n_pupil, n_airy)
    """
    npa = n, n_pupil, n_airy
    i = npa.index(None)
    if i == 0:
        n = round((n_pupil * n_airy)/2.44)
    elif i == 1:
        n_pupil = round(2.44*n/n_airy)
    elif i == 2:
        n_airy = round(2.44*n/n_pupil)
    return n, n_pupil, n_airy


def calc_psf_scaling(pupil_grid, ndim, maxdim):
    """Calculate the input and output grid spacings.

    Args:
        pupil_grid: A RayGrid instance
        ndim: The sampling across the wavefront
        maxdim: The total width of the sampling grid

    Returns:
        delta_x: The linear grid spacing on the entrance pupil
        delta_xp: The linear grid spacing on the image plane
    """
    opt_model = pupil_grid.opt_model
    fod = opt_model['analysis_results']['parax_data'].fod
    wl = opt_model.nm_to_sys_units(pupil_grid.wvl)

    fill_factor = ndim/maxdim
    max_D = 2 * fod.enp_radius / fill_factor
    delta_x = max_D / maxdim
    C = wl/fod.exp_radius

    delta_theta = (fill_factor * C) / 2
    ref_sphere_radius = pupil_grid.fld.ref_sphere[2]
    delta_xp = delta_theta * ref_sphere_radius

    return delta_x, delta_xp


def calc_psf(wavefront, ndim, maxdim):
    """Calculate the point spread function of wavefront W.

    Args:
        wavefront: ndim x ndim Numpy array of wavefront errors. No data
                   condition is indicated by nan
        ndim: The sampling across the wavefront
        maxdim: The total width of the sampling grid

    Returns: AP, the PSF of the input wavefront
    """
    maxdim_by_2 = maxdim//2
    W = np.zeros([maxdim, maxdim])
    nd2 = ndim//2
    W[maxdim_by_2-(nd2-1):maxdim_by_2+(nd2+1),
      maxdim_by_2-(nd2-1):maxdim_by_2+(nd2+1)] = np.nan_to_num(wavefront)

    phase = np.exp(1j*2*np.pi*W)

    for i in range(len(phase)):
        for j in range(len(phase)):
            if phase[i][j] == 1:
                phase[i][j] = 0

    AP = abs(fftshift(fft2(fftshift(phase))))**2
    AP_max = np.nanmax(AP)
    AP = AP/AP_max
    return AP


def update_psf_data(pupil_grid, build='rebuild'):
    pupil_grid.update_data(build=build)
    ndim = pupil_grid.num_rays
    maxdim = pupil_grid.maxdim
    AP = calc_psf(pupil_grid.grid[2], ndim, maxdim)
    return AP
