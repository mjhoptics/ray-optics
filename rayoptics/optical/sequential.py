#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Manager class for a sequential optical model

.. codeauthor: Michael J. Hayford
"""

import itertools
import logging

from . import surface
from . import profiles
from . import gap
from . import medium as m
from . import raytrace as rt
from . import trace as trace
from . import transform as trns
from rayoptics.optical.model_constants import Surf, Gap
from opticalglass import glassfactory as gfact
from opticalglass import glasserror as ge
import numpy as np
import pandas as pd
from math import copysign
from rayoptics.util.misc_math import isanumber


class SequentialModel:
    """ Manager class for a sequential optical model

    A sequential optical model is a sequence of surfaces and gaps.

    The sequential model has this structure
    ::

        IfcObj  Ifc1  Ifc2  Ifc3 ... Ifci-1   IfcImg
             \  /  \  /  \  /             \   /
             GObj   G1    G2              Gi-1

    where

        - Ifc is a :class:`~rayoptics.optical.surface.Interface` instance
        - G   is a :class:`~rayoptics.optical.gap.Gap` instance

    There are N interfaces and N-1 gaps. The initial configuration has an
    object and image Surface and an object gap.

    The Interface API supports implementation of an optical action, such as
    refraction, reflection, scatter, diffraction, etc. The Interface may be
    realized as a physical profile separating the adjacent gaps or an idealized
    object, such as a thin lens or 2 point HOE.
    The Gap class maintains a simple separation (z translation) and the medium
    filling the gap. More complex coordinate transformations are handled
    through the Interface API.
    """

    def __init__(self, opt_model):
        self.opt_model = opt_model
        self.ifcs = []
        self.gaps = []
        self.transforms = []
        self.lcl_tfrms = []
        self.z_dir = []
        self.stop_surface = 1
        self.cur_surface = 0
        self.rndx = []
        self.ifcs.append(surface.Surface('Obj'))
        self.ifcs.append(surface.Surface('Img'))
        self.gaps.append(gap.Gap())

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['opt_model']
        del attrs['transforms']
        del attrs['lcl_tfrms']
        del attrs['z_dir']
        del attrs['rndx']
        return attrs

    def reset(self):
        self.__init__()

    def get_num_surfaces(self):
        return len(self.ifcs)

    def path(self):
        return itertools.zip_longest(self.ifcs, self.gaps)

    def calc_ref_indices_for_spectrum(self, wvls):
        """ returns a |DataFrame| with refractive indices for all **wvls**

        Args:
            wvls: list of wavelengths in nm
        """
        indices = []
        for g in self.gaps:
            ri = []
            mat = g.medium
            for w in wvls:
                ri.append(mat.rindex(w))
            indices.append(ri)
        indices.append(indices[-1])

        return pd.DataFrame(indices, columns=wvls)

    def central_wavelength(self):
        """ returns the central wavelength in nm of the model's ``WvlSpec`` """
        return self.opt_model.optical_spec.spectral_region.central_wvl

    def central_rndx(self, i):
        """ returns the central refractive index of the model's ``WvlSpec`` """
        central_wvl = self.central_wavelength()
        return self.rndx[central_wvl].iloc[i]

    def get_surface_and_gap(self, srf=None):
        if not srf:
            srf = self.cur_surface
        s = self.ifcs[srf]
        if srf == len(self.gaps):
            g = None
        else:
            g = self.gaps[srf]
        return s, g

    def set_cur_surface(self, s):
        self.cur_surface = s

    def set_stop(self):
        """ sets the stop surface to the current surface """
        self.stop_surface = self.cur_surface
        return self.stop_surface

    def __iadd__(self, node):
        if isinstance(node, gap.Gap):
            self.gaps.append(node)
        else:
            self.ifcs.insert(len(self.ifcs)-1, node)
        return self

    def insert(self, surf, gap):
        """ insert surf and gap at the cur_gap edge of the sequential model
            graph """
        self.cur_surface += 1
        self.ifcs.insert(self.cur_surface, surf)
        self.gaps.insert(self.cur_surface, gap)

    def add_surface(self, surf):
        """ add a surface where surf is a list that contains:
            [curvature, thickness, refractive_index, v-number] """
        s, g = self.create_surface_and_gap(surf)
        self.insert(s, g)

    def update_surface(self, surf):
        """ add a surface where surf is a list that contains:
            [curvature, thickness, refractive_index, v-number] """

        self.insert(*self.create_surface_and_gap(surf))

    def sync_to_restore(self, opt_model):
        self.opt_model = opt_model
        if hasattr(self, 'optical_spec'):
            opt_model.optical_spec = self.optical_spec
            delattr(self, 'optical_spec')
        for sg in self.path():
            if hasattr(sg[Surf], 'sync_to_restore'):
                sg[Surf].sync_to_restore(self)
            if sg[Gap]:
                if hasattr(sg[Gap], 'sync_to_restore'):
                    sg[Gap].sync_to_restore(self)

    def update_model(self):
        # delta n across each surface interface must be set to some
        #  reasonable default value. use the index at the central wavelength
        osp = self.opt_model.optical_spec
        ref_wl = osp.spectral_region.reference_wvl

        wvlns = osp.spectral_region.wavelengths
        self.rndx = self.calc_ref_indices_for_spectrum(wvlns)
        n_before = self.rndx.iloc[0]

        self.z_dir = []
        z_dir_before = 1.0

        for i, s in enumerate(self.ifcs):
            z_dir_after = copysign(1.0, z_dir_before)
            n_after = np.copysign(self.rndx.iloc[i], n_before)
            if s.refract_mode == 'REFL':
                n_after = -n_after
                z_dir_after = -z_dir_after

            s.delta_n = n_after.iloc[ref_wl] - n_before.iloc[ref_wl]
            n_before = n_after
            self.rndx.iloc[i] = n_after
            z_dir_before = z_dir_after
            self.z_dir.append(z_dir_after)
            # call update() on the surface interface
            s.update()

        self.transforms = self.compute_global_coords()
        self.lcl_tfrms = self.compute_local_transforms()
        osp.update_model()

        self.set_clear_apertures()

    def insert_surface_and_gap(self):
        s = surface.Surface()
        g = gap.Gap()
        self.insert(s, g)
        return s, g

    def update_surface_and_gap_from_cv_input(self, dlist, idx=None):
        if isinstance(idx, int):
            s, g = self.get_surface_and_gap(idx)
        else:
            s, g = self.insert_surface_and_gap()

        if self.opt_model.radius_mode:
            if dlist[0] != 0.0:
                s.profile.cv = 1.0/dlist[0]
            else:
                s.profile.cv = 0.0
        else:
            s.profile.cv = dlist[0]

        if g:
            g.thi = dlist[1]
            if len(dlist) < 3:
                g.medium = m.Air()
            elif dlist[2].upper() == 'REFL':
                s.refract_mode = 'REFL'
                g.medium = self.gaps[self.cur_surface-1].medium
            elif isanumber(dlist[2]):
                n, v = m.glass_decode(float(dlist[2]))
                g.medium = m.Glass(n, v, '')
            else:
                name, cat = dlist[2].split('_')
                if cat.upper() == 'SCHOTT' and name[:1].upper() == 'N':
                    name = name[:1]+'-'+name[1:]
                try:
                    g.medium = gfact.create_glass(name, cat)
                except ge.GlassNotFoundError as gerr:
                    logging.info('%s glass data type %s not found',
                                 gerr.catalog,
                                 gerr.name)
                    logging.info('Replacing material with air.')
                    g.medium = m.Air()
        else:
            # at image surface, apply defocus to previous thickness
            self.gaps[idx-1].thi += dlist[1]

    def update_surface_profile_cv_input(self, profile_type, idx=None):
        if not isinstance(idx, int):
            idx = self.cur_surface
        cur_profile = self.ifcs[idx].profile
        new_profile = profiles.mutate_profile(cur_profile, profile_type)
        self.ifcs[idx].profile = new_profile
        return self.ifcs[idx].profile

    def update_surface_decenter_cv_input(self, idx=None):
        if not isinstance(idx, int):
            idx = self.cur_surface
        if not self.ifcs[idx].decenter:
            self.ifcs[idx].decenter = surface.DecenterData()
        return self.ifcs[idx].decenter

    def create_surface_and_gap(self, surf_data):
        """ create a surface and gap where surf_data is a list that contains:
            [curvature, thickness, refractive_index, v-number] """
        s = surface.Surface()

        if self.opt_model.radius_mode:
            if surf_data[0] != 0.0:
                s.profile.cv = 1.0/surf_data[0]
            else:
                s.profile.cv = 0.0
        else:
            s.profile.cv = surf_data[0]

        if len(surf_data) > 2:
            if isanumber(surf_data[2]):
                if len(surf_data) < 3:
                    if surf_data[2] == 1.0:
                        mat = m.Air()
                    else:
                        mat = m.Medium(surf_data[2])
                else:
                    mat = m.Glass(surf_data[2], surf_data[3], '')

            else:
                if surf_data[2].upper() == 'REFL':
                    s.refract_mode = 'REFL'
                    mat = self.gaps[self.cur_surface-1].medium
                else:
                    name, cat = surf_data[2], surf_data[3]
                    try:
                        mat = gfact.create_glass(name, cat)
                    except ge.GlassNotFoundError as gerr:
                        logging.info('%s glass data type %s not found',
                                     gerr.catalog,
                                     gerr.name)
                        logging.info('Replacing material with air.')
                        mat = m.Air()

        else:  # only curvature and thickness entered, set material to air
            mat = m.Air()

        g = gap.Gap(surf_data[1], mat)
        return s, g

    def surface_label_list(self):
        """ list of surface labels or surface number, if no label """
        labels = []
        for i, s in enumerate(self.ifcs):
            if len(s.label) == 0:
                if i == self.stop_surface:
                    labels.append('Stop')
                else:
                    labels.append(str(i))
            else:
                labels.append(s.label)
        return labels

    def list_model(self):
        for i, sg in enumerate(self.path()):
            if sg[Gap]:
                print(i, sg[Surf])
                print('    ', sg[Gap])
            else:
                print(i, sg[Surf])

    def list_gaps(self):
        for i, gp in enumerate(self.gaps):
            print(i, gp)

    def list_surfaces(self):
        for i, s in enumerate(self.ifcs):
            print(i, s)

    def list_surface_and_gap(self, i):
        s = self.ifcs[i]
        cvr = s.profile.cv
        if self.opt_model.radius_mode:
            if cvr != 0.0:
                cvr = 1.0/cvr
        sd = s.surface_od()

        if i < len(self.gaps):
            g = self.gaps[i]
            thi = g.thi
            med = g.medium.name()
        else:
            thi = ''
            med = ''
        return [cvr, thi, med, sd]

    def list_decenters(self):
        for i, sg in enumerate(self.path()):
            if sg[Gap]:
                print(i, sg[Gap])
                if sg[Surf].decenter is not None:
                    print(' ', repr(sg[Surf].decenter))
            else:
                if sg[Surf].decenter is not None:
                    print(i, repr(sg[Surf].decenter))

    def list_elements(self):
        for i, gp in enumerate(self.gaps):
            if gp.medium.label.lower() != 'air':
                print(self.ifcs[i].profile,
                      self.ifcs[i+1].profile,
                      gp)

    def trace_fan(self, fct, fi, xy, num_rays=21):
        """ xy determines whether x (=0) or y (=1) fan """
        osp = self.opt_model.optical_spec
        fld = osp.field_of_view.fields[fi]
        wvl = self.central_wavelength()
        foc = 0.0
#        trace.setup_canonical_coords(self, fld, wvl)
        rs_pkg, cr_pkg = trace.setup_pupil_coords(self.opt_model,
                                                  fld, wvl, foc)
        fld.chief_ray = cr_pkg
        fld.ref_sphere = rs_pkg
        img_pt = cr_pkg[0].ray[-1][0]

        wvls = osp.spectral_region
        fans_x = []
        fans_y = []
        fan_start = np.array([0., 0.])
        fan_stop = np.array([0., 0.])
        fan_start[xy] = -1.0
        fan_stop[xy] = 1.0
        fan_def = [fan_start, fan_stop, num_rays]
        max_y_val = 0.0
        rc = []
        for wi, wvl in enumerate(wvls.wavelengths):
            rc.append(wvls.render_colors[wi])

            rs_pkg, cr_pkg = trace.setup_pupil_coords(self.opt_model,
                                                      fld, wvl, foc,
                                                      image_pt=img_pt)
            fld.chief_ray = cr_pkg
            fld.ref_sphere = rs_pkg
            fan = trace.trace_fan(self.opt_model, fan_def, fld, wvl, foc,
                                  img_filter=lambda p, ray_pkg:
                                  fct(p, xy, ray_pkg, fld, wvl))
            f_x = []
            f_y = []
            for p, y_val in fan:
                f_x.append(p[xy])
                f_y.append(y_val)
                if abs(y_val) > max_y_val:
                    max_y_val = abs(y_val)
            fans_x.append(f_x)
            fans_y.append(f_y)
        fans_x = np.array(fans_x)
        fans_y = np.array(fans_y)
        return fans_x, fans_y, max_y_val, rc

    def trace_grid(self, fct, fi, wl=None, num_rays=21, form='grid',
                   append_if_none=True):
        """ fct is applied to the raw grid and returned as a grid  """
        osp = self.opt_model.optical_spec
        wvls = osp.spectral_region
        wvl = self.central_wavelength()
        wv_list = wvls.wavelengths if wl is None else [wvl]
        fld = osp.field_of_view.fields[fi]
        foc = 0.0

        rs_pkg, cr_pkg = trace.setup_pupil_coords(self.opt_model,
                                                  fld, wvl, foc)
        fld.chief_ray = cr_pkg
        fld.ref_sphere = rs_pkg

        grids = []
        origin = -.05
        delta = 0.1
        grid_start = np.array([origin, origin])
        grid_stop = np.array([origin+delta, origin+delta])
        grid_def = [grid_start, grid_stop, num_rays]
        for wi, wvl in enumerate(wv_list):
            grid = trace.trace_grid(self.opt_model, grid_def, fld, wvl, foc,
                                    form=form, append_if_none=append_if_none,
                                    img_filter=lambda p, ray_pkg:
                                    fct(p, wi, ray_pkg, fld, wvl))
            grids.append(grid)
        rc = wvls.render_colors
        return grids, rc

    def trace_wavefront(self, fld, wvl, foc, num_rays=32):

        def wave(p, ray_pkg, fld, wvl):
            x = p[0]
            y = p[1]
            if ray_pkg is not None:
                opd_pkg = rt.wave_abr(fld, wvl, ray_pkg)
                opd = opd_pkg[0]/self.opt_model.nm_to_sys_units(wvl)
            else:
                opd = 0.0
            return np.array([x, y, opd])

        rs_pkg, cr_pkg = trace.setup_pupil_coords(self.opt_model,
                                                  fld, wvl, foc)
        fld.chief_ray = cr_pkg
        fld.ref_sphere = rs_pkg

        grid_start = np.array([-1., -1.])
        grid_stop = np.array([1., 1.])
        grid_def = (grid_start, grid_stop, num_rays)

        grid = trace.trace_grid(self.opt_model, grid_def, fld, wvl, foc,
                                img_filter=lambda p, ray_pkg:
                                wave(p, ray_pkg, fld, wvl), form='grid')
        return grid

    def set_clear_apertures(self):
        def rd(v):
            """ take 2d length of input vector v """
            return np.sqrt(v[0]*v[0]+v[1]*v[1])

        if self.get_num_surfaces() > 2:
            fields_df = trace.trace_all_fields(self.opt_model)
            # a) Select the inc_pt data from the unstacked result and transpose
            #    so that intrfcs are the index
            inc_pts = fields_df.unstack()['inc_pt'].T
            # b) applymap() is used to apply the function rd() to each element
            #    in the dataframe
            inc_pts_rd = inc_pts.applymap(rd)
            # c) apply max() function to each row (i.e. across columns, axis=1)
            semi_ap = inc_pts_rd.max(axis=1)
            for s, max_ap in zip(self.ifcs[1:-1], semi_ap[1:-1]):
                s.set_max_aperture(max_ap)

    def trace(self, pt0, dir0, wvl, **kwargs):
        return rt.trace(self, pt0, dir0, wvl, **kwargs)

    def compute_global_coords(self, glo=1):
        """ Return global surface coordinates (rot, t) wrt surface glo. """
        tfrms = []
        r, t = np.identity(3), np.array([0., 0., 0.])
        prev = r, t
        tfrms.append(prev)
#        print(glo, t, *np.rad2deg(t3d.euler.mat2euler(r)))
        if glo > 0:
            # iterate in reverse over the segments before the
            #  global reference surface
            go = glo
            path = itertools.zip_longest(self.ifcs[glo::-1],
                                         self.gaps[glo-1::-1])
            after = next(path)
            # loop of remaining surfaces in path
            while True:
                try:
                    before = next(path)
                    go -= 1
                    r, t = trns.reverse_transform(before[Surf], after[Gap],
                                                  after[Surf])
                    t = prev[0].dot(t) + prev[1]
                    r = prev[0].dot(r)
#                    print(go, t,
#                          *np.rad2deg(euler2opt(t3d.euler.mat2euler(r))))
                    prev = r, t
                    tfrms.append(prev)
                    after = before
                except StopIteration:
                    break
            tfrms.reverse()
        path = itertools.zip_longest(self.ifcs[glo:], self.gaps[glo:])
        before = next(path)
        prev = np.identity(3), np.array([0., 0., 0.])
        go = glo
        # loop forward over the remaining surfaces in path
        while True:
            try:
                after = next(path)
                go += 1
                r, t = trns.forward_transform(before[Surf], before[Gap],
                                              after[Surf])
                t = prev[0].dot(t) + prev[1]
                r = prev[0].dot(r)
#                print(go, t,
#                      *np.rad2deg(euler2opt(t3d.euler.mat2euler(r))))
                prev = r, t
                tfrms.append(prev)
                before = after
            except StopIteration:
                break
        return tfrms

    def compute_local_transforms(self):
        """ Return forward surface coordinates (r.T, t) for each interface. """
        tfrms = []
        path = self.path()
        before = next(path)
        while before is not None:
            try:
                after = next(path)
            except StopIteration:
                tfrms.append((np.identity(3), np.array([0., 0., 0.])))
                break
            else:
                r, t = trns.forward_transform(before[Surf], before[Gap],
                                              after[Surf])
                rt = r.transpose()
                tfrms.append((rt, t))
                before = after

        return tfrms
