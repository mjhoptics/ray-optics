#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Wireframe rendering routines for 2D lens picture

.. Created on Tue Sep 18 14:23:28 2018

.. codeauthor: Michael J. Hayford
"""

import numpy as np

import rayoptics.optical.elements as ele
import rayoptics.optical.gap as gap
import rayoptics.optical.transform as trns
from rayoptics.optical.trace import trace_boundary_rays_at_field
import rayoptics.optical.model_constants as mc


def shift_start_of_ray_bundle(ray_bundle, start_offset, r, t):
    """ modify ray_bundle so that rays begin "start_offset" from 1st surface

    Args:
        ray_bundle: list of rays in a bundle, i.e. all for one field.
                    ray_bundle[0] is assumed to be the chief ray
        start_offset: z distance rays should start wrt first surface.
                      positive if to left of first surface
        r: transformation rotation
        t: transformation translation
    """

    for ri, ray in enumerate(ray_bundle):
        b4_pt = r.dot(ray[mc.ray][1][mc.p]) - t
        b4_dir = r.dot(ray[mc.ray][0][mc.d])
        if ri == 0:
            # For the chief ray, use the input offset.
            dst = -b4_pt[2]/b4_dir[2]
        else:
            pt0 = ray_bundle[0][mc.ray][0][mc.p]
            dir0 = ray_bundle[0][mc.ray][0][mc.d]
            # Calculate distance along ray to plane perpendicular to
            #  the chief ray.
            dst = -(b4_pt - pt0).dot(dir0)/b4_dir.dot(dir0)
        pt = b4_pt + dst*b4_dir
        ray[mc.ray][0][mc.p] = pt
        ray[mc.ray][0][mc.d] = b4_dir


def setup_shift_of_ray_bundle(seq_model, start_offset):
    """ compute transformation for rays "start_offset" from 1st surface

    Args:
        seq_model: the sequential model
        start_offset: z distance rays should start wrt first surface.
                      positive if to left of first surface
    Returns:
        transformation rotation and translation::
            (r, t)
    """

    s1 = seq_model.ifcs[1]
    s0 = seq_model.ifcs[0]
    g0 = gap.Gap(start_offset, seq_model.gaps[0].medium)
    r, t = trns.reverse_transform(s0, g0, s1)
    return r, t


def shift_start_of_rayset(opt_model, rayset, start_offset):
    """ modify rayset so that rays begin "start_offset" from 1st surface

    Args:
        opt_model: the optical model
        rayset: list of ray bundles, i.e. for a list of fields
        start_offset: z distance rays should start wrt first surface.
                      positive if to left of first surface
    Returns:
        transformation rotation and translation
    """

    r, t = setup_shift_of_ray_bundle(opt_model.seq_model, start_offset)
    for ray_bundle in rayset:
        shift_start_of_ray_bundle(ray_bundle, start_offset, r, t)
    return r, t


def bbox_from_poly(poly):
    minx, miny = np.min(poly, axis=0)
    maxx, maxy = np.max(poly, axis=0)
    return np.array([[minx, miny], [maxx, maxy]])


def scale_bounds(bbox, oversize_factor):
    inc_x = oversize_factor*(bbox[1][0] - bbox[0][0])
    inc_y = oversize_factor*(bbox[1][1] - bbox[0][1])
    incr = max(inc_x, inc_y)
    return np.array([[bbox[0][0]-incr, bbox[0][1]-incr],
                     [bbox[1][0]+incr, bbox[1][1]+incr]])


def transform_poly(tfrm, poly):
        coord_flip = np.array([[0., 1.], [1., 0.]])

        poly = np.matmul(coord_flip, poly.T)
        poly = np.matmul(tfrm[0][1:, 1:], poly).T

        t = np.array([tfrm[1][1], tfrm[1][2]])
        poly += t

        # flip coordinates back to 2D plot coordinates, +y points up
        poly = np.matmul(poly, coord_flip)
        bbox = bbox_from_poly(poly)
        return poly, bbox


def create_optical_element(e):
    if isinstance(e, ele.Element):
        return LensElement(e)
    elif isinstance(e, ele.Mirror):
        return MirrorElement(e)
    else:
        return OpticalElement(e)


class OpticalElement():
    def __init__(self, e):
        self.e = e

    def update_shape(self):
        poly = np.array(self.e.shape())
        return transform_poly(self.e.tfrm, poly)

    def render_color(self):
        return self.e.render_color

    def update_element_shape(self, view):
        poly, bbox = self.update_shape()
        p = view.create_polygon(poly, self.render_color())
        return p, bbox


class LensElement(OpticalElement):
    def __init__(self, e):
        super().__init__(e)
        self.handles = {}

    def update_shape(self):
        self.e.render_handles()
        for key, value in self.e.handles.items():
            poly = np.array(value[0])
            self.handles[key] = transform_poly(self.e.tfrm, poly)
        return self.handles['shape']

    def update_element_shape(self, view):
        poly, bbox = self.update_shape()
        p = view.create_polygon(poly, self.render_color())
        return p, bbox

    def render_color(self):
        return self.e.render_color


class MirrorElement(OpticalElement):
    def __init__(self, e):
        super().__init__(e)

    def update_shape(self):
        poly = np.array(self.e.shape())
        return transform_poly(self.e.tfrm, poly)

    def render_color(self):
        return self.e.render_color


class RayBundle():
    """ class for ray bundle from a single field point """
    def __init__(self, opt_model, fld, wvl, start_offset):
        self.opt_model = opt_model
        self.fld = fld
        self.wvl = wvl
        self.start_offset = start_offset

    def update_shape(self):
        try:
            rayset = trace_boundary_rays_at_field(self.opt_model,
                                                  self.fld, self.wvl)
        except:
            wvl = self.opt_model.optical_spec.spectral_region.central_wvl
            self.wvl = wvl
            rayset = trace_boundary_rays_at_field(self.opt_model,
                                                  self.fld, wvl)
        finally:
            # If the object distance (tfrms[0][1][2]) is greater than the
            #  start_offset, then modify rayset start to match start_offset.
            # Remember object transformation for resetting at the end.
            seq_model = self.opt_model.seq_model
            tfrms = seq_model.gbl_tfrms
            tfrtm0 = tfrms[0]

            try:
                if abs(tfrms[0][1][2]) > self.start_offset:
                    r, t = setup_shift_of_ray_bundle(seq_model,
                                                     self.start_offset)
                    tfrms[0] = (r, t)
                    shift_start_of_ray_bundle(rayset, self.start_offset,
                                              r, t)

                poly1 = []
                for i, r in enumerate(rayset[3][mc.ray][0:]):
                    rot, trns = tfrms[i]
                    p = rot.dot(r[mc.p]) + trns
        #            print(i, r[0], rot, trns, p)
        #            print("r3", i, p[2], p[1])
                    poly1.append([p[2], p[1]])

                poly2 = []
                for i, r in enumerate(rayset[4][mc.ray][0:]):
                    rot, trns = tfrms[i]
                    p = rot.dot(r[mc.p]) + trns
        #            print(i, r[0], rot, trns, p)
        #            print("r4", i, p[2], p[1])
                    poly2.append([p[2], p[1]])

                poly2.reverse()
                poly1.extend(poly2)
                bbox = bbox_from_poly(poly1)

            finally:
                tfrms[0] = tfrtm0

        return poly1, bbox

    def update_ray_fan_shape(self, view):
        rndr_clr = [254, 197, 254, 64]  # magenta, 25%

        poly, bbox = self.update_shape()
        p = view.create_polygon(poly, rndr_clr)

        return p, bbox


# TODO (mjhoptics@gmail.com): additional refactoring needed in clients
class LensLayout():
    def __init__(self, opt_model, **kwargs):
        self.opt_model = opt_model

        ele_model = self.opt_model.ele_model
        if len(ele_model.elements) == 0:
            ele_model.elements_from_sequence(self.opt_model.seq_model)

    def system_length(self, ele_bbox):
        osp = self.opt_model.optical_spec
        img_dist = abs(osp.parax_data[2].img_dist)
        ele_length = ele_bbox[1][0] - ele_bbox[0][0]
        return ele_length+img_dist

    def create_element_model(self, view):
        elements = []
        ele_model = self.opt_model.ele_model
        for e in ele_model.elements:
            oe = create_optical_element(e)
            elements.append((oe.update_element_shape, view))
        return elements

    def create_ray_model(self, view, start_offset):
        ray_bundles = []
        fov = self.opt_model.optical_spec.field_of_view
        wvl = self.opt_model.seq_model.central_wavelength()
        for fld in fov.fields:
            rb = RayBundle(self.opt_model, fld, wvl, start_offset)
            ray_bundles.append((rb.update_ray_fan_shape, view))
        return ray_bundles
