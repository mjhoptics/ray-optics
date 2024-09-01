#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
""" functions to calculate etendue and imager parameters

.. Created on Tue Aug  6 18:08:16 2019

.. codeauthor: Michael J. Hayford
"""
import math

from rayoptics.util import dict2d
from rayoptics.util.dict2d import dict2D


obj_img_set = ['object', 'image']
fld_ape_set = ['field', 'aperture']
fld_labels = ['height', 'angle']
ap_labels = ['epd', 'NA', 'f/#']


def create_etendue_dict():
    """Returns an empty dict2D(fld_ape_set, obj_img_set)."""
    return dict2D(fld_ape_set, obj_img_set)


def na2slp(na: float, n=1.0) -> float:
    """ convert numerical aperture to slope """
    return na/n


def slp2na(slp: float, n=1.0) -> float:
    """ convert a ray slope to numerical aperture """
    return n*slp


def ang2slp(ang: float) -> float:
    """ convert an angle in degrees to a slope """
    return math.tan(math.radians(ang))


def slp2ang(slp: float) -> float:
    """ convert a slope to an angle in degrees """
    return math.degrees(math.atan(slp))


def get_aperture_from_slope(imager, slope, n=1):
    if slope == 0:
        return 0, 0, 0
    fno = -1/(2*slope)
    na = slp2na(slope, n=n)
    epd = imager.f/fno if imager.f is not None else None
    return na, fno, epd


def get_slope_from_aperture(imager, input_cell, n=1):
    if 'NA' in input_cell:
        na = input_cell['NA']
        slope = na2slp(na, n=n)
    elif 'f/#' in input_cell:
        fno = input_cell['f/#']
        slope = -1/(2*fno)
    elif 'epd' in input_cell:
        epd = input_cell['epd']
        if imager.f is None or imager.f == 0:
            slope = 0
        else:
            slope = -epd/(2*imager.f)
    else:
        slope = None
    return slope


def calc_aperture_from_input(conj_type, imager, input_cell, n=1):
    """ conj_type is for the io_input space"""
    slope = get_slope_from_aperture(imager, input_cell, n=n)
    if conj_type == 'finite':
        na, fno, epd = get_aperture_from_slope(imager, slope, n=n)
    else:
        
        na, fno, epd = 0, 0, input_cell['epd']
    return na, fno, epd


def do_etendue_via_imager(conj_type, imager, etendue_inputs, etendue_grid,
                          n_0=1, n_k=1):
    imager_inputs = {}
    li = dict2d.num_items_by_type(etendue_inputs, fld_ape_set, obj_img_set)
    if li['field'] == 1:
        row = dict2d.row(etendue_inputs, 'field')
        obj_img_key = 'object' if len(row['object']) else 'image'
        do_field_via_imager(conj_type, imager, etendue_inputs, obj_img_key,
                            etendue_grid, n_0=n_0, n_k=n_k)

    if li['aperture'] == 1:
        row = dict2d.row(etendue_inputs, 'aperture')
        obj_img_key = 'object' if len(row['object']) else 'image'
        do_aperture_via_imager(conj_type, imager, etendue_inputs, obj_img_key,
                               etendue_grid)

    if li['field'] == 2 or li['aperture'] == 2:
        fld_ape_key = 'field' if li['field'] == 2 else 'aperture'
        # solve for imager
        imager_inputs = do_etendue_to_imager(fld_ape_key, etendue_inputs,
                                             etendue_grid, n_0=n_0, n_k=n_k)

    return imager_inputs


def do_field_via_imager(conj_type, imager, etendue_inputs, obj_img_key,
                        etendue_grid, n_0=1, n_k=1):
    input_cell = etendue_inputs['field'][obj_img_key]
    if obj_img_key == 'object':
        output_cell = etendue_grid['field']['image']
        if 'angle' in input_cell:
            efl = imager.f
            obj_ang = input_cell['angle']
            obj_slp = ang2slp(obj_ang)
            output_cell['height'] = efl*obj_slp
        elif 'height' in input_cell:
            m = imager.m
            obj_ht = input_cell['height']
            output_cell['height'] = m*obj_ht

    elif obj_img_key == 'image':
        output_cell = etendue_grid['field']['object']
        if imager.m == 0:  # infinite conjugate
            efl = imager.f
            if 'height' in input_cell:
                img_ht = input_cell['height']
                obj_slp = img_ht/efl
                obj_ang = slp2ang(obj_slp)
                output_cell['angle'] = obj_ang
        else:  # finite conjugate
            m = imager.m
            if 'height' in input_cell:
                img_ht = input_cell['height']
                output_cell['height'] = img_ht/m


def do_aperture_via_imager(conj_type, imager, etendue_inputs, obj_img_key,
                           etendue_grid, n_0=1, n_k=1):
    inpt = etendue_inputs['aperture'][obj_img_key]
    input_cell = etendue_grid['aperture'][obj_img_key]
    if conj_type == 'infinite':
        efl = imager.f

        if obj_img_key == 'object':
            epd = inpt['epd']
            slpk = (epd/2.0)/efl
            na, fno, expd = get_aperture_from_slope(imager, slpk, n=n_k)
            output_cell = etendue_grid['aperture']['image']
            output_cell['f/#'] = fno
            output_cell['NA'] = na
            output_cell['epd'] = expd

        elif obj_img_key == 'image':
            slpk = get_slope_from_aperture(imager, inpt, n=n_k)
            na, fno, pupil = get_aperture_from_slope(imager, slpk, n=n_k)
            input_cell['f/#'] = fno
            input_cell['NA'] = na
            input_cell['epd'] = pupil

            epd = 2*slpk*efl

            output_cell = etendue_grid['aperture']['object']
            output_cell['epd'] = epd

    else:
        mag = imager.m

        if obj_img_key == 'object':
            slp0 = get_slope_from_aperture(imager, inpt, n=n_0)
            na, fno, pupil = get_aperture_from_slope(imager, slp0, n=n_0)
            input_cell['f/#'] = fno
            input_cell['NA'] = na
            input_cell['epd'] = pupil

            slpk = slp0/mag
            na, fno, pupil = get_aperture_from_slope(imager, slpk, n=n_k)
            output_cell = etendue_grid['aperture']['image']

        elif obj_img_key == 'image':
            slpk = get_slope_from_aperture(imager, inpt, n=n_k)
            na, fno, pupil = get_aperture_from_slope(imager, slpk, n=n_k)
            input_cell['f/#'] = fno
            input_cell['NA'] = na
            input_cell['epd'] = pupil

            slp0 = mag*slpk
            na, fno, pupil = get_aperture_from_slope(imager, slp0, n=n_0)
            output_cell = etendue_grid['aperture']['object']

        output_cell['f/#'] = fno
        output_cell['NA'] = na
        output_cell['epd'] = pupil


def do_etendue_to_imager(fld_ape_key, etendue_inputs, etendue_grid,
                         n_0=1.0, n_k=1.0):
    imager_inputs = None
    if fld_ape_key == 'field':
        obj_cell = etendue_inputs['field']['object']
        img_cell = etendue_inputs['field']['image']
        obj_grid = etendue_grid['field']['object']
        img_grid = etendue_grid['field']['image']
        obj_key = obj_cell.keys()
        img_key = img_cell.keys()
        if 'angle' in obj_key:
            obj_ang = obj_cell['angle']
            obj_slp = ang2slp(obj_ang)
            if 'height' in img_key:
                img_ht = img_cell['height']
                img_grid['height'] = img_ht
                efl = img_ht/obj_slp
                imager_inputs = 'f', efl
            else:
                return None
        elif 'height' in obj_key:
            obj_ht = obj_cell['height']
            if 'height' in img_key:
                img_ht = img_cell['height']
                img_grid['height'] = img_ht
                mag = img_ht/obj_ht
                imager_inputs = 'm', mag

    if fld_ape_key == 'aperture':
        obj_cell = etendue_inputs['aperture']['object']
        img_cell = etendue_inputs['aperture']['image']
        obj_grid = etendue_grid['aperture']['object']
        img_grid = etendue_grid['aperture']['image']
        obj_key = obj_cell.keys()
        img_key = img_cell.keys()
        if 'epd' in obj_key:
            epd = obj_cell['epd']
            if 'f/#' in img_key:
                fno = img_cell['f/#']
                img_grid['f/#'] = fno
                efl = epd * fno
                imager_inputs = 'f', efl
            if 'NA' in img_key:
                na = img_cell['NA']
                img_grid['NA'] = na
                slpk = na2slp(na, n=n_k)
                efl = (epd/2.0)/slpk
                imager_inputs = 'f', efl
        else:
            if 'NA' in obj_key:
                nao = obj_cell['NA']
                obj_grid['NA'] = nao
                slp0 = na2slp(nao, n=n_0)
            elif 'f/#' in obj_key:
                fno = obj_cell['f/#']
                obj_grid['f/#'] = fno
                slp0 = -1/(2*fno)

            if 'NA' in img_key:
                na = img_cell['NA']
                img_grid['NA'] = na
                slpk = na2slp(na, n=n_k)
            elif 'f/#' in img_key:
                fno = img_cell['f/#']
                img_grid['f/#'] = fno
                slpk = -1/(2*fno)
            mag = slp0/slpk
            imager_inputs = 'm', mag

    return imager_inputs


def fill_in_etendue_data(conj_type, imager, fld_ape_key,
                         inputs, values, n=1.0):
    if len(inputs) == 0:
        return

    if fld_ape_key == 'field':
        for key in inputs:
            values[key] = inputs[key]

    if fld_ape_key == 'aperture':
        na, fno, pupil = calc_aperture_from_input(conj_type, imager,
                                                  inputs, n=n)
        values['NA'] = na
        values['f/#'] = fno
        values['epd'] = pupil
