#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
"""

.. Created on Tue Aug  6 18:08:16 2019

.. codeauthor: Michael J. Hayford
"""
import math
from collections import namedtuple

from rayoptics.util import dict2d
from rayoptics.util.dict2d import dict2D


obj_img_set = ['object', 'image']
fld_ape_set = ['field', 'aperture']
fld_labels = ['height', 'angle']
ap_labels = ['pupil', 'NA', 'f/#']


def na2slp(na, n=1.0):
    """ convert numerical aperture to slope """
    return n*math.tan(math.asin(na/n))


def slp2na(slp, n=1.0):
    """ convert a ray slope to numerical aperture """
    return n*math.sin(math.atan(slp/n))


def ang2slp(ang):
    """ convert an angle in degrees to a slope """
    return math.tan(math.radians(ang))


def slp2ang(slp):
    """ convert a slope to an angle in degrees """
    return math.degrees(math.atan(slp))


def etendue_setup(**inputs):
    """ Calculate the ideal imaging properties given two independent parameters

    Given 2 system parameters from the following list, this function
    calculates the remaining parameters.

    Note that if specifying ``tt`` and ``f``, their ratio, tt/f, must be
    greater than or equal to 4. A `ValueError` is raised otherwise.

    For a typical system, the value of ``s`` is negative, i.e. the object is to
    the left of the first principal plane.

    Example::

        In [3]: m1s1 = etendue_setup(m=-0.5, s=-10.0); m1s1
        Out[3]: IdealImager(m=-0.5, s=-10.0, sp=5.0, tt=15.0, f=3.333333333333)

        In [4]: s_inf_efl = etendue_setup(s=-math.inf, f=25.0); s_inf_efl
        Out[4]: IdealImager(m=-0.0, s=-inf, sp=25.0, tt=inf, f=25.0)

    Args:
        m: (lateral) magnification
        s: object distance from first principal plane, P1->Obj
        sp: image distance from second principal plane, P2->Img
        tt: total track length, tt = sp - s
        f: focal length

    Returns:
        :class:`IdealImager` namedtuple

    Raises:
        ValueError: if tt/f < 4
    """
    if 'm' in inputs:
        m = inputs['m']
        if 's' in inputs:
            s = inputs['s']
            sp = m*s
            tt = sp - s
            f = s*sp/(s - sp)
        elif 'sp' in inputs:
            sp = inputs['sp']
            s = sp/m
            tt = sp - s
            f = s*sp/(s - sp)
        elif 'tt' in inputs:
            tt = inputs['tt']
            s = tt/(m - 1)
            sp = m*s
            f = s*sp/(s - sp)
        elif 'f' in inputs:
            f = inputs['f']
            tt = -f*(m - 1)**2/m
            s = tt/(m - 1)
            sp = m*s
        else:
            return None

    elif 's' in inputs:
        # arrange calculations so that s=-inf is handled gracefully
        s = inputs['s']
        if 'sp' in inputs:
            sp = inputs['sp']
            f = 1/(1/sp - 1/s)
            m = sp/s
            tt = sp - s
        elif 'tt' in inputs:
            tt = inputs['tt']
            m = 1 + tt/s
            sp = m*s
            f = s*sp/(s - sp)
        elif 'f' in inputs:
            f = inputs['f']
            m = f/(s + f)
            sp = 1/(1/f + 1/s)
            tt = sp - s
        else:
            return None

    elif 'sp' in inputs:
        # arrange calculations so that sp=inf is handled gracefully
        sp = inputs['sp']
        if 'tt' in inputs:
            tt = inputs['tt']
            m = sp/(sp - tt)
            s = sp/m
            f = s*sp/(s - sp)
        elif 'f' in inputs:
            f = inputs['f']
            m = (f - sp)/f
            s = 1/(1/sp - 1/f)
            tt = sp - s
        else:
            return None

    elif 'tt' in inputs:
        tt = inputs['tt']
        if 'f' in inputs:
            f = inputs['f']
            ttf = tt/f
            # tt/f >= 4, else no solution
            # pick root (+) that gives |s|>=|sp|, i.e. -1 <= m < 0
            m = ((2 - ttf) + math.sqrt(ttf*(ttf - 4)))/2
            s = tt/(m - 1)
            sp = m*s
        else:
            return None

    else:
        return None

    return IdealImager(m, s, sp, tt, f)


item_labels = ['height', 'angle', 'pupil', 'NA', 'f/#']

               
def do_aperture_spec_to_imager(**inputs):
    n_0 = 1.0
    n_k = 1.0
    if 'pupil' in inputs:
        epd = inputs['pupil']
        if 'f/#' in inputs:
            fno = inputs['f/#']
            efl = epd * fno
        if 'NA' in inputs:
            na = inputs['NA']
            slpk = na2slp(na, n=n_k)
            efl = (epd/2.0)/slpk
        if 'nao' in inputs:
            nao = inputs['nao']
            slp0 = n_0*math.tan(math.asin(nao/n_0))
            if 'f/#' in inputs:
                fno = inputs['f/#']
                efl = epd * fno
            if 'na' in inputs:
                na = inputs['na']
                slpk = n_k*math.tan(math.asin(na/n_k))
                efl = (epd/2.0)/slpk
    if 'f/#' in inputs:
        fno = inputs['f/#']
        if 'pupil' in inputs:
            epd = inputs['pupil']
            efl = epd * fno
    if 'NA' in inputs:
        na = inputs['NA']
        slpk = na2slp(na, n=n_k)
        if 'pupil' in inputs:
            epd = inputs['pupil']
            efl = (epd/2.0)/slpk
    return efl


def do_etendue_via_imager(imager_inputs, imager, etendue_inputs, etendue_grid):
    li = dict2d.num_items_by_type(etendue_inputs, fld_ape_set, obj_img_set)
    if li['field'] == 1:
        obj_cell = etendue_inputs['field']['object']
        img_cell = etendue_inputs['field']['image']
        obj_grid = etendue_grid['field']['object']
        img_grid = etendue_grid['field']['image']
        obj_key = obj_cell.keys()
        img_key = img_cell.keys()
        if len(obj_key) > 0:
            if 'angle' in obj_key:
                efl = imager.f
                obj_ang = obj_cell['angle']
                obj_slp = ang2slp(obj_ang)
                if 'height' in img_grid:
                    img_ht = efl*obj_slp
                    img_grid['height'] = img_ht
                else:
                    return None
            elif 'height' in obj_key:
                m = imager.m
                obj_ht = obj_cell['height']
                if 'height' in img_key:
                    img_ht = m*obj_ht
                    img_cell['height'] = img_ht
        elif len(img_key) > 0:
            if 'height' in img_key:
                img_ht = img_cell['height']
                obj_ang = obj_cell['angle']
                obj_slp = ang2slp(obj_ang)
                if 'height' in obj_grid:
                    m = imager.m
                    obj_ht = m/img_ht
                    obj_grid['height'] = obj_ht
                elif 'angle' in obj_grid:
                    efl = imager.f
                    obj_slp = img_ht/efl
                    obj_grid['angle'] = slp2ang(obj_slp)

    elif li['field'] == 2:
        obj_cell = etendue_inputs['field']['object']
        img_cell = etendue_inputs['field']['image']
        obj_key = obj_cell.keys()
        img_key = img_cell.keys()
        if 'angle' in obj_key:
            obj_ang = obj_cell['angle']
            obj_slp = ang2slp(obj_ang)
            if 'height' in img_key:
                img_ht = img_cell['height']
                efl = img_ht/obj_slp
                imager_inputs = 'f', efl
            else:
                return None
        elif 'height' in obj_key:
            obj_ht = obj_cell['height']
            if 'height' in img_key:
                img_ht = img_cell['height']
                mag = img_ht/obj_ht
                imager_inputs = 'm', mag

    elif li['aperture'] >= 1:
        obj_cell = etendue_inputs['aperture']['object']
        img_cell = etendue_inputs['aperture']['image']
        obj_key = obj_cell.keys()
        img_key = img_cell.keys()
        if 'pupil' in obj_key:
            epd = obj_cell['pupil']
            if 'f/#' in img_key:
                fno = img_cell['f/#']
                efl = epd * fno
                imager_inputs = 'f', efl
            if 'NA' in img_key:
                na = img_cell['NA']
                slpk = na2slp(na, n=n_k)
                efl = (epd/2.0)/slpk
                imager_inputs = 'f', efl
        else:
            if 'NA' in obj_key:
                nao = obj_cell['NA']
                slp0 = na2slp(nao, n=n_0)
            elif 'f/#' in obj_key:
                fno = obj_cell['f/#']
                slp0 = -1/(2*fno)

            if 'NA' in img_key:
                na = img_cell['NA']
                slpk = na2slp(na, n=n_k)
            elif 'f/#' in obj_key:
                fno = img_cell['f/#']
                slpk = -1/(2*fno)
            mag = slp0/slpk
            imager_inputs = 'f', efl
    return imager_inputs


def update_aperture_spec(key, cell):
    if 'NA' in key:
        na = cell['NA']
        slp = na2slp(na)
        fno = -1/(2*slp)
        cell['f/#'] = fno
    elif 'f/#' in key:
        fno = cell['f/#']
        slp = -1/(2*fno)
        na = slp2na(slp)
        cell['NA'] = na


def do_etendue_to_imager(inputs, grid, n_0=1.0, n_k=1.0):
    li = dict2d.num_items_by_type(inputs, fld_ape_set, obj_img_set)
    if li['field'] == 2:
        obj_cell = inputs['field']['object']
        img_cell = inputs['field']['image']
        obj_key = obj_cell.keys()
        img_key = img_cell.keys()
        if 'angle' in obj_key:
            obj_ang = obj_cell['angle']
            obj_slp = ang2slp(obj_ang)
            if 'height' in img_key:
                img_ht = img_cell['height']
                efl = img_ht/obj_slp
                imager_inputs = 'f', efl
            else:
                return None
        elif 'height' in obj_key:
            obj_ht = obj_cell['height']
            if 'height' in img_key:
                img_ht = img_cell['height']
                mag = img_ht/obj_ht
                imager_inputs = 'm', mag

    elif li['aperture'] == 2:
        obj_cell = inputs['aperture']['object']
        img_cell = inputs['aperture']['image']
        obj_key = obj_cell.keys()
        img_key = img_cell.keys()
        if 'pupil' in obj_key:
            epd = obj_cell['pupil']
            if 'f/#' in img_key:
                fno = img_cell['f/#']
                efl = epd * fno
                imager_inputs = 'f', efl
            if 'NA' in img_key:
                na = img_cell['NA']
                slpk = na2slp(na, n=n_k)
                efl = (epd/2.0)/slpk
                imager_inputs = 'f', efl
        else:
            if 'NA' in obj_key:
                nao = obj_cell['NA']
                slp0 = na2slp(nao, n=n_0)
            elif 'f/#' in obj_key:
                fno = obj_cell['f/#']
                slp0 = -1/(2*fno)

            if 'NA' in img_key:
                na = img_cell['NA']
                slpk = na2slp(na, n=n_k)
            elif 'f/#' in obj_key:
                fno = img_cell['f/#']
                slpk = -1/(2*fno)
            mag = slp0/slpk
            imager_inputs = 'f', efl
    return imager_inputs


def do_aperture_spec_to_mag(**inputs):
    n_0 = 1.0
    n_k = 1.0
    if 'pupil' in inputs:
        epd = inputs['pupil']
        if 'f/#' in inputs:
            fno = inputs['f/#']
            efl = epd * fno
        if 'na' in inputs:
            na = inputs['na']
            slpk = na2slp(na, n=n_k)
            efl = (epd/2.0)/slpk
    if 'NA' in inputs:
        nao = inputs['NA']
        slp0 = na2slp(nao, n=n_0)
        if 'f/#' in inputs:
            fno = inputs['f/#']
            slpk = -1/(2*fno)
            mag = slp0/slpk
        if 'NA' in inputs:
            na = inputs['NA']
            mag = (nao/na)*(n_k/n_0)
    return mag


def get_aperture_spec_from_mag(imager, **inputs):
    n_0 = 1.0
    n_k = 1.0
    if 'epd' in inputs:
        epd = inputs['epd']
        if 'fno' in inputs:
            fno = imager.f/epd
            aper = ('fno', fno)
        if 'na' in inputs:
            slpk = -(epd/2.0)/imager.f
            na = slp2na(slpk, n=n_k)
            aper = ('na', na)
    if 'nao' in inputs:
        nao = inputs['nao']
        slp0 = na2slp(nao, n=n_0)
        if 'fno' in inputs:
            slpk = imager.m / slp0
            fno = -1/(2.0*slpk)
            aper = ('fno', fno)
        if 'na' in inputs:
            na = imager.m * nao
            aper = ('na', na)
    if 'fno' in inputs:
        fno = inputs['fno']
        if 'epd' in inputs:
            epd = imager.f/fno
            aper = ('epd', epd)
    if 'na' in inputs:
        na = inputs['na']
        slpk = na2slp(na, n=n_k)
        if 'epd' in inputs:
            epd = 2*slpk*imager.f
            aper = ('epd', epd)
    return aper


def get_aperture_spec_from_efl(imager, **inputs):
    n_0 = 1.0
    n_k = 1.0
    if 'epd' in inputs:
        epd = inputs['epd']
        if 'fno' in inputs:
            fno = imager.f/epd
            aper = ('fno', fno)
        if 'na' in inputs:
            slpk = -(epd/2.0)/imager.f
            na = slp2na(slpk, n=n_k)
            aper = ('na', na)
    if 'fno' in inputs:
        fno = inputs['fno']
        if 'epd' in inputs:
            epd = imager.f/fno
            aper = ('epd', epd)
    if 'na' in inputs:
        na = inputs['na']
        slpk = na2slp(na, n=n_k)
        if 'epd' in inputs:
            epd = 2*slpk*imager.f
            aper = ('epd', epd)
    return aper


def do_field_spec_to_efl(**inputs):
    if 'obj_ang' in inputs:
        obj_ang = inputs['obj_ang']
        obj_slp = ang2slp(obj_ang)
        if 'img_ht' in inputs:
            img_ht = inputs['img_ht']
            efl = img_ht/obj_slp
        else:
            return None
    return efl


def do_field_spec_to_mag(**inputs):
    if 'obj_ht' in inputs:
        obj_ht = inputs['obj_ht']
        if 'img_ht' in inputs:
            img_ht = inputs['img_ht']
            mag = img_ht/obj_ht
    elif 'obj_ang' in inputs:
        obj_slp = ang2slp(inputs['obj_ang'])
        if 'img_ang' in inputs:
            img_slp = ang2slp(inputs['img_ang'])
            mag = obj_slp/img_slp
    return mag
