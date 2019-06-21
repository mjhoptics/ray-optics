#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
""" module to facilitate first order definition of an optical model

.. Created on Thu May 16 19:57:47 2019

.. codeauthor: Michael J. Hayford
"""

import math
from collections import namedtuple


from rayoptics.optical.model_enums import (PupilType, FieldType)

ideal_imager_keys = ['m', 's', 'sp', 'tt', 'f']
ideal_imager_labels = ["m", "s", "s'", "tt", "f"]
IdealImager = namedtuple('IdealImager', ideal_imager_keys)
""" tuple grouping together first order specifications

    Attributes:
        m: (lateral) magnification
        s: object distance from first principal plane, P1->Obj
        sp: image distance from second principal plane, P2->Img
        tt: total track length, tt = sp - s
        f: focal length
"""


def ideal_imager_setup(**inputs):
    """ Calculate the ideal imaging properties given two independent parameters

    Given 2 system parameters from the following list, this function
    calculates the remaining parameters.

    Note that if specifying ``tt`` and ``f``, their ratio, tt/f, must be
    greater than or equal to 4. A `ValueError` is raised otherwise.

    For a typical system, the value of ``s`` is negative, i.e. the object is to
    the left of the first principal plane.

    Example::

        In [3]: m1s1 = ideal_imager_setup(m=-0.5, s=-10.0); m1s1
        Out[3]: IdealImager(m=-0.5, s=-10.0, sp=5.0, tt=15.0, f=3.333333333333)

        In [4]: s_inf_efl = ideal_imager_setup(s=-math.inf, f=25.0); s_inf_efl
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


def na2slp(na, n=1.0):
    """ convert numerical aperture to slope """
    return n*math.tan(math.asin(na/n))


def slp2na(slp, n=1.0):
    """ convert a ray slope to numerical aperture """
    return n*math.sin(math.atan(slp/n))


def ang2slp(ang):
    """ convert an angle in degrees to a slope """
    return math.tan(math.radians(ang))

    
class SpecSheet:
    """ Top level description of the optical model

    Attributes:
        conjugates: whether the system works at finite or infinite conjugates,
                    either 'Infinite' or 'Finite'

        power: the optical power of the system

        fld_obj:
        fld_img:
        aper_obj:
        aper_img:

    """

    def __init__(self, conjugates='Infinite', cmd_list=None):
        self.set_conjugates(conjugates)

    def set_conjugates(self, conjugates):
        self.conjugates = conjugates
        if conjugates is 'Infinite':
            self.power_str = 'EFL'
            self.fld_obj_str = 'Object angle'
            self.fld_img_str = 'Image height'
            self.aper_obj_str = 'Ent Pupil Diameter'
            self.aper_img_str = 'f/#'
        elif conjugates is 'Finite':
            self.power_str = 'magnification'
            self.fld_obj_str = 'Object height'
            self.fld_img_str = 'Image height'
            self.aper_obj_str = 'Object NA'
            self.aper_img_str = 'NA'

    def do_infinite_conj(self, efl, fld, aper):
        if aper.pupil_type == PupilType.EPD:
            slp0 = aper.value/obj2enp_dist
        if pupil.pupil_type == PupilType.NAO:
            slp0 = n_0*math.tan(math.asin(pupil.value/n_0))
        if pupil.pupil_type == PupilType.FNO:
            slpk = -1./(2.0*pupil.value)
            slp0 = slpk/red
        if pupil.pupil_type == PupilType.NA:
            slpk = n_k*math.tan(math.asin(pupil.value/n_k))
            slp0 = slpk/red
    
            pass


def create_spec_sheet(**kwargs):
    if 'm' in kwargs:
        m = kwargs['m']
    elif 'f' in kwargs:
        f = kwargs['f']
    elif 'obj ht' in kwargs:
        obj_ht = kwargs('obj ht')
    elif 'obj ang' in kwargs:
        obj_ang = kwargs('obj ang')        
    elif 'img ht' in kwargs:
        img_ht = kwargs('img ht')
    elif 'epd' in kwargs:
        epd = kwargs('epd')
    elif 'fno' in kwargs:
        fno = kwargs('fno')
    elif 'nao' in kwargs:
        nao = kwargs('nao')
    elif 'na' in kwargs:
        na = kwargs('na')

    else:
        return None


def do_aperture_spec_to_efl(aper1, aper2):
    n_0 = 1.0
    n_k = 1.0
    if aper1[0] == PupilType.EPD:
        epd = aper1[1]
        if aper2[0] == PupilType.FNO:
            fno = aper2[1]
            efl = epd * fno
        if aper2[0] == PupilType.NA:
            na = aper2[1]
            slpk = na2slp(na, n=n_k)
            efl = (epd/2.0)/slpk
#        if aper1[0] == PupilType.NAO:
#            nao = aper1[1]
#            slp0 = n_0*math.tan(math.asin(pupil.value/n_0))
#            if aper2[0] == PupilType.FNO:
#                fno = aper2[1]
#                slpk =
#                efl = epd * fno
#            if aper2[0] == PupilType.NA:
#                na = aper2[1]
#                slpk = n_k*math.tan(math.asin(na/n_k))
#                efl = (epd/2.0)/slpk
    if aper1[0] == PupilType.FNO:
        fno = aper1[1]
        if aper2[0] == PupilType.EPD:
            epd = aper2[1]
            efl = epd * fno
    if aper1[0] == PupilType.NA:
        na = aper1[1]
        slpk = na2slp(na, n=n_k)
        if aper2[0] == PupilType.EPD:
            epd = aper2[1]
            efl = (epd/2.0)/slpk
    return efl


def do_aperture_spec_to_mag(aper1, aper2):
    if aper1[0].value > aper2[0].value:
        aper1, aper2 = aper2, aper1
    n_0 = 1.0
    n_k = 1.0
#        if aper1[0] == PupilType.EPD:
#            epd = aper1[1]
#            if aper2[0] == PupilType.FNO:
#                fno = aper2[1]
#                efl = epd * fno
#            if aper2[0] == PupilType.NA:
#                na = aper2[1]
#                slpk = na2slp(na, n=n_k)
#                efl = (epd/2.0)/slpk
    if aper1[0] == PupilType.NAO:
        nao = aper1[1]
        slp0 = na2slp(nao, n=n_0)
        if aper2[0] == PupilType.FNO:
            fno = aper2[1]
            slpk = -1/(2*fno)
            mag = slp0/slpk
        if aper2[0] == PupilType.NA:
            na = aper2[1]
            mag = (nao/na)*(n_k/n_0)
    return mag


def get_aperture_spec_from_mag(imager, aper_in, aper_out):
    n_0 = 1.0
    n_k = 1.0
    if aper_in[0] == PupilType.EPD:
        epd = aper_in[1]
        if aper_out == PupilType.FNO:
            fno = imager.f/epd
            aper = (PupilType.FNO, fno)
        if aper_out == PupilType.NA:
            slpk = -(epd/2.0)/imager.f
            na = slp2na(slpk, n=n_k)
            aper = (PupilType.NA, na)
    if aper_in[0] == PupilType.NAO:
        nao = aper_in[1]
        slp0 = na2slp(nao, n=n_0)
        if aper_out == PupilType.FNO:
            slpk = imager.m / slp0
            fno = -1/(2.0*slpk)
            aper = (PupilType.FNO, fno)
        if aper_out == PupilType.NA:
            na = imager.m * nao
            aper = (PupilType.NA, na)
    if aper_in[0] == PupilType.FNO:
        fno = aper_in[1]
        if aper_out == PupilType.EPD:
            epd = imager.f/fno
            aper = (PupilType.EPD, epd)
    if aper_in[0] == PupilType.NA:
        na = aper_in[1]
        slpk = na2slp(na, n=n_k)
        if aper_out == PupilType.EPD:
            epd = 2*slpk*imager.f
            aper = (PupilType.EPD, epd)
    return aper


def get_aperture_spec_from_efl(efl, aper1, aper2):
    n_0 = 1.0
    n_k = 1.0
    if aper1[0] == PupilType.EPD:
        epd = aper1[1]
        if aper2[0] == PupilType.FNO:
            fno = efl/epd
            aper = (PupilType.FNO, fno)
        if aper2[0] == PupilType.NA:
            slpk = -(epd/2.0)/efl
            na = slp2na(slpk, n=n_k)
            aper = (PupilType.NA, na)
#        if aper1[0] == PupilType.NAO:
#            nao = aper1[1]
#            slp0 = n_0*math.tan(math.asin(nao/n_0))
#            if aper2[0] == PupilType.FNO:
#                fno = aper2[1]
#                slpk = 
#                efl = epd * fno 
#            if aper2[0] == PupilType.NA:
#                na = aper2[1]
#                slpk = n_k*math.tan(math.asin(na/n_k))
#                efl = (epd/2.0)/slpk
    if aper1[0] == PupilType.FNO:
        fno = aper1[1]
        if aper2[0] == PupilType.EPD:
            epd = efl/fno
            aper = (PupilType.EPD, epd)
    if aper1[0] == PupilType.NA:
        na = aper1[1]
        slpk = na2slp(na, n=n_k)
        if aper2[0] == PupilType.EPD:
            epd = 2*slpk*efl
            aper = (PupilType.EPD, epd)
    return aper


def do_field_spec_to_efl(fld1, fld2):
    if fld1[0].value > fld2[0].value:
        fld1, fld2 = fld2, fld1
    if fld1[0] == FieldType.OBJ_ANG:
        obj_slp = ang2slp(fld1[1])
        if fld1[0] == FieldType.IMG_HT:
            img_ht = fld2[1]
            efl = img_ht/obj_slp
    return efl


def do_field_spec_to_mag(fld1, fld2):
    if fld1[0].value > fld2[0].value:
        fld1, fld2 = fld2, fld1
    if fld1[0] == FieldType.OBJ_HT:
        obj_ht = fld1[1]
        if fld1[0] == FieldType.IMG_HT:
            img_ht = fld2[1]
            mag = img_ht/obj_ht
    elif fld1[0] == FieldType.OBJ_ANG:
        obj_slp = ang2slp(fld1[1])
        if fld2[0] == FieldType.IMG_ANG:
            img_slp = ang2slp(fld2[1])
            mag = obj_slp/img_slp
    return mag


#dispatch = {
#  (PupilType.EPD, PupilType.FNO): SpecSheet.do_aperture_spec_to_efl,
#  (PupilType.EPD, PupilType.NA): SpecSheet.do_aperture_spec_to_efl,
#  (PupilType.NAO, PupilType.FNO): SpecSheet.do_aperture_spec_to_efl,
#  (PupilType.NAO, PupilType.NA): SpecSheet.do_aperture_spec_to_efl,
#  (FieldType.OBJ_ANG, FieldType.IMG_HT): SpecSheet.do_field_spec_to_efl,
#  (FieldType.OBJ_HT, FieldType.IMG_HT): SpecSheet.do_field_spec_to_efl,
#  }
