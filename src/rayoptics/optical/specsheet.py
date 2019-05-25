#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
""" module to facilitate first order definition of an optical model

.. Created on Thu May 16 19:57:47 2019

.. codeauthor: Michael J. Hayford
"""

from collections import namedtuple
import math

from rayoptics.optical.model_enums import (PupilType, FieldType)

FiniteConj = namedtuple('FiniteConj', ['m', 's', 'sp', 'tt', 'f'])
""" tuple grouping together first order finite conjugate specifications

    Attributes:
        m: (lateral) magnification
        s: object distance from first principal plane, P1->Obj
        sp: image distance from second principal plane, P2->Img
        tt: total track length, tt = sp - s
        f: focal length
"""

InfiniteConj = namedtuple('InfiniteConj', ['m', 's', 'sp', 'tt', 'f'])
""" tuple grouping together first order infinite conjugate specifications

    Attributes:
        sp: image distance from second principal plane, P2->Img
        f: focal length
"""


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

    Lightweight manager class to manage connections between a model and
    a ui view.

    The main function of AppManager is refresh_gui(). This is called after
    a user input to the gui to update the model and call a refresh function
    for each ui view of that model.

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

    def get_aperture_spec_from_mag(fconj, aper1, aper2):
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
        if aper1[0] == PupilType.NAO:
            nao = aper1[1]
            slp0 = na2slp(nao, n=n_0)
            if aper2[0] == PupilType.FNO:
                efl = epd * fno
                slpk = mag / slp0
                fno = -1/(2.0*slpk)
            if aper2[0] == PupilType.NA:
                na = fconj.m * nao
                slpk = na2slp(na, n=n_k)
                efl = (epd/2.0)/slpk
        if aper1[0] == PupilType.FNO:
            fno = aper1[1]
            if aper2[0] == PupilType.EPD:
                epd = efl/fno
                aper = (PupilType.EPD, epd)
        if aper1[0] == PupilType.NA:
            na = aper1[1]
            slpk = n_k*math.tan(math.asin(na/n_k))
            if aper2[0] == PupilType.EPD:
                epd = 2*slpk*efl
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


def do_finite_setup(arg1, arg2):
    """ Calculate the first order system parameters for finite conjugates

    Given 2 system parameters from the following list, this function
    calculates the remaining parameters.

    Note that if specifying ``tt`` and ``f``, their ratio, tt/f, must be
    greater than or equal to 4. A `ValueError` is raised otherwise.

    For a typical system, the value of ``s`` is negative, i.e. the object is to
    the left of the first principal plane.

    First order parameters:
        - **m**: (lateral) magnification
        - **s**: object distance from first principal plane, P1->Obj
        - **sp**: image distance from second principal plane, P2->Img
        - **tt**: total track length, tt = sp - s
        - **f**: focal length

    Example::

        In [1]: m1 = ('m', -0.5)
        In [2]: s1 = ('s', -10.)
        In [3]: m1s1 = do_finite_setup(m1, s1); m1s1
        Out[3]: FiniteConj(m=-0.5, s=-10.0, sp=5.0, tt=15.0, f=3.333333333333)

    Args:
        arg1: first tuple of parameter key and value
        arg2: secondtuple of parameter key and value

    Returns:
        :class:`FiniteConj` namedtuple

    Raises:
        ValueError: if tt/f < 4

    """
    finite_parms = ['m', 's', 'sp', 'tt', 'f']
    if finite_parms.index(arg1[0]) > finite_parms.index(arg2[0]):
        arg1, arg2 = arg2, arg1

    if arg1[0] == 'm':
        m = arg1[1]
        if arg2[0] == 's':
            s = arg2[1]
            sp = m*s
            tt = sp - s
            f = s*sp/(s - sp)
        elif arg2[0] == 'sp':
            sp = arg2[1]
            s = sp/m
            tt = sp - s
            f = s*sp/(s - sp)
        elif arg2[0] == 'tt':
            tt = arg2[1]
            s = tt/(m - 1)
            sp = m*s
            f = s*sp/(s - sp)
        else:  # arg2[0] == 'f':
            f = arg2[1]
            tt = -f*(m - 1)**2/m
            s = tt/(m - 1)
            sp = m*s

    elif arg1[0] == 's':
        # arrange calculations so that s=-inf is handled gracefully
        s = arg1[1]
        if arg2[0] == 'sp':
            sp = arg2[1]
            f = 1/(1/sp - 1/s)
            m = sp/s
            tt = sp - s
        elif arg2[0] == 'tt':
            tt = arg2[1]
            m = 1 + tt/s
            sp = m*s
            f = s*sp/(s - sp)
        elif arg2[0] == 'f':
            f = arg2[1]
            m = f/(s + f)
            sp = 1/(1/f + 1/s)
            tt = sp - s

    elif arg1[0] == 'sp':
        # arrange calculations so that sp=inf is handled gracefully
        sp = arg1[1]
        if arg2[0] == 'tt':
            tt = arg2[1]
            m = sp/(sp - tt)
            s = sp/m
            f = s*sp/(s - sp)
        elif arg2[0] == 'f':
            f = arg2[1]
            m = (f - sp)/f
            s = 1/(1/sp - 1/f)
            tt = sp - s
    elif arg1[0] == 'tt':
        tt = arg1[1]
        if arg2[0] == 'f':
            f = arg2[1]
            ttf = tt/f
            # tt/f >= 4, else no solution
            # pick root (+) that gives |s|>=|sp|, i.e. -1 <= m < 0
            m = ((2 - ttf) + math.sqrt(ttf*(ttf - 4)))/2
            s = tt/(m - 1)
            sp = m*s

    return FiniteConj(m, s, sp, tt, f)

dispatch = {
  (PupilType.EPD, PupilType.FNO): SpecSheet.do_aperture_spec_to_efl,
  (PupilType.EPD, PupilType.NA): SpecSheet.do_aperture_spec_to_efl,
  (PupilType.NAO, PupilType.FNO): SpecSheet.do_aperture_spec_to_efl,
  (PupilType.NAO, PupilType.NA): SpecSheet.do_aperture_spec_to_efl,
  (FieldType.OBJ_ANG, FieldType.IMG_HT): SpecSheet.do_field_spec_to_efl,
  (FieldType.OBJ_HT, FieldType.IMG_HT): SpecSheet.do_field_spec_to_efl,
  }
