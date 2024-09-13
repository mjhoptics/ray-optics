#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
""" module to setup an ideal imager

.. Created on Thu May 16 19:57:47 2019

.. codeauthor: Michael J. Hayford
"""

import math
from collections import namedtuple


ideal_imager_keys = ['m', 's', 'sp', 'tt', 'f']
ideal_imager_labels = ["m", "s", "s'", "tt", "f"]
IdealImager = namedtuple('IdealImager', ideal_imager_keys)
IdealImager.m.__doc__ = "(lateral) magnification"
IdealImager.s.__doc__ = "object distance from first principal plane, P1->Obj"
IdealImager.sp.__doc__ = "image distance from second principal plane, P2->Img"
IdealImager.tt.__doc__ = "total track length, tt = sp - s"
IdealImager.f.__doc__ = "focal length"

""" tuple grouping together first order specifications

    Attributes:
        m: (lateral) magnification
        s: object distance from first principal plane, P1->Obj
        sp: image distance from second principal plane, P2->Img
        tt: total track length, tt = sp - s
        f: focal length
"""


def ideal_imager_setup(**inputs) -> IdealImager:
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
            s = -f*(m - 1)/m  # = tt/(m - 1)
            sp = m*s
        else:
            return IdealImager(m, None, None, None, None)

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
            return IdealImager(None, s, None, None, None)

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
            return IdealImager(None, None, sp, None, None)

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
            return IdealImager(None, None, None, tt, None)

    elif 'f' in inputs:
        f = inputs['f']
        return IdealImager(None, None, None, None, f)

    else:
        return IdealImager(None, None, None, None, None)

    return IdealImager(m, s, sp, tt, f)
