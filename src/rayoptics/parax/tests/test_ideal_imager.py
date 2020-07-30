#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" tests for ideal imager

.. Created on Wed Jun 12 17:57:06 2019

.. codeauthor: Michael J. Hayford
"""

import math
import pytest
from rayoptics.parax.idealimager import ideal_imager_setup


def test_ideal_imager_setup(**kwargs):
    m1 = -0.5
    s1 = -10.
    m1s1 = ideal_imager_setup(m=m1, s=s1)
    sp1 = m1s1.sp
    tt1 = m1s1.tt
    f1 = m1s1.f
    assert ideal_imager_setup(m=m1, s=s1) == ideal_imager_setup(s=s1, m=m1)
    assert ideal_imager_setup(m=m1, sp=sp1) == ideal_imager_setup(sp=sp1, m=m1)
    assert ideal_imager_setup(m=m1, tt=tt1) == ideal_imager_setup(tt=tt1, m=m1)
    assert ideal_imager_setup(m=m1, f=f1) == ideal_imager_setup(f=f1, m=m1)

    tt2 = 4.0
    f2 = 1.0
    tt2f2 = ideal_imager_setup(tt=tt2, f=f2)
    assert tt2f2.m == -1
    assert tt2f2.s == -2
    assert tt2f2.sp == 2
    f3 = 2.0
    with pytest.raises(ValueError):
        ideal_imager_setup(tt=tt2, f=f3)

    s_inf = -math.inf
    efl = 25.0
    s_inf_efl = ideal_imager_setup(s=s_inf, f=efl)
    assert s_inf_efl.m == -0
    assert s_inf_efl.s == -math.inf
    assert s_inf_efl.sp == efl
