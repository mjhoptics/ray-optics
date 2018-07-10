#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
"""
Created on Tue Jun 26 16:11:58 2018

@author: Michael J. Hayford
"""

from optical.model_constants import ax, pr
from optical.model_constants import ht, slp


def __decode_lens__(lens_package):
    """ decode the lens_package tuple into its constituents

        lens_package is a tuple or a paraxial lens model
        if it's a tuple:
            the first element is a paraxial lens model
            the second element is a tuple with the begining and ending indicies
                into lens model lists
        else its just a lens model
            return it and some range defaults
    """
    if type(lens_package) is tuple:
        lens = lens_package[0]
        bgn, end = lens_package[1]
    else:
        lens = lens_package
        bgn, end = (1, -1)
    return lens, bgn, end


def separation_ratio(lens_package):
    """ the ratio of the backfocus to the mirror separation

        lens_package is a tuple or a paraxial lens model
        if it's a tuple:
            the first element is a paraxial lens model
            the second element is a tuple with the begining and ending indicies
                into lens model lists
    """
    lens, bgn, end = __decode_lens__(lens_package)
    ax_ray = lens[ax]
#    a = ax_ray[ht][end]/ax_ray[slp][end]
#    B = (ax_ray[ht][end] - ax_ray[ht][bgn])/ax_ray[slp][bgn]
#    s = a/B
    m = ax_ray[slp][bgn]/ax_ray[slp][end]
    s = m*ax_ray[ht][end-1]/(ax_ray[ht][bgn] - ax_ray[ht][end-1])
    return s


def mag(lens_package):
    """ the magnification of the input paraxial lens over the specified range

        lens_package is a tuple or a paraxial lens model
        if it's a tuple:
            the first element is a paraxial lens model
            the second element is a tuple with the begining and ending indicies
                into lens model lists
    """
    lens, bgn, end = __decode_lens__(lens_package)
    ax_ray = lens[ax]
    m = ax_ray[slp][bgn]/ax_ray[slp][end]
    return m


def cassegrain(lens_package):
    """ function to calculate the conic constants for an input paraxial lens
    """
    m = mag(lens_package)
    k_pri = -1.0
    k_sec = -4.0*m/(m - 1.0)**2 - 1.0
    return k_pri, k_sec


def dall_kirkham(lens_package):
    """ function to calculate the conic constants for an input paraxial lens
    """
    m = mag(lens_package)
    s = separation_ratio(lens_package)
    k_pri = (s*(m - 1.)*(m + 1.)**2)/((m + s)*m**3) - 1.
    k_sec = 0.0
    return k_pri, k_sec


def ritchey_chretien(lens_package):
    """ function to calculate the conic constants for an input paraxial lens
    """
    m = mag(lens_package)
    s = separation_ratio(lens_package)
    k_pri = -2.0*s/m**3 - 1.0
    k_sec = -(4.0*m*(m - 1.0) + 2.0*(m + s))/(m - 1.0)**3 - 1.0
    return k_pri, k_sec


def spheres(lens_package):
    """ function to revert the conic constants to spherical surfaces
    """
    k_pri = 0.0
    k_sec = 0.0
    return k_pri, k_sec
