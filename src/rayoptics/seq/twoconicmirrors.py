#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" calculate conic constants for different 2 mirror configurations

.. Created on Tue Jun 26 16:11:58 2018
.. |ParaxialModel| replace:: :class:`~.paraxialdesign.ParaxialModel`

.. codeauthor: Michael J. Hayford
"""

from rayoptics.optical.model_constants import ht, slp


def __decode_lens__(lens_package):
    """ decode the lens_package tuple into its constituents

        Args:
            lens_package: a tuple or a |ParaxialModel|. If it's a tuple:

               - the first element is a |ParaxialModel|
               - the second element is a tuple with the begining and ending
                 indicies into lens model lists

        Returns:
            |ParaxialModel|, beginning and ending indicies into paraxial model
    """
    if type(lens_package) is tuple:
        parax_model = lens_package[0]
        bgn, end = lens_package[1]
    else:
        parax_model = lens_package
        bgn, end = (1, -1)
    return parax_model, bgn, end


def _separation_ratio(lens_package):
    """ the ratio of the backfocus to the mirror separation

        Args:
            see description in __decode_lens__
    """
    parax_model, bgn, end = __decode_lens__(lens_package)
    ax_ray = parax_model.ax
#    a = ax_ray[ht][end]/ax_ray[slp][end]
#    B = (ax_ray[ht][end] - ax_ray[ht][bgn])/ax_ray[slp][bgn]
#    s = a/B
    m = ax_ray[bgn][slp]/ax_ray[end][slp]
    s = m*ax_ray[end-1][ht]/(ax_ray[bgn][ht] - ax_ray[end-1][ht])
    return s


def _mag(lens_package):
    """ the magnification of the input paraxial lens over the specified range

        Args:
            see description in __decode_lens__
    """
    parax_model, bgn, end = __decode_lens__(lens_package)
    ax_ray = parax_model.ax
    m = ax_ray[bgn][slp]/ax_ray[end][slp]
    return m


def cassegrain(lens_package):
    """ calculate the conic constants for a cassegrain telescope

        Args:
            lens_package: a tuple or a |ParaxialModel|. If it's a tuple:

               - the first element is a |ParaxialModel|
               - the second element is a tuple with the begining and ending
                 indicies into lens model lists

        Returns:
            the conic constants of the primary and secondary mirrors
    """
    m = _mag(lens_package)
    k_pri = -1.0
    k_sec = -4.0*m/(m - 1.0)**2 - 1.0
    return k_pri, k_sec


def dall_kirkham(lens_package):
    """ calculate the conic constants for a dall-kirkham telescope

        Args:
            lens_package: a tuple or a |ParaxialModel|. If it's a tuple:

               - the first element is a |ParaxialModel|
               - the second element is a tuple with the begining and ending
                 indicies into lens model lists

        Returns:
            the conic constants of the primary and secondary mirrors
    """
    m = _mag(lens_package)
    s = _separation_ratio(lens_package)
    k_pri = (s*(m - 1.)*(m + 1.)**2)/((m + s)*m**3) - 1.
    k_sec = 0.0
    return k_pri, k_sec


def ritchey_chretien(lens_package):
    """ calculate the conic constants for a ritchey-chretien telescope

        Args:
            lens_package: a tuple or a |ParaxialModel|. If it's a tuple:

               - the first element is a |ParaxialModel|
               - the second element is a tuple with the begining and ending
                 indicies into lens model lists

        Returns:
            the conic constants of the primary and secondary mirrors
    """
    m = _mag(lens_package)
    s = _separation_ratio(lens_package)
    k_pri = -2.0*s/m**3 - 1.0
    k_sec = -(4.0*m*(m - 1.0) + 2.0*(m + s))/(m - 1.0)**3 - 1.0
    return k_pri, k_sec


def spheres(lens_package):
    """ function to revert the conic constants to spherical surfaces

        Args:
            lens_package: a tuple or a |ParaxialModel|. If it's a tuple:

               - the first element is a |ParaxialModel|
               - the second element is a tuple with the begining and ending
                 indicies into lens model lists

        Returns:
            the conic constants of the primary and secondary mirrors
    """
    k_pri = 0.0
    k_sec = 0.0
    return k_pri, k_sec
