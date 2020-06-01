#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
""" Support for spectral line data

.. codeauthor: Michael J. Hayford
"""


spectral_lines = [[2325.42, '', 'infrared mercury line', 'Hg'],
                  [1970.09, '', 'infrared mercury line', 'Hg'],
                  [1529.582, '', 'infrared mercury line', 'Hg'],
                  [1060.0, '', 'neodymium glass laser', 'Nd'],
                  [1013.98, 't', 'infrared mercury line', 'Hg'],
                  [852.11, 's', 'infrared cesium line', 'Cs'],
                  [706.5188, 'r', 'red helium line', 'He'],
                  [656.2725, 'C', 'red hydrogen line', 'H'],
                  [643.8469, "C'", 'red cadmium line', 'Cd'],
                  [632.8, '', 'helium-neon-gas-laser', 'He-Ne'],
                  [589.2938, 'D', 'center of double sodium line', 'Na'],
                  [587.5618, 'd', 'yellow helium line', 'He'],
                  [546.074, 'e', 'green mercury line', 'Hg'],
                  [486.1327, 'F', 'blue hydrogen line', 'H'],
                  [479.9914, "F'", 'blue cadmium line', 'Cd'],
                  [435.8343, 'g', 'blue mercury line', 'Hg'],
                  [404.6561, 'h', 'violet mercury line', 'Hg'],
                  [365.014, 'i', 'ultraviolet mercury line', 'Hg'],
                  [334.1478, '', 'ultraviolet mercury line', 'Hg'],
                  [312.5663, '', 'ultraviolet mercury line', 'Hg'],
                  [296.7278, '', 'ultraviolet mercury line', 'Hg'],
                  [280.4, '', 'ultraviolet mercury line', 'Hg'],
                  [248.3, '', 'ultraviolet mercury line', 'Hg']]


spectra = {'Nd': 1060.0,
           't': 1013.98,
           's': 852.11,
           'r': 706.5188,
           'C': 656.2725,
           "C'": 643.8469,
           'He-Ne': 632.8,
           'D': 589.2938,
           'd': 587.5618,
           'e': 546.074,
           'F': 486.1327,
           "F'": 479.9914,
           'g': 435.8343,
           'h': 404.6561,
           'i': 365.014}


spectra_uc = {key.upper(): val for key, val in spectra.items()}


def get_wavelength(wvl):
    """Return wvl in nm, where wvl can be a spectral line

    Example::

        In [1]: from rayoptics.util.spectral_lines import *

        In [2]: wl_e = get_wavelength('e'); wl_e
        Out[2]: 546.074

        In [3]: wl_HeNe = get_wavelength('He-Ne'); wl_HeNe
        Out[3]: 632.8

        In [4]: wl_550 = get_wavelength(550); wl_550
        Out[4]: 550.0

        In [5]: wl_fl = get_wavelength(555.0); wl_fl
        Out[5]: 555.0

        In [6]: wl_f = get_wavelength('f'); wl_f
        Out[6]: 486.1327

    Args:
        wvl: either the wavelength in nm or a string with a spectral line
             identifier. Case insensitive.

    Returns:
        float: the wavelength in nm

    Raises:
        KeyError: if ``wvl`` is not in the spectra dictionary
    """
    if isinstance(wvl, float):
        return wvl
    elif isinstance(wvl, int):
        return float(wvl)
    else:
        return spectra_uc[wvl.upper()]
