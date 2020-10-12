#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2017 Michael J. Hayford
""" Module for simple optical media definitions

.. Created on Fri Sep 15 17:06:17 2017

.. codeauthor: Michael J. Hayford
"""

from scipy.interpolate import interp1d

from opticalglass.spectral_lines import get_wavelength
import opticalglass.buchdahl as buchdahl


def glass_encode(n, v):
    return str(1000*round((n - 1), 3) + round(v/100, 3))


def glass_decode(gc):
    return round(1.0 + (int(gc)/1000), 3), round(100.0*(gc - int(gc)), 3)


class Medium:
    """ Constant refractive index medium. """

    def __init__(self, nd, lbl, cat=''):
        self.label = lbl
        self.n = nd
        self._catalog_name = cat

    def __repr__(self):
        return 'Medium(' + str(self.n) + ', ' + f"'{self.label}'" + ')'

    def name(self):
        return self.label

    def catalog_name(self):
        return self._catalog_name

    def rindex(self, wv_nm):
        """ returns the interpolated refractive index at wv_nm

        Args:
            wv_nm: the wavelength in nm for the refractive index query

        Returns:
            float: the refractive index at wv_nm
        """
        return self.n


class Air(Medium):
    """ Optical definition for air (low fidelity definition) """

    def __init__(self):
        self.label = 'air'
        self.n = 1.0

    def __repr__(self):
        return 'Air()'


class Glass(Medium):
    """ Optical medium defined by a glass code, i.e. index - V number pair """

    def __init__(self, nd=1.5168, vd=64.17, mat='N-BK7', cat=''):
        self.label = mat
        self._catalog_name = cat
        if mat == 'N-BK7':
            self.n = 1.5168
            self.v = 64.17
        else:
            self.n = nd
            self.v = vd
        self.bdhl_model = buchdahl.Buchdahl2(self.n, self.v)

    def __str__(self):
        return 'Glass ' + self.label + ': ' + glass_encode(self.n, self.v)

    def __repr__(self):
        return ('Glass(nd=' + str(self.n) +
                ', vd=' + str(self.v) +
                ', mat=' + f"'{self.label}'" + ')')

    def sync_to_restore(self):
        if not hasattr(self, 'bdhl_model'):
            self.bdhl_model = buchdahl.Buchdahl2(self.n, self.v)

    def glass_code(self):
        return str(1000*round((self.n - 1), 3) + round(self.v/100, 3))

    def name(self):
        if self.label == '':
            return glass_encode(self.n, self.v)
        else:
            return self.label

    def rindex(self, wv_nm):
        return self.bdhl_model.rindex(wv_nm)

    def update(self, nd, vd):
        self.n = nd
        self.v = vd
        self.bdhl_model.update_model(nd, vd)


class InterpolatedGlass():
    """ Optical medium defined by a list of wavelength/index pairs

    Attributes:
        label: required string identifier for the material
        wvls: list of wavelenghts in nm, used as x axis
        rndx: list of refractive indices corresponding to the values in wvls
        rindex_interp: the interpolation function
    """

    def __init__(self, label, pairs=None, rndx=None, wvls=None):
        self.label = label
        if pairs is not None:
            self.wvls = []
            self.rndx = []
            for w, n in pairs:
                self.wvls.append(w)
                self.rndx.append(n)
        else:
            self.wvls = wvls
            self.rndx = rndx
        self.update()

    def __repr__(self):
        return ('InterpolatedGlass(' + repr(self.label) +
                ', wvls=' + repr(self.wvls) +
                ', rndx=' + repr(self.rndx) + ')')

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['rindex_interp']
        return attrs

    def sync_to_restore(self):
        """ rebuild interpolating function """
        self.update()

    def update(self):
        self.rindex_interp = interp1d(self.wvls, self.rndx, kind='cubic',
                                      assume_sorted=False)

    def glass_code(self):
        nd = self.rindex('d')
        nF = self.rindex('F')
        nC = self.rindex('C')
        vd = (nd - 1)/(nF - nC)
        return str(glass_encode(nd, vd))

    def name(self):
        if self.label == '':
            return self.glass_code()
        else:
            return self.label

    def rindex(self, wv_nm):
        """ returns the interpolated refractive index at wv_nm

        Args:
            wvl: either the wavelength in nm or a string with a spectral line
                 identifier. for the refractive index query

        Returns:
            float: the refractive index at wv_nm

        Raises:
            KeyError: if ``wvl`` is not in the spectra dictionary
        """
        return float(self.rindex_interp(get_wavelength(wv_nm)))
