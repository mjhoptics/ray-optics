#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2016 Christian Hill
""" converting a spectrum to a colour

code from the web blog:
    https://scipython.com/blog/converting-a-spectrum-to-a-colour/

@author: Christian Hill
"""

import numpy as np
from scipy.constants import h, c, k
from rayoptics.util.line_intersection import intersect_with_3lines
from pathlib import Path


def xyz_from_xy(x, y):
    """Return the vector (x, y, 1-x-y)."""
    return np.array((x, y, 1-x-y))


class ColourSystem:
    """A class representing a colour system.

    A colour system defined by the CIE x, y and z=1-x-y coordinates of
    its three primary illuminants and its "white point".

    TODO: Implement gamma correction

    """

    # The CIE colour matching function for 380 - 780 nm in 5 nm intervals
    path = Path(__file__).resolve().parent
    cmf = np.loadtxt(path / 'cie-cmf.txt', usecols=(1, 2, 3))

    def __init__(self, red, green, blue, white):
        """Initialise the ColourSystem object.

        Pass vectors (ie NumPy arrays of shape (3,)) for each of the
        red, green, blue  chromaticities and the white illuminant
        defining the colour system.

        """

        # Chromaticities
        self.red, self.green, self.blue = red, green, blue
        self.white = white
        # The chromaticity matrix (rgb -> xyz) and its inverse
        self.M = np.vstack((self.red, self.green, self.blue)).T
        self.MI = np.linalg.inv(self.M)
        # White scaling array
        self.wscale = self.MI.dot(self.white)
        # xyz -> rgb transformation matrix
        self.T = self.MI / self.wscale[:, np.newaxis]

    def xyz_to_rgb(self, xyz, out_fmt=None):
        """Transform from xyz to rgb representation of colour.

        The output rgb components are normalized on their maximum
        value. If xyz is out the rgb gamut, it is desaturated until it
        comes into gamut.

        By default, fractional rgb components are returned; if
        out_fmt='html', the HTML hex string '#rrggbb' is returned.

        """

        rgb = self.T.dot(xyz)
        if np.any(rgb < 0):
            # We're not in the RGB gamut: apply relative colorimetric mapping
            rgb1 = self.relative_colorimetric_gamut_mapping(xyz)
            if np.any(rgb1 < -1.0e-14):
                # We're still not in the RGB gamut: approximate by desaturating
                print("still not in gamut: xyz, rgb0, rgb1", xyz, rgb, rgb1)
                w = - np.min(rgb1)
                rgb1 += w
            rgb = rgb1
        if not np.all(rgb == 0):
            # Normalize the rgb vector
            rgb /= np.max(rgb)

        if out_fmt == 'html':
            return self.rgb_to_hex(rgb)
        return rgb

    def rgb_to_hex(self, rgb):
        """Convert from fractional rgb values to HTML-style hex string."""

        hex_rgb = (255 * rgb).astype(int)
        return '#{:02x}{:02x}{:02x}'.format(*hex_rgb)

    def spec_to_xyz(self, spec):
        """Convert a spectrum to an xyz point.

        The spectrum must be on the same grid of points as the colour-matching
        function, self.cmf: 380-780 nm in 5 nm steps.

        """

        XYZ = np.sum(spec[:, np.newaxis] * self.cmf, axis=0)
        den = np.sum(XYZ)
        if den == 0.:
            return XYZ
        return XYZ / den

    def wvl_to_xyz(self, wv):
        """Convert a wavelength (nm) to an xyz point.

        The wavelength must be on the same grid of points as the
        colour-matching function, self.cmf: 380-780 nm in 5 nm steps.

        """
        index = int((wv-380.0)/5.0)
        if index < 0:
            return np.array(0., 0., 0.)

        try:
            XYZ = self.cmf[index]
        except IndexError:
            XYZ = np.array(0., 0., 0.)
        finally:
            den = np.sum(XYZ)
            if den == 0.:
                return XYZ
            return XYZ / den

    def spec_to_rgb(self, spec, out_fmt=None):
        """Convert a spectrum to an rgb value."""

        xyz = self.spec_to_xyz(spec)
        return self.xyz_to_rgb(xyz, out_fmt)

    def wvl_to_rgb(self, wv, out_fmt=None):
        """Convert a wavelength to an rgb value."""

        xyz = self.wvl_to_xyz(wv)
        return self.xyz_to_rgb(xyz, out_fmt)

    def relative_colorimetric_gamut_mapping(self, xyz):
        """ Apply relative colorimetric gamut mapping to XYZ and return RGB """
        xy_mapped = intersect_with_3lines(xyz[:2], self.white[:2],
                                          (self.blue[:2], self.green[:2]),
                                          (self.green[:2], self.red[:2]),
                                          (self.red[:2], self.blue[:2]))
        rgb = self.T.dot(xyz_from_xy(*xy_mapped))
        return rgb


illuminant_D65 = xyz_from_xy(0.3127, 0.3291)
cs_hdtv = ColourSystem(red=xyz_from_xy(0.67, 0.33),
                       green=xyz_from_xy(0.21, 0.71),
                       blue=xyz_from_xy(0.15, 0.06),
                       white=illuminant_D65)

cs_smpte = ColourSystem(red=xyz_from_xy(0.63, 0.34),
                        green=xyz_from_xy(0.31, 0.595),
                        blue=xyz_from_xy(0.155, 0.070),
                        white=illuminant_D65)

cs_srgb = ColourSystem(red=xyz_from_xy(0.64, 0.33),
                       green=xyz_from_xy(0.30, 0.60),
                       blue=xyz_from_xy(0.15, 0.06),
                       white=illuminant_D65)


def planck(lam, T):
    """ Returns the spectral radiance of a black body at temperature T.

    Returns the spectral radiance, B(lam, T), in W.sr-1.m-2 of a black body
    at temperature T (in K) at a wavelength lam (in nm), using Planck's law.

    """

    lam_m = lam / 1.e9
    fac = h*c/lam_m/k/T
    B = 2*h*c**2/lam_m**5 / (np.exp(fac) - 1)
    return B
