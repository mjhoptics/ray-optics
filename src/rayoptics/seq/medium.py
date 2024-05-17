#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2017 Michael J. Hayford
""" Module building on :mod:`opticalglass` for ray-optics material support

.. Created on Fri Sep 15 17:06:17 2017

.. codeauthor: Michael J. Hayford
"""
import json
import difflib
import logging

import deprecation
import rayoptics

from scipy.interpolate import interp1d

from rayoptics.util.misc_math import isanumber

from opticalglass import glass as cat_glass
from opticalglass import glassfactory as gfact
from opticalglass.glassfactory import create_glass
from opticalglass import glasserror
from opticalglass import util
from opticalglass.spectral_lines import get_wavelength
from opticalglass.buchdahl import Buchdahl2

from opticalglass import opticalmedium as om
from opticalglass import modelglass as mg
from opticalglass import rindexinfo as rii

logger = logging.getLogger(__name__)


def glass_encode(n: float, v: float) -> str:
    return f'{int(1000*round((n - 1), 3)):3d}.{int(round(10*v, 3)):3d}'


def glass_decode(gc: float) -> float:
    return round(1.0 + (int(gc)/1000), 3), round(100.0*(gc - int(gc)), 3)


def decode_medium(*inputs, **kwargs):
    """ Input utility for parsing various forms of glass input. 

    The **inputs** can have several forms:
        
        - **refractive_index, v-number**: float -> :class:`opticalglass.modelglass.ModelGlass`
        - **refractive_index** only: float -> :class:`opticalglass.opticalmedium.ConstantIndex`
        - **glass_name, catalog_name** as 1 or 2 strings
        - an instance with a `rindex` attribute
        - **air**: str -> :class:`opticalglass.opticalmedium.Air`
        - blank -> defaults to :class:`opticalglass.opticalmedium.Air`

    """
    if 'cat_list' in kwargs:
        cat = kwargs['cat_list']
    else:
        cat = None
    mat = None

    logger.debug(f"num inputs = {len(inputs)}, inputs[0] = {inputs[0]}, "
                  f"{type(inputs[0])}")
    if isanumber(inputs[0]):  # assume all args are numeric
        if len(inputs) <= 1:
            if inputs[0] == 1.0:
                mat = om.Air()
            else:
                mat = om.ConstantIndex(inputs[0], f"n:{inputs[0]:.3f}")
        else:
            if inputs[0] == 1.0:
                mat = om.Air()
            elif inputs[1] == '':
                mat = om.ConstantIndex(inputs[0], f"n:{inputs[0]:.3f}")
            else:
                mat = mg.ModelGlass(inputs[0], inputs[1], '')

    elif isinstance(inputs[0], str):  # string args
        num_str_args = 0
        for tkn in inputs:
            if isinstance(tkn, str) and len(tkn) > 0:
                num_str_args += 1
        if num_str_args == 2:
            name, cat = inputs[0].strip(), inputs[1].strip()
        elif num_str_args == 1:
            if inputs[0].upper() == 'AIR':
                mat = om.Air()
            elif cat is not None:
                name = inputs[0]
            else:
                name, cat = inputs[0].split(',')
                name, cat = name.strip(), cat.strip()
        elif num_str_args == 0:
            mat = om.Air()

        if mat is None:
            try:
                mat = gfact.create_glass(name, cat)
            except glasserror.GlassNotFoundError as gerr:
                logger.info('%s glass data type %s not found',
                             gerr.catalog,
                             gerr.name)
                logger.info('Replacing material with air.')
                mat = om.Air()
    # glass instance args. if they respond to `rindex`, they're in
    elif hasattr(inputs[0], 'rindex'):
        mat = inputs[0]

    if mat:
        logger.info(f"mat = {mat.name()}, {mat.catalog_name()}, {type(mat)}")
    return mat


# --- glass finder base class
class GlassHandlerBase():
    """Base class for glass matching capability.

    This class is used to match catalog glasses to input glass names. It is
    implemented as a class for ease of use by file importers. If the
    glass can be matched up with an existing :mod:`opticalglass` catalog, the
    glass is instantiated and entered into the model. If the glass cannot be
    found, a search for a .smx file of the same name as the model file is made.
    If found, it is a JSON file with a dict that provides an eval() string to
    create an instance to replace the missing glass name. If this file
    isn't found, it is created and contains a JSON template of a dict that has
    the missing glass names as keys; the values are the number of times the
    glass occurs in the file. These values should be replaced with the desired
    eval() string to create a replacement glass.

    Subclasses, e.g. used for different importers, should implement a single
    method that can be called during the import process to return a glass
    instance given an input string.
    """

    def __init__(self, filename):
        self.glass_catalogs = []
        self.glasses_not_found = util.Counter()
        self.track_contents = util.Counter()
        self.filename = None
        if filename:
            self.filename = filename.with_suffix('.smx')
            self.glasses_not_found = self.load_replacements(self.filename)
        self.no_replacements = not self.glasses_not_found

    def load_replacements(self, filename):
        glasses_not_found = util.Counter()
        if filename.exists():
            with filename.open('r') as file:
                glasses_not_found = json.load(file)
        return glasses_not_found

    def save_replacements(self):
        """If unfound glasses, write smx template file. """
        if self.track_contents['glass not found'] > 0:
            if self.glasses_not_found and self.filename:
                fname = self.filename.name.rsplit('.', 1)[0]
                fname += '_tmpl.smx'
                self.filename = self.filename.with_name(fname)
                with self.filename.open('w') as file:
                    json.dump(self.glasses_not_found, file)

    def find_glass(self, name, catalog, always=True):
        """ find `name` glass or a substitute or, if always is True, n=1.5 """

        try:
            if catalog is None or len(catalog) == 0:
                catalog = gfact._cat_names
            medium = create_glass(name, catalog)
        except glasserror.GlassNotFoundError:
            pass
        else:
            self.track_contents['glass found'] += 1
            return medium

        if self.no_replacements:
            medium = self.find_substitute_glass(name)
            if medium is not None:
                self.track_contents['glass substituted'] += 1
                return medium

        medium = self.handle_glass_not_found(name)
        if medium is None and always is True:
            self.track_contents['glass not found'] += 1
            medium = om.ConstantIndex(1.5, 'not '+name)

        return medium

    def find_6_digit_code(self, name):
        """ process `name` as a 6 digit glass code"""
        if isanumber(name):
            if len(name) == 6:
                # process as a 6 digit code, no decimal point
                nd = 1 + float(name[:3])/1000
                vd = float(name[3:])/10
                medium = mg.ModelGlass(nd, vd, name)
                self.track_contents['6 digit code'] += 1
                return medium
        else:
            return None

    def find_substitute_glass(self, name):
        """Try to find a similar glass to ``name``."""

        # create a list of catalogs
        # the original lookup didn't find anything so
        # look in all of our catalogs
        cat_names = [gc.upper() for gc in self.glass_catalogs]
        for gc in gfact._cat_names_uc:
            if gc not in cat_names:
                cat_names.append(gc)

        # create a glass list for the given catalogs
        glist = []
        for glass_cat in cat_names:
            try:
                glass_cat = gfact.get_glass_catalog(glass_cat)
            except glasserror.GlassCatalogNotFoundError:
                pass
            else:
                glist += glass_cat.glass_list

        # Add legacy glasses
        glist += cat_glass.Robb1983Catalog().glass_list

        # decode the input name
        gn_decode = cat_glass.decode_glass_name(name)
        # build an uppercase version for comparisons
        gn_decode_uc = gn_decode[0][0].upper() + gn_decode[0][1]

        subs_glasses = []
        for g in glist:
            gn_decode, gn, gc = g
            if gn_decode_uc == gn_decode[0][0].upper() + gn_decode[0][1]:
                subs_glasses.append(g)

        if len(subs_glasses):
            possibilities = [gn for gn_decode, gn, gc in subs_glasses]
            matches = difflib.get_close_matches(name, possibilities)
            if len(matches) > 0:
                gn = matches[0]
                gc = next((g[2] for g in subs_glasses if g[1] == gn), None)
            else:
                gn_decode, gn, gc = subs_glasses[0]
            medium = create_glass(gn, gc)
            self.glasses_not_found[name] = gn, gc
            return medium
        else:
            return None

    def handle_glass_not_found(self, name):
        """Record missing glasses or create new replacement glass instances."""

        if self.no_replacements:                # track the number of times
            self.glasses_not_found[name] += 1   # each missing glass is used
            return None

        else:  # create a new instance of the replacement glass
            if name in self.glasses_not_found:
                val = self.glasses_not_found[name]
                if len(val) == 2:  # call create_glass with name and catalog
                    gn, gc = val
                    return create_glass(gn, gc)
                else:  # eval code to create a new glass instance
                    mat = eval(self.glasses_not_found[name])
                    if hasattr(mat, 'convert_to_OG'):
                        mat = mat.convert_to_OG()
                    return mat
            else:
                return None


# --- material definitions
class Medium:
    """ Constant refractive index medium. """

    @deprecation.deprecated(deprecated_in="0.8.4", removed_in="0.9.0",
                            current_version=rayoptics.__version__,
                            details="Use opticalglass.opticalmedium."
                                    "ConstantIndex instead")
    def __init__(self, nd, lbl, cat=''):
        self.label = lbl
        self.n = nd
        self._catalog_name = cat

    def __repr__(self):
        return ('Medium(' + str(self.n) + ', ' + f"'{self.label}'" +
                ', cat=' + f"'{self._catalog_name}'" + ')')

    def convert_to_OG(self):
        """ returns an equivalent material from the :mod:`opticalglass` package. """
        if not hasattr(self, '_catalog_name'):
            self._catalog_name = ''
        return om.ConstantIndex(self.n, self.label, cat=self._catalog_name)

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

    @deprecation.deprecated(deprecated_in="0.8.4", removed_in="0.9.0",
                            current_version=rayoptics.__version__,
                            details="Use opticalglass.opticalmedium."
                                    "Air instead")
    def __init__(self):
        self.label = 'air'
        self.n = 1.0

    def __repr__(self):
        return 'Air()'

    def convert_to_OG(self):
        """ returns an equivalent material from the :mod:`opticalglass` package. """
        return om.Air()

    def name(self):
        return self.label

    def catalog_name(self):
        return ''


class Glass(Medium):
    """ Optical medium defined by a glass code, i.e. index - V number pair """

    @deprecation.deprecated(deprecated_in="0.8.4", removed_in="0.9.0",
                            current_version=rayoptics.__version__,
                            details="Use opticalglass.modelglass."
                                    "ModelGlass instead")
    def __init__(self, nd=1.5168, vd=64.17, mat='', cat=''):
        self.label = mat
        self._catalog_name = cat
        self.n = nd
        self.v = vd
        self.bdhl_model = Buchdahl2(self.n, self.v)

    def __str__(self):
        return 'Glass ' + self.label + ': ' + glass_encode(self.n, self.v)

    def __repr__(self):
        return ('Glass(nd=' + str(self.n) +
                ', vd=' + str(self.v) +
                ', mat=' + f"'{self.label}'" +
                ', cat=' + f"'{self._catalog_name}'" + ')')

    def convert_to_OG(self):
        """ returns an equivalent material from the :mod:`opticalglass` package. """
        if not hasattr(self, '_catalog_name'):
            self._catalog_name = ''
        return mg.ModelGlass(self.n, self.v, 
                             self.label, cat=self._catalog_name)

    def sync_to_restore(self):
        if not hasattr(self, 'bdhl_model'):
            self.bdhl_model = Buchdahl2(self.n, self.v)

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

    @deprecation.deprecated(deprecated_in="0.8.4", removed_in="0.9.0",
                            current_version=rayoptics.__version__,
                            details="Use opticalglass.opticalmedium."
                                    "InterpolatedMedium instead")
    def __init__(self, label, pairs=None, rndx=None, wvls=None, cat=''):
        self.label = label
        self._catalog = cat
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
        return ('InterpolatedGlass(' + f"'{self.label}'" +
                ', cat=' + f"'{self._catalog}'" +
                ', wvls=' + repr(self.wvls) +
                ', rndx=' + repr(self.rndx) + ')')

    def convert_to_OG(self):
        """ returns an equivalent material from the :mod:`opticalglass` package. """
        if not hasattr(self, '_catalog_name'):
            self._catalog_name = ''
        return om.InterpolatedMedium(self.label, 
                                     rndx=self.rndx, wvls=self.wvls, 
                                     cat=self._catalog)

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

    def catalog_name(self):
        """ returns the glass catalog name """
        return self._catalog

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
