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

from rayoptics.util.misc_math import isanumber

import opticalglass as og
from opticalglass import glass as cat_glass
from opticalglass import glassfactory as gfact
from opticalglass import opticalmedium as om
from opticalglass import modelglass as mg
from opticalglass import rindexinfo as rii
from opticalglass import glasserror

logger = logging.getLogger(__name__)


def glass_encode(n: float, v: float) -> str:
    return f'{int(1000*round((n - 1), 3)):3d}.{int(round(10*v, 3)):3d}'


def glass_decode(gc: float) -> tuple[float, float]:
    return round(1.0 + (int(gc)/1000), 3), round(100.0*(gc - int(gc)), 3)


def decode_medium(*inputs, **kwargs) -> om.OpticalMedium:
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
        if len(inputs) == 1:
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
    glass is instantiated and entered into the model. Searching includes the 
    custom_glass_registry from the `opticalglass` package. If the glass cannot
    be found, a search for a .smx file of the same name as the model file is 
    made. If found, it is a JSON file with a dict that provides an eval() 
    string to create an instance to replace the missing glass name. If this 
    file isn't found, it is created and contains a JSON template of a dict that 
    has the missing glass names as keys; the values are the number of times the 
    glass occurs in the file. These values should be replaced with the desired
    eval() string to create a replacement glass.

    Subclasses, e.g. used for different importers, should implement a single
    method that can be called during the import process to return a glass
    instance given an input string.
    """

    def __init__(self, filename):
        self.glass_catalogs = []
        self.glasses_not_found = og.util.Counter()
        self.track_contents = og.util.Counter()
        self.filename = None
        if filename:
            self.filename = filename.with_suffix('.smx')
            self.glasses_not_found = self.load_replacements(self.filename)
        self.no_replacements = not self.glasses_not_found
        self.custom_glass_registry = {key[0]: (key[1], value) 
                                      for key, value in gfact._custom_glass_registry.items()}

    def load_replacements(self, filename):
        glasses_not_found = og.util.Counter()
        if filename.exists():
            with filename.open('r') as file:
                glasses_not_found = json.load(file)
            for name, val in glasses_not_found.items():
                if len(val) == 2:  # call create_glass with name and catalog
                    gn, gc = val
                    mat = gfact.create_glass(gn, gc)
                else:  # eval code to create a new glass instance
                    mat = eval(val)
                gfact.register_glass(mat)
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

    def find_glass(self, name, catalog, always=True) -> om.OpticalMedium|None:
        """ find `name` glass or a substitute or, if always is True, n=1.5 
        
        Include searching the custom_glass_registry from the 
        opticalglass package.
        """

        try:
            if catalog is None or len(catalog) == 0:
                if name in self.custom_glass_registry:
                    catalog = self.custom_glass_registry[name][0]
                else:
                    catalog = gfact._cat_names
            medium = gfact.create_glass(name, catalog)
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

    def find_6_digit_code(self, name) -> om.OpticalMedium|None:
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

    def find_substitute_glass(self, name) -> om.OpticalMedium|None:
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
            medium = gfact.create_glass(gn, gc)
            self.glasses_not_found[name] = gn, gc
            return medium
        else:
            return None

    def handle_glass_not_found(self, name) -> om.OpticalMedium|None:
        """Record missing glasses or create new replacement glass instances."""

        if self.no_replacements:                # track the number of times
            self.glasses_not_found[name] += 1   # each missing glass is used
            return None

        else:  # create a new instance of the replacement glass
            if name in self.glasses_not_found:
                val = self.glasses_not_found[name]
                if len(val) == 2:  # call create_glass with name and catalog
                    gn, gc = val
                    mat = gfact.create_glass(gn, gc)
                else:  # eval code to create a new glass instance
                    mat = eval(self.glasses_not_found[name])
                gfact.register_glass(mat)
                return mat
            else:
                return None
