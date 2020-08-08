#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2020 Michael J. Hayford
"""

.. Created on Fri Jul 31 15:40:21 2020

.. codeauthor: Michael J. Hayford
"""
import logging
import math
import json

from rayoptics.optical.opticalmodel import OpticalModel
from rayoptics.optical.model_enums import DimensionType as dt
from rayoptics.optical.model_enums import DecenterType as dec
from rayoptics.elem.surface import (DecenterData, Circular, Rectangular,
                                    Elliptical)
from rayoptics.elem import profiles
from rayoptics.seq.medium import Medium, Air, Glass, InterpolatedGlass
from rayoptics.raytr.opticalspec import Field
from rayoptics.util.misc_math import isanumber

from opticalglass import glassfactory as gfact
from opticalglass import glasserror

_glass_handler = None


def read_lens(filename, **kwargs):
    ''' given a Zemax .zmx filename, return an OpticalModel  '''
    global _glass_handler
    logging.basicConfig(filename='zmx_read_lens.log',
                        filemode='w',
                        level=logging.DEBUG)
    opt_model = OpticalModel(do_init=False)

    with filename.open() as file:
        inpt = file.read()
    input_lines = inpt.splitlines()

    _glass_handler = GlassHandler(filename)

    for i, line in enumerate(input_lines):
        process_line(opt_model, line, i+1)

    post_process_input(opt_model)
    _glass_handler.save_replacements()

    opt_model.update_model()

    return opt_model


def process_line(opt_model, line, line_no):
    global _glass_handler
    sm = opt_model.seq_model
    osp = opt_model.optical_spec
    cur = sm.cur_surface
    if not line.strip():
        return
    line = line.strip().split(" ", 1)
    cmd = line[0]
    inputs = len(line) == 2 and line[1] or ""
    if cmd == "UNIT":
        dim = inputs.split()[0]
        if dim == 'MM':
            dim = dt.MM
        elif dim == 'IN' or dim == 'INCH':
            dim = dt.IN
        opt_model.system_spec.dimensions = dim
    elif cmd == "NAME":
        opt_model.system_spec.title = inputs.strip("\"")
    elif cmd == "NOTE":
        opt_model.note = inputs.strip("\"")
    elif cmd == "SURF":
        s, g = sm.insert_surface_and_gap()
    elif cmd == "CURV":
        s = sm.ifcs[cur]
        s.profile.cv = float(inputs.split()[0])
    elif cmd == "DISZ":
        g = sm.gaps[cur]
        g.thi = float(inputs)

    elif _glass_handler(sm, cur, cmd, inputs):
        pass

    elif cmd == "DIAM":
        s = sm.ifcs[cur]
        s.set_max_aperture(float(inputs.split()[0]))
    elif cmd == "STOP":
        sm.set_stop()

    elif cmd == "WAVM":  # WAVM 1 0.55000000000000004 1
        sr = osp.spectral_region
        inputs = inputs.split()
        new_wvl = float(inputs[1])*1e+3
        if new_wvl not in sr.wavelengths:
            sr.wavelengths.append(new_wvl)
            sr.spectral_wts.append(float(inputs[2]))  # needs check
    elif cmd == "WAVL":
        s.wavelengths = [float(i)*1e-6 for i in inputs.split() if i]

    elif pupil_data(opt_model, cmd, inputs):
        pass

    elif field_spec_data(opt_model, cmd, inputs):
        pass

    elif handle_profiles(opt_model, cur, cmd, inputs):
        pass

    elif cmd in ("OPDX",  # opd
                 "RAIM",  # ray aiming
                 "CONF",  # configurations
                 "ENPD", "PUPD",  # pupil
                 "EFFL",  # focal lengths
                 "VERS",  # version
                 "MODE",  # mode
                 "NOTE",  # note
                 "HIDE",  # surface hide
                 "MIRR",  # surface is mirror
                 "PARM",  # aspheric parameters
                 "SQAP",  # square aperture?
                 "XDAT", "YDAT",  # xy toroidal data
                 "OBNA",  # object na
                 "PKUP",  # pickup
                 "MAZH", "CLAP", "PPAR", "VPAR", "EDGE", "VCON",
                 "UDAD", "USAP", "TOLE", "PFIL", "TCED", "FNUM",
                 "TOL", "MNUM", "MOFF", "FTYP", "SDMA", "GFAC",
                 "PUSH", "PICB", "ROPD", "PWAV", "POLS", "GLRS",
                 "BLNK", "COFN", "NSCD", "GSTD", "DMFS", "ISNA",
                 "VDSZ", "ENVD", "ZVDX", "ZVDY", "ZVCX", "ZVCY",
                 "ZVAN", "XFLN", "YFLN", "VDXN", "VDYN", "VCXN",
                 "VCYN", "VANN", "FWGT", "FWGN", "WWGT", "WWGN",
                 "WAVN", "XFLD", "YFLD", "MNCA", "MNEA",
                 "MNCG", "MNEG", "MXCA", "MXCG", "RGLA", "TRAC",
                 "FLAP", "TCMM", "FLOA", "PMAG", "TOTR", "SLAB",
                 "POPS", "COMM", "PZUP", "LANG", "FIMP", "COAT",
                 ):
        logging.info('Line %d: Command %s not supported', line_no, cmd)
    else:
        print(cmd, "not handled", inputs)


def post_process_input(opt_model):
    sm = opt_model.seq_model
    sm.gaps.pop()
    conj_type = 'finite'
    if math.isinf(sm.gaps[0].thi):
        sm.gaps[0].thi = 1e10
        conj_type = 'infinite'

    osp = opt_model.optical_spec
    sr = osp.spectral_region
    if len(sr.wavelengths) > 1:
        sr.wavelengths.pop()
        sr.spectral_wts.pop()
    sr.reference_wvl = len(sr.wavelengths)//2

    fov = osp.field_of_view
    if conj_type == 'finite':
        fov.key = 'field', 'object', 'height'
    elif conj_type == 'infinite':
        fov.key = 'field', 'object', 'angle'

    max_fld, max_fld_idx = fov.max_field()
    fov.fields = [f for f in fov.fields[:max_fld_idx+1]]
    # switch vignetting definition to asymmetric vly, vuy style
    # need to verify this is how this works
    for f in fov.fields:
        f.vlx = f.vcx + f.vdx
        f.vux = f.vcx - f.vdx
        f.vly = f.vcy + f.vdy
        f.vuy = f.vcy - f.vdy


def log_cmd(label, cmd, inputs):
    logging.debug("%s: %s %s", label, cmd, str(inputs))


def handle_profiles(optm, cur, cmd, inputs):
    sm = optm.seq_model
    if cmd == "TYPE":
        typ = inputs.split()[0]
        if typ == 'EVENASPH':
            cur_profile = sm.ifcs[cur].profile
            new_profile = profiles.mutate_profile(cur_profile,
                                                  'EvenPolynomial')
            sm.ifcs[cur].profile = new_profile
    elif cmd == "CONI":
        s = sm.ifcs[cur]
        cur_profile = s.profile
        if not hasattr(cur_profile, 'cc'):
            s.profile = profiles.mutate_profile(cur_profile, 'Conic')
        s.profile.cc = float(inputs.split()[0])
    elif cmd == "PARM":
        s = sm.ifcs[cur]
        i, coef = inputs.split()
        coef = float(coef)
        s.profile.coefs.append(coef)
    else:
        return False
    return True


def pupil_data(optm, cmd, inputs):
    # FNUM 2.1 0
    # OBNA 1.5E-1 0
    pupil = optm.optical_spec.pupil
    if cmd == 'FNUM':
        pupil.key = 'aperture', 'image', 'f/#'
    elif cmd == 'OBNA':
        pupil.key = 'aperture', 'object', 'NA'
    else:
        return False

    pupil.value = float(inputs.split()[0])

    log_cmd("pupil_data", cmd, inputs)

    return True


def field_spec_data(optm, cmd, inputs):
    # XFLN 0 0 0 0 0 0 0 0 0 0 0 0
    # YFLN 0 8.0 1.36E+1 0 0 0 0 0 0 0 0 0
    # FWGN 1 1 1 1 1 1 1 1 1 1 1 1
    # VDXN 0 0 0 0 0 0 0 0 0 0 0 0
    # VDYN 0 0 0 0 0 0 0 0 0 0 0 0
    # VCXN 0 0 0 0 0 0 0 0 0 0 0 0
    # VCYN 0 0 0 0 0 0 0 0 0 0 0 0
    # VANN 0 0 0 0 0 0 0 0 0 0 0 0
    fov = optm.optical_spec.field_of_view
    if cmd == 'XFLN' or cmd == 'YFLN':
        fov.key = 'field', 'object', 'angle'
        attr = cmd[0].lower()
    elif cmd == 'VDXN' or cmd == 'VDYN':
        attr = 'vd' + cmd[2].lower()
    elif cmd == 'VCXN' or cmd == 'VCYN':
        attr = 'vc' + cmd[2].lower()
    elif cmd == 'VANN':
        attr = 'van'
    elif cmd == 'FWGN':
        attr = 'wt'
    else:
        return False

    inputs = inputs.split()

    if len(fov.fields) != len(inputs):
        fov.fields = [Field() for f in range(len(inputs))]

    for i, f in enumerate(fov.fields):
        f.__setattr__(attr, float(inputs[i]))

    log_cmd("field_spec_data", cmd, inputs)

    return True


class GlassHandler():
    """Handle glass restoration during Zemax import.

    This class handles the GCAT and GLAS commands found in .zmx files. If the
    glass can be matched up with an existing :mod:`opticalglass` catalog, the
    glass is instantiated and entered into the model. If the glass cannot be
    found, a search for a .smx file of the same name as the .zmx file is made.
    If found, it is a JSON file with a dict that provides an eval() string to
    create an instance to replace the missing Zemax glass name. If this file
    isn't found, it is created and contains a JSON template of a dict that has
    the missing glass names as keys; the values are the number of times the
    glass occurs in the file. Thes values should be replaced with the desired
    eval() string to create a replacement glass.
    """

    def __init__(self, filename):
        self.glass_catalogs = []
        self.glasses_not_found = {}
        self.filename = filename.with_suffix('.smx')
        self.glasses_not_found = self.load_replacements(self.filename)
        self.no_replacements = not self.glasses_not_found

    def load_replacements(self, filename):
        glasses_not_found = {}
        if filename.exists():
            with filename.open('r') as file:
                glasses_not_found = json.load(file)
        return glasses_not_found

    def save_replacements(self):
        with self.filename.open('w') as file:
            json.dump(self.glasses_not_found, file)

    def __call__(self, sm, cur, cmd, inputs):
        """ process GLAS command for fictitious, catalog glass or mirror"""
        if cmd == "GCAT":
            inputs = inputs.split()
            self.glass_catalogs = inputs
            return True
        elif cmd == "GLAS":
            g = sm.gaps[cur]
            inputs = inputs.split()
            name = inputs[0]
            medium = None
            if name == 'MIRROR':
                sm.ifcs[cur].interact_mode = 'reflect'
                g.medium = sm.gaps[cur-1].medium
                return True
            else:
                try:
                    medium = gfact.create_glass(name, self.glass_catalogs)
                except glasserror.GlassNotFoundError:
                    pass
                else:
                    g.medium = medium
                    return True

                medium = self.handle_glass_not_found(name)
                if medium is None:
                    medium = Medium(1.5, 'glass')
                g.medium = medium
                return True
        else:
            return False

    def handle_glass_not_found(self, name):
        """Record missing glasses or create new replacement glass instances."""

        """Import all supported glass catalogs."""
        from opticalglass.cdgm import CDGMGlass
        from opticalglass.hoya import HoyaGlass
        from opticalglass.ohara import OharaGlass
        from opticalglass.schott import SchottGlass

        if self.no_replacements:                # track the number of times
            if name in self.glasses_not_found:  # each missing glass is used
                self.glasses_not_found[name] += 1
            else:
                self.glasses_not_found[name] = 1
            return None

        else:  # create a new instance of the replacement glass
            if name in self.glasses_not_found:
                return eval(self.glasses_not_found[name])
            else:
                return None
