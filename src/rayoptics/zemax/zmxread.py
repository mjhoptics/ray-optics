#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2020 Michael J. Hayford
"""

.. Created on Fri Jul 31 15:40:21 2020

.. codeauthor: Michael J. Hayford
"""
import logging
import math
import requests

from rayoptics.elem.surface import (DecenterData, Circular, Rectangular,
                                    Elliptical)
from rayoptics.elem import profiles
from rayoptics.seq.medium import GlassHandlerBase
from rayoptics.raytr.opticalspec import Field
from rayoptics.util.misc_math import isanumber
import rayoptics.zemax.zmx2ro as zmx2ro
from rayoptics.oprops import doe
import rayoptics.oprops.thinlens as thinlens

from opticalglass import glassfactory as gfact
from opticalglass import modelglass as mg
from opticalglass import opticalmedium as om
from opticalglass import glasserror
from opticalglass import util

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
_fh = logging.FileHandler('zmx_read_lens.log', mode='w', delay=True)
_fh.setLevel(logging.INFO)
logger.addHandler(_fh)

_glass_handler = None
_cmd_not_handled = None
_track_contents = None


def read_lens_file(filename, **kwargs):
    """ given a Zemax .zmx filename, return an OpticalModel

    It appears that Zemax .zmx files are written in UTF-16 encoding. Test
    against what seem to be common encodings. If other encodings are used in
    'files in the wild', please add them to the list.
    
    Args:
        filename (pathlib.Path): a Zemax .zmx file path
        kwargs (dict): keyword args passed to the reader functions

    Returns:
        an OpticalModel instance and a info tuple
    """
    global _track_contents
    if 'encoding' in kwargs:
        encodings = kwargs['encoding']
        if isinstance(encodings, str):
            encodings = encodings,
    else:  # default list of encodings to try
        encodings = ['utf-16', 'utf-8', 'utf-8-sig', 'iso-8859-1']

    for decode in encodings:
        try:
            with filename.open(encoding=decode) as file:
                inpt = file.read()
        except UnicodeError:
            pass
        else:
            break

    opt_model, info = read_lens(filename, inpt, **kwargs)
    _track_contents['encoding'] = decode

    return opt_model, info


def read_lens_url(url, **kwargs):
    ''' given a url to a Zemax file, return an OpticalModel  '''
    global _track_contents
    r = requests.get(url, allow_redirects=True)

    apparent_encoding = r.apparent_encoding
    r.encoding = r.apparent_encoding
    inpt = r.text

    opt_model, info = read_lens(None, inpt, **kwargs)
    _track_contents['encoding'] = apparent_encoding

    return opt_model, info


def read_lens(filename, inpt, **kwargs):
    ''' given inpt str of a Zemax .zmx file, return an OpticalModel  '''
    import rayoptics.optical.opticalmodel as opticalmodel
    global _glass_handler, _cmd_not_handled, _track_contents
    _cmd_not_handled = util.Counter()
    _track_contents = util.Counter()

    # create an empty optical model; all surfaces will come from .zmx file
    opt_model = opticalmodel.OpticalModel(do_init=False)

    input_lines = inpt.splitlines()

    _glass_handler = ZmxGlassHandler(filename)

    for i, line in enumerate(input_lines):
        process_line(opt_model, line, i+1)

    post_process_input(opt_model, filename, **kwargs)
    _glass_handler.save_replacements()
    _track_contents.update(_glass_handler.track_contents)

    opt_model.update_model()

    info = _track_contents, _glass_handler.glasses_not_found
    return opt_model, info


def process_line(opt_model, line, line_no):
    global _glass_handler, _cmd_not_handled, _track_contents
    sm = opt_model['seq_model']
    osp = opt_model['optical_spec']
    cur = sm.cur_surface
    if not line.strip():
        return
    line = line.strip().split(" ", 1)
    cmd = line[0]
    inputs = len(line) == 2 and line[1] or ""
    if cmd == "UNIT":
        dim = inputs.split()[0]
        if dim == 'MM':
            dim = 'mm'
        elif dim == 'IN' or dim == 'INCH':
            dim = 'inches'
        opt_model['system_spec'].dimensions = dim
    elif cmd == "NAME":
        opt_model['system_spec'].title = inputs.strip("\"")
    elif cmd == "NOTE":
        opt_model.note = inputs.strip("\"")
    elif cmd == "VERS":
        _track_contents["VERS"] = inputs.strip("\"")
    elif cmd == "SURF":
        s, g = sm.insert_surface_and_gap()
        # set type to Standard, some files don't have a Type command
        s.z_type = 'STANDARD'
    elif cmd == "CURV":
        s = sm.ifcs[cur]
        if hasattr(s, 'profile'):
            s.profile.cv = float(inputs.split()[0])
    elif cmd == "DISZ":
        g = sm.gaps[cur]
        g.thi = float(inputs)

    elif _glass_handler(sm, cur, cmd, inputs):
        pass

    elif cmd == "STOP":
        sm.set_stop()

    elif cmd == "WAVM":  # WAVM 1 0.55000000000000004 1
        sr = osp.spectral_region
        inputs = inputs.split()
        new_wvl = float(inputs[1])*1e+3
        if new_wvl not in sr.wavelengths:
            sr.wavelengths.append(new_wvl)
            sr.spectral_wts.append(float(inputs[2]))  # needs check
    # WAVL 0.4861327 0.5875618 0.6562725
    # WWGT 1 1 1
    elif cmd == "WAVL":
        sr = osp.spectral_region
        sr.wavelengths = [float(i)*1e+3 for i in inputs.split() if i]
    elif cmd == "WWGT":
        sr = osp.spectral_region
        sr.spectral_wts = [float(i)*1e+3 for i in inputs.split() if i]

    elif pupil_data(opt_model, cmd, inputs):
        pass

    elif field_spec_data(opt_model, cmd, inputs):
        pass

    elif handle_types_and_params(opt_model, cur, cmd, inputs):
        pass

    elif handle_aperture_data(opt_model, cur, cmd, inputs):
        pass

    elif cmd in ("OPDX",  # opd
                 "RAIM",  # ray aiming
                 "CONF",  # configurations
                 "PUPD",  # pupil
                 "EFFL",  # focal lengths
                 "MODE",  # mode
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
                 "ZVAN", "WWGN",
                 "WAVN", "MNCA", "MNEA",
                 "MNCG", "MNEG", "MXCA", "MXCG", "RGLA", "TRAC",
                 "TCMM", "FLOA", "PMAG", "TOTR", "SLAB",
                 "POPS", "COMM", "PZUP", "LANG", "FIMP", "COAT",
                 ):
        logger.info('Line %d: Command %s not supported', line_no, cmd)
    else:
        # don't recognize this cmd, record # of times encountered
        _cmd_not_handled[cmd] += 1


def post_process_input(opt_model, filename, **kwargs):
    global _track_contents
    sm = opt_model['seq_model']
    sm.gaps.pop()
    sm.z_dir.pop()
    sm.rndx.pop()

    if opt_model['system_spec'].title == '' and filename is not None:
        fname_full = filename.resolve()
        cat, fname = fname_full.parts[-2:]
        title = "{:s}: {:s}".format(cat, fname)
        opt_model['system_spec'].title = title

    conj_type = 'finite'
    if math.isinf(sm.gaps[0].thi):
        sm.gaps[0].thi = 1e10
        conj_type = 'infinite'
    _track_contents['conj type'] = conj_type

    sm.ifcs[0].label = 'Obj'
    sm.ifcs[0].interact_mode = 'dummy'
    sm.ifcs[-1].label = 'Img'
    sm.ifcs[-1].interact_mode = 'dummy'
    _track_contents['# surfs'] = len(sm.ifcs)

    # if DIAM records, turn off sm aperture setting
    if _track_contents.get('# clear ap', 0) > 0:
        sm.do_apertures = False

    do_post_processing = kwargs.get('do_postprocess', False)
    if do_post_processing:  # everything is on by default
        if kwargs.get('do_bend', True):
            zmx2ro.apply_fct_to_sm(opt_model, zmx2ro.convert_to_bend)
        if kwargs.get('do_dar', True):
            zmx2ro.apply_fct_to_sm(opt_model, zmx2ro.convert_to_dar)
        if kwargs.get('do_remove_null_sg', True):
            zmx2ro.apply_fct_to_sm(opt_model, zmx2ro.remove_null_sg)
        if kwargs.get('do_collapse_cb', True):
            zmx2ro.apply_fct_to_sm(opt_model, zmx2ro.collapse_coordbrk)

    osp = opt_model['optical_spec']
    sr = osp['wvls']
    if len(sr.wavelengths) > 1:
        if sr.wavelengths[-1] == 550.0:
            sr.wavelengths.pop()
            sr.spectral_wts.pop()
    sr.reference_wvl = len(sr.wavelengths)//2
    _track_contents['# wvls'] = len(sr.wavelengths)

    fov = osp['fov']
    _track_contents['fov'] = fov.key

    max_fld, max_fld_idx = fov.max_field()
    fov.value = max_fld
    fov.fields = [f for f in fov.fields[:max_fld_idx+1]]
    _track_contents['# fields'] = len(fov.fields)
    # switch vignetting definition to asymmetric vly, vuy style
    # need to verify this is how this works
    for f in fov.fields:
        # I think this is probably "has one has all", but we'll test all TBS
        if hasattr(f, 'vcx') and hasattr(f, 'vdx'):
            f.vlx = f.vcx + f.vdx
            f.vux = f.vcx - f.vdx
        if hasattr(f, 'vcy') and hasattr(f, 'vdy'):
            f.vly = f.vcy + f.vdy
            f.vuy = f.vcy - f.vdy
    fov.is_wide_angle = fov.check_is_wide_angle()


def log_cmd(label, cmd, inputs):
    logger.debug("%s: %s %s", label, cmd, str(inputs))


def handle_types_and_params(optm, cur, cmd, inputs):
    global _track_contents
    if cmd == "TYPE":
        ifc = optm.seq_model.ifcs[cur]
        typ = inputs.split()[0]
        # useful to remember the Type of Zemax surface
        ifc.z_type = typ
        _track_contents[typ] += 1
        if typ == 'EVENASPH':
            cur_profile = ifc.profile
            new_profile = profiles.mutate_profile(cur_profile,
                                                  'EvenPolynomial')
            ifc.profile = new_profile
        elif typ == 'TOROIDAL':
            cur_profile = ifc.profile
            new_profile = profiles.mutate_profile(cur_profile,
                                                  'YToroid')
            ifc.profile = new_profile
        elif typ == 'XOSPHERE':
            cur_profile = ifc.profile
            new_profile = profiles.mutate_profile(cur_profile,
                                                  'RadialPolynomial')
            ifc.profile = new_profile
        elif typ == 'COORDBRK':
            ifc.interact_mode = 'dummy'
            ifc.decenter = DecenterData('decenter')
        elif typ == 'PARAXIAL':
            ifc = thinlens.ThinLens()
            ifc.z_type = typ
            optm.seq_model.ifcs[cur] = ifc
        elif typ == 'DGRATING':
            ifc.phase_element = doe.DiffractionGrating()
            ifc.z_type = typ
            optm.seq_model.ifcs[cur] = ifc
    elif cmd == "CONI":
        _track_contents["CONI"] += 1
        ifc = optm.seq_model.ifcs[cur]
        cur_profile = ifc.profile
        if not hasattr(cur_profile, 'cc'):
            ifc.profile = profiles.mutate_profile(cur_profile, 'Conic')
        ifc.profile.cc = float(inputs.split()[0])
    elif cmd == "PARM":
        ifc = optm.seq_model.ifcs[cur]
        i, param_val = inputs.split()
        i = int(i)
        param_val = float(param_val)
        if ifc.z_type == 'COORDBRK':
            if i == 1:
                ifc.decenter.dec[0] = param_val
            elif i == 2:
                ifc.decenter.dec[1] = param_val
            elif i == 3:
                ifc.decenter.euler[0] = param_val
            elif i == 4:
                ifc.decenter.euler[1] = param_val
            elif i == 5:
                ifc.decenter.euler[2] = param_val
            elif i == 6:
                if param_val != 0:
                    ifc.decenter.dtype = 'reverse'
            ifc.decenter.update()
        elif ifc.z_type == 'DGRATING':
            if i == 1:
                ifc.phase_element.grating_freq_um = param_val
            elif i == 2:
                ifc.phase_element.order = param_val
        elif ifc.z_type == 'EVENASPH':
            ifc.profile.coefs[i-1] = param_val
        elif ifc.z_type == 'PARAXIAL':
            if i == 1:
                ifc.optical_power = 1/param_val
        elif ifc.z_type == 'TOROIDAL':
            if i == 1:
                ifc.profile.rR = param_val
            elif i > 1:
                ifc.profile.coefs.append(param_val)
    elif cmd == "XDAT":
        ifc = optm.seq_model.ifcs[cur]
        inputs = inputs.split()
        i = int(inputs[0])
        param_val = float(inputs[1])
        if ifc.z_type == 'XOSPHERE':
            if i == 1:
                num_terms = param_val
                ifc.profile.coefs = []
            elif i == 2:
                normalizing_radius = param_val
                if normalizing_radius != 1.0:
                    logger.info('Normalizing radius not supported on extended surfaces')
            elif i >= 3:
                ifc.profile.coefs.append(param_val)
    else:
        return False
    return True


def handle_aperture_data(optm, cur, cmd, inputs):
    # DIAM 7.5 1 0 0 1 ""
    # FLAP 0 7.5 0
    # CLAP 0 25.399999999999999 0
    # OBDC 0.000000000000E+00 1.906000000000E+02
    global _track_contents
    sm = optm.seq_model
    items = inputs.split()
    if cmd == "DIAM":
        ifc = sm.ifcs[cur]
        ca_val = float(items[0])
        if ca_val == 0.0:
            ca_val = 1.0
            logger.info(f"Surf {cur}: zero value on DIAM input.")
        ca_type = int(items[1])
        if hasattr(ifc, 'clear_apertures'):
            ca_list = ifc.clear_apertures
            if len(ca_list) == 0:
                ca = None
                if ca_type == 0:
                    ca = Circular()
                elif ca_type == 1:
                    ca = Circular()
                elif ca_type == 4:
                    ca = Rectangular()
                    _track_contents['non_circular_ca_type'] += 1
                elif ca_type == 6:
                    ca = Elliptical()
                    _track_contents['non_circular_ca_type'] += 1
                else:
                    _track_contents['ca_type_not_recognized'] += 1
                    # print('ca_type', cur, ca_type, items[1])
                    return True
    
                if ca:
                    _track_contents['# clear ap'] += 1
                    ca_list.append(ca)
            else:
                ca = ca_list[-1]
    
            ca.radius = ca_val
        ifc.set_max_aperture(ca_val)
    elif cmd == "OBDC":
        # appears to be aperture offsets, x and y
        ifc = sm.ifcs[cur]
        ca = ifc.clear_apertures[0]
        ca.x_offset = float(items[0])
        ca.y_offset = float(items[1])
    elif cmd == "FLAP":
        # Don't really understand how this is used...
        pass
    else:
        return False

    return True


def pupil_data(optm, cmd, inputs):
    # FNUM 2.1 0
    # OBNA 1.5E-1 0
    # ENPD 20
    global _track_contents
    pupil = optm.optical_spec.pupil
    if cmd == 'FNUM':
        pupil.key = 'image', 'f/#'
    elif cmd == 'OBNA':
        pupil.key = 'object', 'NA'
    elif cmd == 'ENPD':
        pupil.key = 'object', 'epd'
    else:
        return False

    _track_contents['pupil'] = pupil.key

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

    # older files (perhaps?)
    # XFLD 0 0 0
    # YFLD 0 35 50
    # FWGT 1 1 1
    global _track_contents

    fov = optm.optical_spec.field_of_view
    if cmd == 'XFLN' or cmd == 'YFLN' or cmd == 'XFLD' or cmd == 'YFLD':
        attr = cmd[0].lower()
    elif cmd == 'FTYP':
        ftyp = int(inputs.split()[0])
        _track_contents["FTYP"] = inputs
        if ftyp == 0:
            fov.key = 'object', 'angle'
        elif ftyp == 1:
            fov.key = 'object', 'height'
        elif ftyp == 2:
            fov.key = 'image', 'height'
        elif ftyp == 3:
            fov.key = 'image', 'real height'
        return True
    elif cmd == 'VDXN' or cmd == 'VDYN':
        attr = 'vd' + cmd[2].lower()
    elif cmd == 'VCXN' or cmd == 'VCYN':
        attr = 'vc' + cmd[2].lower()
    elif cmd == 'VANN':
        attr = 'van'
    elif cmd == 'FWGN' or cmd == 'FWGT':
        attr = 'wt'
    else:
        return False

    inputs = inputs.split()

    if len(fov.fields) < len(inputs):
        fov.fields += [Field() for f in range(len(fov.fields), len(inputs))]

    for i, inpt in enumerate(inputs):
        setattr(fov.fields[i], attr, float(inpt))

    log_cmd("field_spec_data", cmd, inputs)

    return True


class ZmxGlassHandler(GlassHandlerBase):
    """Handle glass restoration during Zemax zmx import.

    This class relies on GlassHandlerBase to provide most of the functionality
    needed to find the requested glass or a substitute.
    """

    def __call__(self, sm, cur, cmd, inputs):
        """ process GLAS command for fictitious, catalog glass or mirror"""

        if cmd == "GCAT":
            inputs = inputs.split()
            # Check catalog names, only add those we recognize
            for gc in inputs:
                try:
                    gfact.get_glass_catalog(gc)
                except glasserror.GlassCatalogNotFoundError:
                    continue
                else:
                    self.glass_catalogs.append(gc)
            # If no catalogs were recognized, use the default set
            if len(self.glass_catalogs) == 0:
                self.glass_catalogs = gfact._cat_names
            self.track_contents["GCAT"] = inputs
            return True
        elif cmd == "GLAS":
            g = sm.gaps[cur]
            inputs = inputs.split()
            name = inputs[0]
            medium = None
            if name == 'MIRROR':
                sm.ifcs[cur].interact_mode = 'reflect'
                g.medium = sm.gaps[cur-1].medium
                self.track_contents[name] += 1
                return True
            elif name == '___BLANK':
                nd = float(inputs[3])
                vd = float(inputs[4])
                if vd == 0:
                    # Zemax treats Vd=0 as constant index
                    g.medium = om.ConstantIndex(nd, f"n:{nd:.3f}")
                else:
                    g.medium = mg.ModelGlass(nd, vd, om.glass_encode(nd, vd))
                self.track_contents[name] += 1
                return True
            elif isanumber(name):
                # process as a 6 digit code, no decimal point
                m = self.find_6_digit_code(name)
                if m is not None:
                    g.medium = m
                    self.track_contents['6 digit code'] += 1
                    return True
            else:  # must be a glass type
                medium = self.find_glass(name, '')
                g.medium = medium
                return True

        else:
            return False
