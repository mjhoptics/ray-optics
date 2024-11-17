#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Functions to read a CODE V .seq file and populate a sequential model

.. Created on Tue Jan 16 10:14:12 2018

.. codeauthor: Michael J. Hayford
"""
import logging
import math

from . import tla
from . import reader as cvr

from rayoptics.elem.surface import (DecenterData, Circular, Rectangular,
                                    Elliptical)
from rayoptics.elem import profiles
from rayoptics.oprops import doe
from rayoptics.seq.medium import GlassHandlerBase
from rayoptics.raytr.opticalspec import Field
from rayoptics.util.misc_math import isanumber

from opticalglass import util

from opticalglass import opticalmedium as om
from opticalglass import modelglass as mg

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
_fh = logging.FileHandler('cv_cmd_proc.log', mode='w', delay=True)
_fh.setLevel(logging.DEBUG)
logger.addHandler(_fh)

_tla = tla.MapTLA()

# Support for CODE V private catalog materials.
# Should this be in OpticalModel?
_reading_private_catalog = False
_private_catalog_wvls = None
_private_catalog_glasses = {}

_glass_handler = None
_track_contents = None


def fictitious_glass_decode(gc):
    ''' glass code parser, allowing for greater precision of n and v '''
    rindex_part = int(gc)
    magnitude = int(math.floor(math.log10(rindex_part))) + 1
    n = 1.0 + rindex_part/10**magnitude
    v = round(100.0*(gc - int(gc)), 6)
    return n, v


def read_lens(filename, **kwargs):
    """ given a CODE V .seq filename, return an OpticalModel
    
    Args:
        filename (pathlib.Path): a CODE V .seq file path
        kwargs (dict): keyword args passed to the reader functions

    Returns:
        an OpticalModel instance and a info tuple
    """
    import rayoptics.optical.opticalmodel as opticalmodel
    global _glass_handler, _track_contents
    global _reading_private_catalog
    _reading_private_catalog = False
    _track_contents = util.Counter()
    opt_model = opticalmodel.OpticalModel(do_init=False)
    _glass_handler = CVGlassHandler(filename)
    cmds = cvr.read_seq_file(filename)
    for i, c in enumerate(cmds):
        cmd_fct, tla, qlist, dlist = process_command(c)
        if cmd_fct:
            eval_str = cmd_fct + '(opt_model, tla, qlist, dlist)'
            eval(eval_str)
        else:
            logger.info('Line %d: Command %s not supported', i+1, c[0])

    post_process_input(opt_model, filename, **kwargs)
    _glass_handler.save_replacements()
    _track_contents.update(_glass_handler.track_contents)

    opt_model.update_model()

    info = _track_contents, _glass_handler.glasses_not_found
    return opt_model, info


def process_command(cmd):
    global _reading_private_catalog
    CmdFct, IndxQuals, DataType, Quals = range(4)
    tla = cmd[0][:3].upper()
    qlist = []
    dlist = []
    cmd_fct = None
    cmd_def = _tla.find(tla)
    if cmd_def:
        cmd_fct = cmd_def[CmdFct]
        iquals = cmd_def[IndxQuals]
        quals = cmd_def[Quals]
        data_type = cmd_def[DataType]
        data_found = False
        for t in cmd[1:]:
            if not data_found:
                if t in quals:
                    qlist.append((t,))
                elif t[:1] in iquals:
                    qlist.append((t[:1], t[1:]))
                else:
                    data_found = True
            if data_found:
                if data_type == 'String':
                    dlist.append(t)
                elif data_type == 'Double':
                    dlist.append(float(t))
                elif data_type == 'Integer':
                    dlist.append(int(t))
                elif data_type == 'Boolean':
                    if t[:1].upper() == 'N':
                        dlist.append(False)
                    else:
                        dlist.append(True)
    elif tla[:1] == 'S':
        cmd_fct = 'surface_cmd'
        qlist.append((tla[:1], tla[1:]))
        tla = 'S'
        dlist.append(float(cmd[1]))  # radius/curvature
        dlist.append(float(cmd[2]))  # thickness
        cmd_len = len(cmd)
        if cmd_len > 3:
            dlist.append(cmd[3])  # glass
        if cmd_len > 4:
            dlist.append(cmd[4])  # rmd
    elif _reading_private_catalog and isinstance(cmd[0], str):
        global _private_catalog_wvls, _private_catalog_glasses
        label = cmd[0]
        if isanumber(cmd[1]):
            for t in cmd[1:]:
                dlist.append(float(t))
            prv_glass = om.InterpolatedMedium(label, 
                                            wvls=_private_catalog_wvls, 
                                            rndx=dlist, cat='CV private catalog')
            _private_catalog_glasses[label] = prv_glass
        else:
            logger.debug(f"Unsupported PRV glass def: {cmd[0]} {cmd[1]}")

    return cmd_fct, tla, qlist, dlist


def log_cmd(label, tla, qlist, dlist):
    logger.debug("%s: %s %s %s", label, tla, str(qlist), str(dlist))


def post_process_input(opt_model, filename, **kwargs):
    global _track_contents
    sm = opt_model['seq_model']
    osp = opt_model['optical_spec']

    # retrieve image thickness and set defocus
    sm.z_dir.pop()
    sm.rndx.pop()
    gi = sm.gaps.pop()
    osp['focus'].focus_shift = gi.thi

    if opt_model.system_spec.title == '' and filename is not None:
        fname_full = filename.resolve()
        cat, fname = fname_full.parts[-2:]
        title = "{:s}: {:s}".format(cat, fname)
        opt_model.system_spec.title = title

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

    if _track_contents.get('# clear ap', 0) > 0:
        sm.do_apertures = False

    _track_contents['# wvls'] = len(osp['wvls'].wavelengths)
    
    fov = osp['fov']
    _track_contents['fov'] = fov.key
    _track_contents['# fields'] = len(fov.fields)
    max_fld, max_fld_idx = fov.max_field()
    fov.value = max_fld
    fov.is_wide_angle = fov.check_is_wide_angle()

def wvl_spec_data(opm, tla, qlist, dlist):
    osp = opm['optical_spec']
    if tla == "WL":
        osp['wvls'].wavelengths = dlist
        osp['wvls'].calc_colors()
    elif tla == "WTW":
        osp['wvls'].spectral_wts = dlist
    elif tla == "REF":
        osp['wvls'].reference_wvl = dlist[0]-1
    elif tla == "CWL":
        osp['wvls'].coating_wvl = dlist


def pupil_spec_data(opm, tla, qlist, dlist):
    pupil = opm['optical_spec']['pupil']
    if tla == "EPD":
        pupil.key = 'object', 'epd'
    elif tla == "NAO":
        pupil.key = 'object', 'NA'
    elif tla == "NA":
        pupil.key = 'image', 'NA'
    elif tla == "FNO":
        pupil.key = 'image', 'f/#'

    pupil.value = dlist[0]
    logger.debug("pupil_spec_data: %s %f", tla, dlist[0])


def field_spec_data(opm, tla, qlist, dlist):
    fov = opm['optical_spec']['fov']
    if tla == 'XOB' or tla == 'YOB':
        fov.key = 'object', 'height'
    elif tla == 'XAN' or tla == 'YAN':
        fov.key = 'object', 'angle'
    elif tla == 'XIM' or tla == 'YIM':
        fov.key = 'image', 'height'
    elif tla == 'XRI' or tla == 'YRI':
        fov.key = 'image', 'real height'
    elif tla == 'WID':
        fov.is_wide_angle = dlist[0]

    if len(fov.fields) != len(dlist):
        fov.fields = [Field() for f in range(len(dlist))]

    if tla[0] == 'V':
        attr = tla.lower()
    elif tla[0] == 'X' or tla[0] == 'Y':
        attr = tla[0].lower()
    elif tla == 'WTF':
        attr = 'wt'

    for i, f in enumerate(fov.fields):
        setattr(f, attr, dlist[i])

    log_cmd("field_spec_data", tla, qlist, dlist)


def spec_data(opm, tla, qlist, dlist):
    if tla == "LEN":
        pass
    elif tla == "RDM":
        if len(dlist) == 0:
            opm.radius_mode = True
        else:
            opm.radius_mode = dlist[0]
    elif tla == "TIT":
        opm['system_spec'].title = dlist[0]
    elif tla == "INI":
        opm['system_spec'].initials = dlist[0]
    elif tla == "DIM":
        dim = dlist[0].upper()
        if dim == 'M' or dim == 'MM':
            dim = 'mm'
        elif dim == 'C' or dim == 'CM':
            dim = 'cm'
        elif dim == 'I' or dim == 'IN':
            dim = 'inches'
        opm['system_spec'].dimensions = dim
    elif tla == "TEM":
        opm['system_spec'].temperature = dlist[0]
    elif tla == "PRE":
        opm['system_spec'].pressure = dlist[0]
    log_cmd("spec_data", tla, qlist, dlist)


def get_index_qualifier(seq_model, qtype, qlist):
    def num_or_alpha(idx):
        if idx.isdigit():
            return int(idx)
        elif idx.isalpha():
            if idx == 'O':
                return 0
            elif idx == 'I':
                return seq_model.get_num_surfaces()-1
            elif idx == 'S':
                return seq_model.stop_surface
        elif idx.isascii():
            if qtype == 'L':
                return idx
        else:
            return None
    for q in qlist:
        if q[0] == qtype:
            if q[1].find('..') > 0:
                i = q[1].find('..')
                return (num_or_alpha(q[1][:i]), num_or_alpha(q[1][i+2:]))
            else:
                return num_or_alpha(q[1]),


def surface_cmd(opt_model, tla, qlist, dlist):
    seq_model = opt_model.seq_model
    idx, = get_index_qualifier(seq_model, 'S', qlist)
    update_surface_and_gap(opt_model, dlist, idx)


def update_surface_and_gap(opt_model, dlist, idx=None):
    global _glass_handler
    seq_model = opt_model.seq_model
    s, g = seq_model.insert_surface_and_gap()

    if opt_model.radius_mode:
        if dlist[0] != 0.0:
            s.profile.cv = 1.0/dlist[0]
        else:
            s.profile.cv = 0.0
    else:
        s.profile.cv = dlist[0]

    if g:
        g.thi = dlist[1]

        if len(dlist) < 3:
            g.medium = om.Air()
        else:
            if dlist[2].upper() == 'REFL':
                s.interact_mode = 'reflect'
                g.medium = seq_model.gaps[seq_model.cur_surface-1].medium
            else:
                g.medium = _glass_handler.process_glass_data(dlist[2])

    else:
        # at image surface, apply defocus to previous thickness
        seq_model.gaps[idx-1].thi += dlist[1]


def private_catalog(optm, tla, qlist, dlist):
    global _reading_private_catalog, _private_catalog_wvls
    if tla == "PRV":
        _reading_private_catalog = True
    elif tla == "PWL":
        _private_catalog_wvls = dlist
    elif tla == "END":
        _reading_private_catalog = False
        _private_catalog_wvls = None

    log_cmd("private_catalog", tla, qlist, dlist)


def surface_data(opm, tla, qlist, dlist):
    seq_model = opm['seq_model']
    idx = get_index_qualifier(seq_model, 'S', qlist)
    if not idx:
        idx = seq_model.cur_surface
    if tla == 'SLB':
        seq_model.ifcs[idx].label = dlist[0]
    elif tla == 'STO':
        seq_model.stop_surface = idx
    elif tla == 'THI':
        seq_model.gaps[idx].thi = dlist[0]
    elif tla == 'SPH':
        update_surface_profile(seq_model, 'Spherical', idx)
    elif tla == 'CON':
        update_surface_profile(seq_model, 'Conic', idx)
    elif tla == 'ASP':
        update_surface_profile(seq_model, 'EvenPolynomial', idx)
    elif tla == 'XTO':
        update_surface_profile(seq_model, 'XToroid', idx)
    elif tla == 'YTO':
        update_surface_profile(seq_model, 'YToroid', idx)
    log_cmd("surface_data", tla, qlist, dlist)


def update_surface_profile(seq_model, profile_type, idx=None):
    if not isinstance(idx, int):
        idx = seq_model.cur_surface
    cur_profile = seq_model.ifcs[idx].profile
    new_profile = profiles.mutate_profile(cur_profile, profile_type)
    seq_model.ifcs[idx].profile = new_profile
    return seq_model.ifcs[idx].profile


def profile_data(opm, tla, qlist, dlist):
    seq_model = opm['seq_model']
    idx = get_index_qualifier(seq_model, 'S', qlist)
    if not idx:
        idx = seq_model.cur_surface
    if tla == 'CUX':
        seq_model.ifcs[idx].profile.cv = dlist[0]
    elif tla == 'CUY':
        seq_model.ifcs[idx].profile.cv = dlist[0]
    elif tla == 'RDX':
        if dlist[0] != 0.0:
            seq_model.ifcs[idx].profile.cv = 1.0/dlist[0]
        else:
            seq_model.ifcs[idx].profile.cv = 0.0
    elif tla == 'RDY':
        if dlist[0] != 0.0:
            seq_model.ifcs[idx].profile.cv = 1.0/dlist[0]
        else:
            seq_model.ifcs[idx].profile.cv = 0.0
    elif tla == 'K':
        seq_model.ifcs[idx].profile.cc = dlist[0]
    elif tla == 'A':
        seq_model.ifcs[idx].profile.coef4 = dlist[0]
    elif tla == 'B':
        seq_model.ifcs[idx].profile.coef6 = dlist[0]
    elif tla == 'C':
        seq_model.ifcs[idx].profile.coef8 = dlist[0]
    elif tla == 'D':
        seq_model.ifcs[idx].profile.coef10 = dlist[0]
    elif tla == 'E':
        seq_model.ifcs[idx].profile.coef12 = dlist[0]
    elif tla == 'F':
        seq_model.ifcs[idx].profile.coef14 = dlist[0]
    elif tla == 'G':
        seq_model.ifcs[idx].profile.coef16 = dlist[0]
    elif tla == 'H':
        seq_model.ifcs[idx].profile.coef18 = dlist[0]
    elif tla == 'J':
        seq_model.ifcs[idx].profile.coef20 = dlist[0]
    log_cmd("profile_data", tla, qlist, dlist)


def aperture_data(opm, tla, qlist, dlist):
    """ add aperture data, either creating a new aperture or modifying the last """
    seq_model = opm['seq_model']
    idx = get_index_qualifier(seq_model, 'S', qlist)
    lbl = get_index_qualifier(seq_model, 'L', qlist)
    if not idx:
        idx = seq_model.cur_surface
    ca_type = tla[0]
    data_type = tla[2]
    ifc = seq_model.ifcs[idx]
    ca_list = ifc.clear_apertures

    for q in qlist:
        if q[0] == 'EDG':
            ca_list = ifc.edge_apertures
        elif q[0] == 'HOL':
            log_cmd("aperture_data", tla, qlist, dlist)
#            ca_list = ifc.holes
            return
        elif q[0] == 'OBS':
            log_cmd("aperture_data", tla, qlist, dlist)
            return

    if len(ca_list) == 0 or ca_type != type(ca_list[-1]).__name__[0]:
        if ca_type == 'C':
            ca = Circular()
        elif ca_type == 'R':
            ca = Rectangular()
        elif ca_type == 'E':
            ca = Elliptical()
        ca_list.append(ca)
    else:
        ca = ca_list[-1]

    if data_type == 'R':
        ca.radius = dlist[0]
    elif data_type == 'X':
        ca.x_half_width = dlist[0]
    elif data_type == 'Y':
        ca.y_half_width = dlist[0]

    if ca:
        _track_contents['# clear ap'] += 1

    log_cmd("aperture_data", tla, qlist, dlist)


def aperture_data_general(opm, tla, qlist, dlist):
    """ handle the general aperture commands, add to end of list """
    seq_model = opm['seq_model']
    idx = get_index_qualifier(seq_model, 'S', qlist)
    if not idx:
        idx = seq_model.cur_surface
    ca_type = tla[0]

    if ca_type == 'C':
        ca = Circular(radius=dlist[0], x_offset=dlist[1], y_offset=dlist[2],
                      rotation=dlist[3])
    elif ca_type == 'R':
        ca = Rectangular(x_half_width=dlist[0], y_half_width=dlist[1],
                         x_offset=dlist[2], y_offset=dlist[3],
                         rotation=dlist[3])
    elif ca_type == 'E':
        ca = Elliptical(x_half_width=dlist[0], y_half_width=dlist[1],
                        x_offset=dlist[2], y_offset=dlist[3],
                        rotation=dlist[3])

    seq_model.ifcs[idx].clear_apertures.append(ca)

    log_cmd("aperture_data_general", tla, qlist, dlist)


def aperture_offset(opm, tla, qlist, dlist):
    """ handle the aperture offset commands, assume last aperture in list """
    seq_model = opm['seq_model']
    idx = get_index_qualifier(seq_model, 'S', qlist)
    if not idx:
        idx = seq_model.cur_surface
    ca = seq_model.ifcs[idx].clear_apertures[-1]
    offset_type = tla[2]

    if offset_type == 'X':
        ca.x_offset = dlist[0]
    elif offset_type == 'Y':
        ca.y_offset = dlist[0]
    elif offset_type == 'O':
        ca.rotation = dlist[0]

    log_cmd("aperture_offset", tla, qlist, dlist)


def decenter_data(opm, tla, qlist, dlist):
    seq_model = opm['seq_model']
    idx = get_index_qualifier(seq_model, 'S', qlist)
    if not idx:
        idx = seq_model.cur_surface
    ifc = seq_model.ifcs[idx]

    if not ifc.decenter:
        ifc.decenter = DecenterData('decenter')

    decenter = ifc.decenter
    if tla == 'XDE':
        decenter.dec[0] = dlist[0]
    elif tla == 'YDE':
        decenter.dec[1] = dlist[0]
    elif tla == 'ZDE':
        decenter.dec[2] = dlist[0]
    elif tla == 'ADE':
        decenter.euler[0] = dlist[0]
    elif tla == 'BDE':
        decenter.euler[1] = dlist[0]
    elif tla == 'CDE':
        decenter.euler[2] = dlist[0]
    elif tla == 'DAR':
        decenter.dtype = 'dec and return'
    elif tla == 'BEN':
        decenter.dtype = 'bend'
    elif tla == 'REV':
        decenter.dtype = 'reverse'

    decenter.update()

    log_cmd("decenter_data", tla, qlist, dlist)


#  DIF DOE
#  HOR 1.0
#  HWL 587.5618; HCT R
#  HCO C1 -0.001807322521767816; HCC C1 100
def diffractive_optic(opm, tla, qlist, dlist):
    seq_model = opm['seq_model']
    idx = get_index_qualifier(seq_model, 'S', qlist)
    if not idx:
        idx = seq_model.cur_surface
    ifc = seq_model.ifcs[idx]

    if tla == "DIF":
        for q in qlist:
            if "DOE" == q[0]:
                ifc.phase_element = doe.DiffractiveElement()
    elif tla == "HOR":
        if hasattr(ifc, 'phase_element'):
            ifc.phase_element.order = dlist[0]
    elif tla == "HWL":
        if hasattr(ifc, 'phase_element'):
            ifc.phase_element.ref_wl = dlist[0]
    elif tla == "HCT":
        if hasattr(ifc, 'phase_element'):
            for q in qlist:
                if "R" == q[0]:
                    ifc.phase_element.phase_fct = doe.radial_phase_fct
    elif tla == "HCO":
        cidx = get_index_qualifier(seq_model, 'C', qlist)[0]
        if hasattr(ifc, 'phase_element'):
            num_coefs = len(ifc.phase_element.coefficients)
            coefs = ifc.phase_element.coefficients
            if cidx <= num_coefs:
                coefs[cidx-1] = dlist[0]
            elif cidx == num_coefs + 1:
                coefs.append(dlist[0])
            else:
                coefs.extend([0.]*(cidx - num_coefs))
                coefs[cidx-1] = dlist[0]

    log_cmd("diffractive_optic", tla, qlist, dlist)


class CVGlassHandler(GlassHandlerBase):
    """Handle glass restoration during CODEV import.

    This class relies on GlassHandlerBase to provide most of the functionality
    needed to find the requested glass or a substitute.
    """

    def process_glass_data(self, glass_data):
        if isanumber(glass_data):
            # process as a 6 digit code, no decimal point
            medium = self.find_6_digit_code(glass_data)
            if medium is not None:
                self.track_contents['6 digit code'] += 1
                return medium
            else:  # process as fictitious glass code
                nd, vd = fictitious_glass_decode(float(glass_data))
                medium = mg.ModelGlass(nd, vd, glass_data)
                self.track_contents['fictitious glass'] += 1
                return medium
        else:  # look for glass name and optional catalog
            name_cat = glass_data.split('_')
            if len(name_cat) == 2:
                name, catalog = name_cat
            elif len(name_cat) == 1:
                name, catalog = glass_data, None

            if catalog is not None:
                if catalog.upper() == 'SCHOTT' and name[:1].upper() == 'N':
                    name = name[:1]+'-'+name[1:]
                elif catalog.upper() == 'OHARA' and name[:1].upper() == 'S':
                    name = name[:1]+'-'+name[1:]
                    if not name[-2:].isdigit() and name[-1].isdigit():
                        name = name[:-1]+' '+name[-1]

            medium = self.find_glass(name, catalog, always=False)
            if medium:
                return medium
            else:  # name with no data. default to crown glass
                global _private_catalog_glasses
                if name in _private_catalog_glasses:
                    medium = _private_catalog_glasses[name]
                else:
                    medium = om.ConstantIndex(1.5, 'not '+name)
                return medium
