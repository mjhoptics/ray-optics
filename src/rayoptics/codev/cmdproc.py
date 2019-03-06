#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Functions to read a CODE V .seq file and populate a sequential model

.. Created on Tue Jan 16 10:14:12 2018

.. codeauthor: Michael J. Hayford
"""
import logging

from . import tla
from . import reader as cvr

from rayoptics.optical.model_enums import PupilType
from rayoptics.optical.model_enums import DimensionType as dt
from rayoptics.optical.model_enums import DecenterType as dec
from rayoptics.optical.surface import (DecenterData, Circular, Rectangular,
                                       Elliptical)
from rayoptics.optical import profiles
from rayoptics.optical import medium as m
from rayoptics.util.misc_math import isanumber

from opticalglass import glassfactory as gfact
from opticalglass import glasserror as ge

_tla = tla.MapTLA()


def read_lens(optm, filename):
    logging.basicConfig(filename='cv_cmd_proc.log',
                        filemode='w',
                        level=logging.DEBUG)
    cmds = cvr.read_seq_file(filename)
    for i, c in enumerate(cmds):
        cmd_fct, tla, qlist, dlist = process_command(c)
        if cmd_fct:
            eval_str = cmd_fct + '(optm, tla, qlist, dlist)'
            eval(eval_str)
        else:
            logging.info('Line %d: Command %s not supported', i+1, c[0])

    optm.update_model()


def process_command(cmd):
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
                    dlist.append(eval("float("+t+")"))
                elif data_type == 'Integer':
                    dlist.append(eval("int("+t+")"))
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

    return cmd_fct, tla, qlist, dlist


def log_cmd(label, tla, qlist, dlist):
    logging.debug("%s: %s %s %s", label, tla, str(qlist), str(dlist))


def wvl_spec_data(optm, tla, qlist, dlist):
    osp = optm.optical_spec
    if tla == "WL":
        osp.spectral_region.wavelengths = dlist
        osp.spectral_region.calc_colors()
    elif tla == "WTW":
        osp.spectral_region.spectral_wts = dlist
    elif tla == "REF":
        osp.spectral_region.reference_wvl = dlist[0]-1
    elif tla == "CWL":
        osp.spectral_region.coating_wvl = dlist


def pupil_spec_data(optm, tla, qlist, dlist):
    osp = optm.optical_spec
    osp.pupil.pupil_type = PupilType[tla]
    osp.pupil.value = dlist[0]
    logging.debug("pupil_spec_data: %s %f", tla, dlist[0])


def field_spec_data(optm, tla, qlist, dlist):
    osp = optm.optical_spec
    if tla == "WID":
        osp.field_of_view.wide_angle = dlist[0]
    else:
        osp.field_of_view.update_fields_cv_input(tla, dlist)

    log_cmd("field_spec_data", tla, qlist, dlist)


def spec_data(optm, tla, qlist, dlist):
    if tla == "LEN":
        optm.reset()
    elif tla == "RDM":
        optm.radius_mode = dlist[0]
    elif tla == "TIT":
        optm.system_spec.title = dlist[0]
    elif tla == "INI":
        optm.system_spec.initials = dlist[0]
    elif tla == "DIM":
        dim = dlist[0].upper()
        if dim == 'M':
            dim = dt.MM
        elif dim == 'C':
            dim = dt.CM
        elif dim == 'I':
            dim = dt.IN
        optm.system_spec.dimensions = dim
    elif tla == "TEM":
        optm.system_spec.temperature = dlist[0]
    elif tla == "PRE":
        optm.system_spec.pressure = dlist[0]
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
    seq_model = opt_model.seq_model
    if isinstance(idx, int):
        s, g = seq_model.get_surface_and_gap(idx)
        seq_model.set_cur_surface(idx)
    else:
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
            g.medium = m.Air()
        elif dlist[2].upper() == 'REFL':
            s.refract_mode = 'REFL'
            g.medium = seq_model.gaps[seq_model.cur_surface-1].medium
        elif isanumber(dlist[2]):
            n, v = m.glass_decode(float(dlist[2]))
            g.medium = m.Glass(n, v, '')
        else:
            name, cat = dlist[2].split('_')
            if cat.upper() == 'SCHOTT' and name[:1].upper() == 'N':
                name = name[:1]+'-'+name[1:]
            try:
                g.medium = gfact.create_glass(name, cat)
            except ge.GlassNotFoundError as gerr:
                logging.info('%s glass data type %s not found',
                             gerr.catalog,
                             gerr.name)
                logging.info('Replacing material with air.')
                g.medium = m.Air()
    else:
        # at image surface, apply defocus to previous thickness
        seq_model.gaps[idx-1].thi += dlist[1]


def surface_data(optm, tla, qlist, dlist):
    seq_model = optm.seq_model
    idx = get_index_qualifier(seq_model, 'S', qlist)
    if not idx:
        idx = seq_model.cur_surface
    if tla == 'STO':
        seq_model.stop_surface = idx
    elif tla == 'THI':
        seq_model.gaps[idx].thi = dlist[0]
    elif tla == 'SPH':
        update_surface_profile(seq_model, 'Spherical', idx)
    elif tla == 'CON':
        update_surface_profile(seq_model, 'Conic', idx)
    elif tla == 'ASP':
        update_surface_profile(seq_model, 'EvenPolynomial', idx)
    log_cmd("surface_data", tla, qlist, dlist)


def update_surface_profile(seq_model, profile_type, idx=None):
    if not isinstance(idx, int):
        idx = seq_model.cur_surface
    cur_profile = seq_model.ifcs[idx].profile
    new_profile = profiles.mutate_profile(cur_profile, profile_type)
    seq_model.ifcs[idx].profile = new_profile
    return seq_model.ifcs[idx].profile


def profile_data(optm, tla, qlist, dlist):
    seq_model = optm.seq_model
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
    seq_model = opm.seq_model
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

    log_cmd("aperture_data", tla, qlist, dlist)


def aperture_data_general(opm, tla, qlist, dlist):
    """ handle the general aperture commands, add to end of list """
    seq_model = opm.seq_model
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
    seq_model = opm.seq_model
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


def decenter_data(optm, tla, qlist, dlist):
    seq_model = optm.seq_model
    idx = get_index_qualifier(seq_model, 'S', qlist)
    if not idx:
        idx = seq_model.cur_surface
    ifc = seq_model.ifcs[idx]

    if not ifc.decenter:
        ifc.decenter = DecenterData(dec.LOCAL)

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
        decenter.dtype = dec.DAR
    elif tla == 'BEN':
        decenter.dtype = dec.BEND
    elif tla == 'REV':
        decenter.dtype = dec.REV

    decenter.update()

    log_cmd("decenter_data", tla, qlist, dlist)
