#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Functions to read a CODE V .seq file and populate a sequential model

Created on Tue Jan 16 10:14:12 2018

@author: Michael J. Hayford
"""
import logging

from . import tla
from . import reader as cvr


_tla = tla.MapTLA()


def read_lens(ldm, filename):
    logging.basicConfig(filename='cv_cmd_proc.log',
                        filemode='w',
                        level=logging.DEBUG)
    cmds = cvr.read_seq_file(filename)
    for i, c in enumerate(cmds):
        cmd_fct, tla, qlist, dlist = process_command(c)
        if cmd_fct:
            eval_str = cmd_fct + '(ldm, tla, qlist, dlist)'
            eval(eval_str)
        else:
            logging.info('Line %d: Command %s not supported', i+1, c[0])

    ldm.update_model()


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


def wvl_spec_data(ldm, tla, qlist, dlist):
    if tla == "WL":
        ldm.global_spec.spectral_region.wavelengths = dlist
    elif tla == "WTW":
        ldm.global_spec.spectral_region.spectral_wts = dlist
    elif tla == "REF":
        ldm.global_spec.spectral_region.reference_wvl = dlist[0]-1
    elif tla == "CWL":
        ldm.global_spec.spectral_region.coating_wvl = dlist


def pupil_spec_data(ldm, tla, qlist, dlist):
    ldm.global_spec.pupil.type = tla
    ldm.global_spec.pupil.value = dlist[0]
    logging.debug("pupil_spec_data: %s %f", tla, dlist[0])


def field_spec_data(ldm, tla, qlist, dlist):
    if tla == "WID":
        ldm.global_spec.field_of_view.wide_angle = dlist[0]
    else:
        ldm.global_spec.field_of_view.update_fields_cv_input(tla, dlist)

    log_cmd("field_spec_data", tla, qlist, dlist)


def spec_data(ldm, tla, qlist, dlist):
    if tla == "LEN":
        ldm.reset()
    elif tla == "RDM":
        ldm.radius_mode = dlist[0]
    elif tla == "TIT":
        ldm.global_spec.specs.title = dlist[0]
    elif tla == "INI":
        ldm.global_spec.specs.initials = dlist[0]
    elif tla == "DIM":
        ldm.global_spec.specs.dimensions = dlist[0]
    elif tla == "CA":
        ldm.global_spec.specs.aperture_override = ''
    elif tla == "TEM":
        ldm.global_spec.specs.temperature = dlist[0]
    elif tla == "PRE":
        ldm.global_spec.specs.pressure = dlist[0]
    log_cmd("spec_data", tla, qlist, dlist)


def get_index_qualifier(ldm, qtype, qlist):
    def num_or_alpha(idx):
        if idx.isdigit():
            return int(idx)
        elif idx.isalpha():
            if qtype == 'L':
                return idx
            if idx == 'O':
                return 0
            elif idx == 'I':
                return ldm.get_num_surfaces()-1
            elif idx == 'S':
                return ldm.stop_surface
        else:
            return None
    for q in qlist:
        if q[0] == qtype:
            if q[1].find('..') > 0:
                i = q[1].find('..')
                return (num_or_alpha(q[1][:i]), num_or_alpha(q[1][i+2:]))
            else:
                return num_or_alpha(q[1]),


def surface_cmd(ldm, tla, qlist, dlist):
    idx, = get_index_qualifier(ldm, 'S', qlist)
    ldm.update_surface_and_gap_from_cv_input(dlist, idx)


def surface_data(ldm, tla, qlist, dlist):
    idx = get_index_qualifier(ldm, 'S', qlist)
    if not idx:
        idx = ldm.cur_surface
    if tla == 'STO':
        ldm.stop_surface = idx
    elif tla == 'THI':
        ldm.gaps[idx].thi = dlist[0]
    elif tla == 'SPH':
        ldm.update_surface_profile_cv_input('Spherical', idx)
    elif tla == 'CON':
        ldm.update_surface_profile_cv_input('Conic', idx)
    elif tla == 'ASP':
        ldm.update_surface_profile_cv_input('EvenPolynomial', idx)
    log_cmd("surface_data", tla, qlist, dlist)


def profile_data(ldm, tla, qlist, dlist):
    idx = get_index_qualifier(ldm, 'S', qlist)
    if not idx:
        idx = ldm.cur_surface
    if tla == 'CUX':
        ldm.surfs[idx].profile.cv = dlist[0]
    elif tla == 'CUY':
        ldm.surfs[idx].profile.cv = dlist[0]
    elif tla == 'RDX':
        if dlist[0] != 0.0:
            ldm.surfs[idx].profile.cv = 1.0/dlist[0]
        else:
            ldm.surfs[idx].profile.cv = 0.0
    elif tla == 'RDY':
        if dlist[0] != 0.0:
            ldm.surfs[idx].profile.cv = 1.0/dlist[0]
        else:
            ldm.surfs[idx].profile.cv = 0.0
    elif tla == 'K':
        ldm.surfs[idx].profile.cc = dlist[0]
    elif tla == 'A':
        ldm.surfs[idx].profile.coef4 = dlist[0]
    elif tla == 'B':
        ldm.surfs[idx].profile.coef6 = dlist[0]
    elif tla == 'C':
        ldm.surfs[idx].profile.coef8 = dlist[0]
    elif tla == 'D':
        ldm.surfs[idx].profile.coef10 = dlist[0]
    elif tla == 'E':
        ldm.surfs[idx].profile.coef12 = dlist[0]
    elif tla == 'F':
        ldm.surfs[idx].profile.coef14 = dlist[0]
    elif tla == 'G':
        ldm.surfs[idx].profile.coef16 = dlist[0]
    elif tla == 'H':
        ldm.surfs[idx].profile.coef18 = dlist[0]
    elif tla == 'J':
        ldm.surfs[idx].profile.coef20 = dlist[0]
    log_cmd("profile_data", tla, qlist, dlist)


def aperture_data(ldm, tla, qlist, dlist):
    log_cmd("aperture_data", tla, qlist, dlist)


def decenter_data(ldm, tla, qlist, dlist):
    idx = get_index_qualifier(ldm, 'S', qlist)
    if not idx:
        idx = ldm.cur_surface
    decenter = ldm.update_surface_decenter_cv_input(idx)
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
        decenter.type = 'DAR'
    elif tla == 'BEN':
        decenter.type = 'BEN'
    elif tla == 'REV':
        decenter.type = 'REV'

    decenter.update()

    log_cmd("decenter_data", tla, qlist, dlist)
