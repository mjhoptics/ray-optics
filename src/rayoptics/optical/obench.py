#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Import files from `OpticalBenchHub` web page

    This module implements lens import from the `OpticalBenchHub <https://www.photonstophotos.net/GeneralTopics/Lenses/OpticalBench/OpticalBenchHub.htm>`_
    portion of Bill Claff's `PhotonsToPhotos <https://www.photonstophotos.net/>`_
    website.

    To import a file from the website, navigate to the lens you wish to import 
    and select the entire web address of the page. Paste this into the url
    argument of the :func:`~.read_obench_url` function.

.. Created on Sat Jul 24 21:34:49 2021

.. codeauthor: Michael J. Hayford
"""
import requests

from rayoptics.optical.opticalmodel import OpticalModel
from rayoptics.elem.profiles import (EvenPolynomial, RadialPolynomial)
from rayoptics.oprops import doe
from rayoptics.oprops.doe import DiffractiveElement
from rayoptics.raytr.opticalspec import WvlSpec
from rayoptics.util.misc_math import isanumber, is_kinda_big, is_fuzzy_zero

from opticalglass import util
import opticalglass.modelglass as mg

_track_contents = None

def read_obench_url(url, **kwargs) -> tuple["OpticalModel", tuple[dict, dict]]:
    ''' given a url to a OpticalBench file, return an OpticalModel and info. '''
    global _track_contents
    _track_contents = util.Counter()

    obench_input, obench_dict = read_url(url, **kwargs)

    opt_model = read_lens(obench_dict, **kwargs)

    opt_model.update_model()

    _track_contents['obench input'] = obench_input   # type: ignore
    return opt_model, (_track_contents, {})   # type: ignore


def read_url(url, **kwargs) -> tuple[list, dict]:
    ''' given a url to a OpticalBench file, return an OpticalModel and info. '''
    global _track_contents

    _track_contents = (util.Counter() if _track_contents is None 
                       else _track_contents)
    
    url1 = url.replace('OpticalBench.htm#', '')
    url2 = url1.partition(',')[0]
    r = requests.get(url2, allow_redirects=True)

    apparent_encoding = r.apparent_encoding
    r.encoding = r.apparent_encoding
    inpt = r.text

    lines = inpt.splitlines()
    inpt = [l.split('\t') for l in lines]
    obench_dict = {}
    for line in inpt:
        if line[0][0] == '[':
            # process new section header, initialize input list
            key = line[0][1:-1]
            obench_dict[key] = []
        else:
            # add input to the currect section's list of inputs
            obench_dict[key].append(line)

    _track_contents['encoding'] = apparent_encoding   # type: ignore
    _track_contents['obench db'] = obench_dict   # type: ignore

    return lines, obench_dict


def read_lens(inpts, opt_model=None) -> OpticalModel:
    global _track_contents
    def read_float(s):
        if s == 'Infinity':
            return float('inf')
        elif isanumber(s):
            return float(s)
        else:
            if s == 'undefined':
                return float('nan')
            elif s == 'AS': # aperture stop
                return 0.
            elif s == 'FS': # field stop
                return 0.
            elif s == 'CG': # cover glass
                return 0.
            elif s == '':
                return 0.
            else:
                try:
                    return float(read_float(var_dists[s][0]))
                except:
                    return 0.

    _track_contents = (util.Counter() if _track_contents is None 
                       else _track_contents)
    if "constants" in inpts:
        constants_inpt = inpts['constants']
        constants = {c_item[0]: c_item[1:] for c_item in constants_inpt}
    else:
        constants = {}
    if "variable distances" in inpts:
        var_dists_inpt = inpts['variable distances']
        var_dists = {var_dist[0]: var_dist[1:] for var_dist in var_dists_inpt}
    else:
        var_dists = {}

    thi_obj = 1.0e10
    if 'd0' in var_dists:
        thi_obj = read_float(var_dists['d0'][0])
        if thi_obj == float('inf'):
            thi_obj = 1.0e10

    conj_type = 'finite'
    if is_kinda_big(thi_obj):
        conj_type = 'infinite'
    _track_contents['conj type'] = conj_type

    if opt_model is None:
        opt_model = OpticalModel(do_init=True)
    opt_model.radius_mode = True

    sm = opt_model['seq_model']
    sm.do_apertures = False
    sm.gaps[0].thi = thi_obj

    if 'ObjectGlass' in constants:
        nd_obj = read_float(constants['ObjectGlass'][0])
        vd_obj = read_float(constants['ObjectGlass'][1])
        sm.gaps[0].medium = mg.ModelGlass(nd_obj, vd_obj, 'ObjectGlass')

    osp = opt_model['optical_spec']
    if 'F-Number' in var_dists:
        osp['pupil'].key = ('image', 'f/#')
        osp['pupil'].value = read_float(var_dists['F-Number'][0])
    elif 'NA' in var_dists:
        osp['pupil'].key = ('object', 'NA')
        osp['pupil'].value = read_float(var_dists['NA'][0])

    if 'Image Height' in var_dists:
        img_ht = read_float(var_dists['Image Height'][0])/2
    if 'Angle of View' in var_dists:
        angle_of_view = read_float(var_dists['Angle of View'][0])
        osp['fov'].is_wide_angle = True if angle_of_view/2 > 45. else False
        osp['fov'].key = ('image', 'real height')
        osp['fov'].value = img_ht
    if 'Magnification' in var_dists:
        mag = read_float(var_dists['Magnification'][0])
        if not is_fuzzy_zero(mag):
            osp['fov'].key = ('object', 'height')
            osp['fov'].value = img_ht/mag
    osp['fov'].is_relative = True
    osp['fov'].set_from_list([0., .707, 1.])

    osp['wvls'] = WvlSpec(wlwts=[('F', .5), ('d', 1.), ('C', .5)], ref_wl=1)

    surfaces = {}   # map surface ID to surface num, surface ID is not always a sequence number
    if 'lens data' in inpts:
        input_lines = inpts['lens data']
        _track_contents['# surfs'] = len(input_lines)
        for line in input_lines:
            surf_id    = line[0]
            radius     = read_float(line[1])
            thickness  = read_float(line[2])
            nd         = line[3]
            vd         = line[5] if len(line) > 5 else ''
            diam       = read_float(line[4])
            glass_name = line[6] if len(line) > 6 else ''
            glass_make = line[7] if len(line) > 7 else ''
            inpt = [radius, thickness]
            # --- Determine optical material fields ---
            if nd == '':
                # Case 1: missing nd → push two blanks
                inpt.extend(['', ''])
            elif glass_name and glass_make:
                # Case 2: we have explicit glass name/make
                inpt.extend([glass_name,glass_make])
            else:
                # Case 3: numeric nd / vd
                inpt.append(read_float(nd))  # nd
                if vd != '':
                    inpt.append(read_float(vd))  # vd
            sm.add_surface(inpt, sd=diam/2)
            if line[1] == 'AS':
                sm.set_stop()
            surfaces[surf_id] = sm.cur_surface  # note surface num
    if 'aspherical data' in inpts:
        if 'AsphericalOddCount' in constants:
            typ = 'AsphericalOddCount'
        elif 'AsphericalA2' in constants:
            typ = 'AsphericalA2'
        else:
            typ = 'Aspherical'

        input_lines = inpts['aspherical data']
        _track_contents[typ] = len(input_lines)
        for line in input_lines:
            surf_id = line[0]
            if typ == 'AsphericalOddCount':
                asp_coefs = [read_float(item) for item in line[3:]]
                asp_coefs = [0., 0.] + asp_coefs
            elif typ == 'AsphericalA2':
                asp_coefs = [read_float(item) for item in line[3:]]
            else:
                asp_coefs = [read_float(item) for item in line[2:]]
                asp_coefs[0] = 0.
            if typ == 'AsphericalOddCount':
                asp = RadialPolynomial(r=read_float(line[1]),
                                       cc=read_float(line[2]),
                                       coefs=asp_coefs)
            else:
                asp = EvenPolynomial(r=read_float(line[1]),
                                     cc=read_float(line[2]),
                                     coefs=asp_coefs)
            idx = surfaces[surf_id] # lookup surface number using surface ID
            sm.ifcs[idx].profile = asp

    if 'diffractive data' in inpts:
        input_lines = inpts['diffractive data']
        _track_contents['# doe'] = len(input_lines)
        for line in input_lines:
            surf_id = line[0]
            coefs = [read_float(item) for item in line[3:]]
            dif_elem = DiffractiveElement(coefficients=coefs,
                                          ref_wl=read_float(line[1]),
                                          order=read_float(line[2]),
                                          phase_fct=doe.radial_phase_fct)
            idx = surfaces[surf_id] # lookup surface number using surface ID
            sm.ifcs[idx].phase_element = dif_elem
    if 'descriptive data' in inpts:
        input_lines = inpts['descriptive data']
        descripts = {input_line[0]: input_line[1:] 
                     for input_line in input_lines}

        if 'title' in descripts:
            opt_model['sys'].title = descripts['title'][0]

    return opt_model
