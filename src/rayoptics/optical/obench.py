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

import rayoptics.parax.specsheet as ss
from rayoptics.optical.opticalmodel import OpticalModel
from rayoptics.elem.profiles import (EvenPolynomial, RadialPolynomial)
from rayoptics.oprops import doe
from rayoptics.oprops.doe import DiffractiveElement
from rayoptics.raytr.opticalspec import WvlSpec
from rayoptics.util.misc_math import isanumber

from opticalglass import util

_track_contents = None


def read_obench_url(url, **kwargs):
    ''' given a url to a OpticalBench file, return an OpticalModel and info. '''
    global _track_contents
    url1 = url.replace('OpticalBench.htm#', '')
    url2 = url1.partition(',')[0]
    r = requests.get(url2, allow_redirects=True)

    apparent_encoding = r.apparent_encoding
    r.encoding = r.apparent_encoding
    inpt = r.text
    lines = inpt.splitlines()
    inpt = [l.split('\t') for l in lines]
    inpt_dict = {}
    for line in inpt:
        if line[0][0] == '[':
            # process new section header, initialize input list
            key = line[0][1:-1]
            inpt_dict[key] = []
        else:
            # add input to the currect section's list of inputs
            inpt_dict[key].append(line)

    opt_model = read_lens(inpt_dict, **kwargs)
    _track_contents['encoding'] = apparent_encoding

    return opt_model, _track_contents


def read_lens(inpts):
    global _track_contents
    def read_float(s):
        if s == 'Infinity':
            return float('inf')
        elif isanumber(s):
            return float(s)
        else:
            if s == 'undefined':
                return float('nan')
            elif s == 'AS':
                return 0.
            elif s == '':
                return 0.
            else:
                try:
                    return float(read_float(var_dists[s][0]))
                except:
                    return 0.

    _track_contents = util.Counter()
    constants_inpt = inpts['constants']
    constants = {c_item[0]: c_item[1:] for c_item in constants_inpt}
    var_dists_inpt = inpts['variable distances']
    var_dists = {var_dist[0]: var_dist[1:] for var_dist in var_dists_inpt}

    thi_obj = 0.
    if 'd0' in var_dists:
        thi_obj = read_float(var_dists['d0'][0])
        if thi_obj == float('inf'):
            thi_obj = 1.0e10
    conj_type = 'finite'
    if thi_obj > 1.0e8:
        conj_type = 'infinite'
    _track_contents['conj type'] = conj_type

    specsheet = ss.create_specsheet(conj_type)

    imager_inputs = dict(specsheet.imager_inputs)
    if conj_type == 'finite':
        imager_inputs['m'] = read_float(var_dists['Magnification'][0])
        specsheet.frozen_imager_inputs = [False]*5
    else:  # conj_type == 'infinite'
        imager_inputs['s'] = -float('inf')
        efl = read_float(var_dists['Focal Length'][0])
        if efl != 0:
            imager_inputs['f'] = efl
        specsheet.frozen_imager_inputs = [True, True, True, True, False]

    etendue_inputs = specsheet.etendue_inputs
    ape_key = ('aperture', 'image', 'f/#')
    ape_value = read_float(var_dists['F-Number'][0])
    etendue_inputs[ape_key[0]][ape_key[1]][ape_key[2]] = ape_value

    fld_key = ('field', 'image', 'height')
    fld_value = read_float(var_dists['Image Height'][0])/2
    etendue_inputs[fld_key[0]][fld_key[1]][fld_key[2]] = fld_value
    specsheet.generate_from_inputs(imager_inputs, etendue_inputs)

    opt_model = OpticalModel(do_init=True, radius_mode=True, specsheet=specsheet)
    sm = opt_model['sm']
    sm.do_apertures = False
    sm.gaps[0].thi = thi_obj

    osp = opt_model['osp']
    osp['fov'].is_relative = True
    osp['fov'].set_from_list([0., .707, 1.])
    osp['wvls'] = WvlSpec(wlwts=[('F', .5), ('d', 1.), ('C', .5)], ref_wl=1)

    if 'lens data' in inpts:
        input_lines = inpts['lens data']
        _track_contents['# surfs'] = len(input_lines)
        for line in input_lines:
            inpt = []
            inpt.append(read_float(line[1]))  # radius
            inpt.append(read_float(line[2]))  # thi
            if line[3] == '':
                inpt.append('')
                inpt.append('')
            else:
                inpt.append(read_float(line[3]))  # nd
                inpt.append(read_float(line[5]))  # vd
            inpt.append(read_float(line[4])/2)  # sd
            sm.add_surface(inpt)
            if line[1] == 'AS':
                sm.set_stop()
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
            idx = int(line[0])
            sm.ifcs[idx].profile = asp

    if 'diffractive data' in inpts:
        input_lines = inpts['diffractive data']
        _track_contents['# doe'] = len(input_lines)
        for line in input_lines:
            coefs = [read_float(item) for item in line[3:]]
            dif_elem = DiffractiveElement(coefficients=coefs,
                                          ref_wl=read_float(line[1]),
                                          order=read_float(line[2]),
                                          phase_fct=doe.radial_phase_fct)
            idx = int(line[0])
            sm.ifcs[idx].phase_element = dif_elem
    if 'descriptive data' in inpts:
        input_lines = inpts['descriptive data']
        descripts = {input_line[0]: input_line[1:] for input_line in input_lines}
        if 'title' in descripts:
            opt_model['sys'].title = descripts['title'][0]
    
    opt_model.update_model()
    return opt_model
