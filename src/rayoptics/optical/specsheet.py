#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
""" module to facilitate first order definition of an optical model

.. Created on Thu May 16 19:57:47 2019

.. codeauthor: Michael J. Hayford
"""
from rayoptics.optical.idealimager import IdealImager, ideal_imager_setup

from rayoptics.util import dict2d
from rayoptics.util.dict2d import dict2D

from rayoptics.optical import etendue
from rayoptics.optical.etendue import (obj_img_set, fld_ape_set,
                                       fld_labels, ap_labels)


class SpecSheet():
    def __init__(self, conjugate_type, imager=None, imager_inputs=None,
                 etendue_inputs=None, etendue_values=None):
        self.conjugate_type = conjugate_type

        if imager is None:
            imager = IdealImager(None, None, None, None, None)
        self.imager = imager

        self.imager_inputs = imager_inputs if imager_inputs else {}

        self.etendue_inputs = (etendue_inputs if etendue_inputs
                               else dict2D(fld_ape_set, obj_img_set))
        self.etendue_values = (etendue_values if etendue_values
                               else dict2D(fld_ape_set, obj_img_set))

        if conjugate_type is 'infinite':
            self.etendue_values['field']['object'] = dict([
                    (fld_labels[1], None)])
            self.etendue_values['aperture']['object'] = dict([
                    (ap_labels[0], None)])
            self.etendue_values['field']['image'] = dict([
                    (fld_labels[0], None)])
            self.etendue_values['aperture']['image'] = dict([
                    (ap_labels[2], None), (ap_labels[1], None)])

        else:
            self.etendue_values['field']['object'] = dict([
                    (fld_labels[0], None)])
            self.etendue_values['aperture']['object'] = dict([
                    (ap_labels[2], None), (ap_labels[1], None)])
            self.etendue_values['field']['image'] = dict([
                    (fld_labels[0], None)])
            self.etendue_values['aperture']['image'] = dict([
                    (ap_labels[2], None), (ap_labels[1], None)])

    def generate_from_inputs(self, imgr_inputs, etendue_inputs):
        """ compute imager and etendue values given input dicts """
        conj_type = self.conjugate_type
        imager_inputs = {}
        # fill in imager_inputs with any previous calculations for m or f
        if conj_type == 'finite':
            if self.imager.m is not None:
                imager_inputs['m'] = self.imager.m
        else:
            if self.imager.f is not None:
                imager_inputs['f'] = self.imager.f

        # update imager_inputs with user entries
        imager_inputs.update(imgr_inputs)
        imager_inputs = {k: v for (k, v) in imager_inputs.items()
                         if v is not None}

        # calculate an ideal imager for imager_inputs
        imager = ideal_imager_setup(**imager_inputs)

        if conj_type == 'finite':
            imager_defined = True if imager.m is not None else False
        else:
            imager_defined = True if imager.f is not None else False

        etendue_values = self.etendue_values

        for fld_ape_key, fld_ape_value in etendue_inputs.items():
            for obj_img_key, input in fld_ape_value.items():
                for key in input:
                    etendue_values[fld_ape_key][obj_img_key][key] = \
                        etendue_inputs[fld_ape_key][obj_img_key][key]

        li = dict2d.num_items_by_type(etendue_inputs, fld_ape_set, obj_img_set)
        num_field_inputs = li['field']
        num_aperture_inputs = li['aperture']
        if imager_defined:
            if num_field_inputs >= 1 and num_aperture_inputs >= 1:
                # we have enough data to calculate all of the etendue grid
                ii = etendue.do_etendue_via_imager(conj_type, imager,
                                                   etendue_inputs,
                                                   etendue_values)

                if ii:
                    imager_inputs[ii[0]] = ii[1]
                    imager = ideal_imager_setup(**imager_inputs)
                    etendue.do_etendue_via_imager(conj_type, imager,
                                                  etendue_inputs,
                                                  etendue_values)
            elif num_field_inputs == 1:
                # we have enough data to calculate all of the etendue grid
                row = dict2d.row(etendue_inputs, 'field')
                obj_img_key = 'object' if len(row['object']) else 'image'
                etendue.do_field_via_imager(conj_type, imager, etendue_inputs,
                                            obj_img_key, etendue_values)
            elif num_aperture_inputs == 1:
                # we have enough data to calculate all of the etendue grid
                row = dict2d.row(etendue_inputs, 'aperture')
                obj_img_key = 'object' if len(row['object']) else 'image'
                etendue.do_aperture_via_imager(conj_type, imager,
                                               etendue_inputs, obj_img_key,
                                               etendue_values)
        else:  # imager not specified
            if num_field_inputs == 2 or num_aperture_inputs == 2:
                fld_ape_key = 'field' if num_field_inputs == 2 else 'aperture'
                # solve for imager
                ii = etendue.do_etendue_to_imager(fld_ape_key, etendue_inputs,
                                                  etendue_values)
                imager_inputs[ii[0]] = ii[1]
                imager = ideal_imager_setup(**imager_inputs)
                # update etendue grid
                etendue.do_etendue_via_imager(conj_type, imager,
                                              etendue_inputs, etendue_values)

        self.imager = imager
        self.etendue_values = etendue_values
        return imager, etendue_values
