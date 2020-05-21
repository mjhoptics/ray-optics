#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
""" module to facilitate first order definition of an optical model

.. Created on Thu May 16 19:57:47 2019

.. codeauthor: Michael J. Hayford
"""
import math

from rayoptics.optical import firstorder
from rayoptics.optical.idealimager import IdealImager, ideal_imager_setup

from rayoptics.util import dict2d

from rayoptics.optical import etendue
from rayoptics.optical.etendue import (obj_img_set, fld_ape_set,
                                       fld_labels, ap_labels,
                                       create_etendue_dict)


def create_specsheet(conjugate_type, **inputs):
    if conjugate_type == 'finite':
        # setup finite conjugate defaults
        fev = create_etendue_dict()
        fev['field']['object'] = dict([(fld_labels[0], None)])
        fev['aperture']['object'] = dict([(ap_labels[2], None),
                                          (ap_labels[1], None)])
        fev['field']['image'] = dict([(fld_labels[0], None)])
        fev['aperture']['image'] = dict([(ap_labels[2], None),
                                         (ap_labels[1], None)])
        fss = SpecSheet('finite', etendue_values=fev, **inputs)
        return fss

    elif conjugate_type == 'infinite':
        # setup infinite conjugate defaults
        imager_inputs = {'s': -math.inf}
        imager = IdealImager(None, -math.inf, None, None, None)

        iev = create_etendue_dict()
        iev['field']['object'] = dict([(fld_labels[1], None)])
        iev['aperture']['object'] = dict([(ap_labels[0], None)])
        iev['field']['image'] = dict([(fld_labels[0], None)])
        iev['aperture']['image'] = dict([(ap_labels[2], None),
                                         (ap_labels[1], None)])

        ifss = SpecSheet('infinite', imager=imager,
                         imager_inputs=imager_inputs,
                         frozen_imager_inputs=[True, True, True, True, False],
                         etendue_values=iev, **inputs)
        return ifss

    else:
        print('create_specsheet: conjugate_type not recognized',
              conjugate_type)
        return None


def create_specsheets():
    specsheets = {}

    # setup finite conjugate defaults
    specsheets['finite'] = create_specsheet('finite')

    # setup infinite conjugate defaults
    specsheets['infinite'] = create_specsheet('infinite')

    return specsheets


def create_specsheet_from_model(opt_model, specsheets=None):
    if specsheets is None:
        specsheets = create_specsheets()

    specsheet = opt_model.specsheet
    if specsheet is None:
        conj_type = 'finite'
        if opt_model.seq_model.gaps[0].thi > 10e8:
            conj_type = 'infinite'
        specsheet = specsheets[conj_type]
    firstorder.specsheet_from_parax_data(opt_model, specsheet)
    opt_model.specsheet = specsheet
    return specsheet


class SpecSheet():
    """ First order optical specification for :class:`~.OpticalModel`

    Attributes:
        conjugate_type: one of `infinite`, `finite`
        imager: instance of `IdealImager`
        imager_inputs: dict of inputs to `ideal_imager_setup`
        frozen_imager_inputs: list of booleans, if True the parameter is frozen
        etendue_inputs: field and aperture inputs used to define the etendue
        etendue_values: dict2D of aperture/field vs object/image
        partitions: 'imager', 'field', and 'aperture'; number of items in each
    """

    def __init__(self, conjugate_type,
                 imager=None, imager_inputs=None, frozen_imager_inputs=None,
                 etendue_inputs=None, etendue_values=None):
        self.conjugate_type = conjugate_type

        if imager is None:
            imager = IdealImager(None, None, None, None, None)
        self.imager = imager

        self.imager_inputs = imager_inputs if imager_inputs else {}

        self.frozen_imager_inputs = (frozen_imager_inputs
                                     if frozen_imager_inputs
                                     else [False]*5)

        self.etendue_inputs = (etendue_inputs if etendue_inputs
                               else create_etendue_dict())
        self.etendue_values = (etendue_values if etendue_values
                               else create_etendue_dict())

        self.partition_defined()

    def __json_encode__(self):
        attrs = dict(vars(self))
        if hasattr(self, 'partitions'):
            del attrs['partitions']
        return attrs

    def __str__(self):
        return ("{!s} conjugates:\nimager: {}\n"
                "imager inputs: {}\n"
                "frozen imager inputs: {}\n"
                "etendue inputs:\n"
                "  field:    {}\n"
                "  aperture: {}\n"
                "etendue values:\n"
                "  field:    {}\n"
                "  aperture:\n"
                "    object: {}\n"
                "    image:  {}"
                .format(self.conjugate_type,
                        self.imager,
                        self.imager_inputs,
                        self.frozen_imager_inputs,
                        self.etendue_inputs['field'],
                        self.etendue_inputs['aperture'],
                        self.etendue_values['field'],
                        self.etendue_values['aperture']['object'],
                        self.etendue_values['aperture']['image']))

    def __repr__(self):
        return ("{!s}({!s}, imager={},"
                "imager_inputs={},"
                "frozen_imager_inputs={},"
                "etendue_inputs={},"
                "etendue_values={})".format(type(self).__name__,
                                            repr(self.conjugate_type),
                                            repr(self.imager),
                                            repr(self.imager_inputs),
                                            repr(self.frozen_imager_inputs),
                                            repr(self.etendue_inputs),
                                            repr(self.etendue_values)))

    def sync_to_restore(self, opt_model):
        # imager is exported as a list. convert back to an IdealImager
        self.imager = IdealImager(*self.imager)

    def imager_defined(self):
        """True if the imager is completely specified. """
        if self.conjugate_type == 'finite':
            imager_defined = 'm' if self.imager.m is not None else False
        else:
            imager_defined = 'f' if self.imager.f is not None else False
        return imager_defined

    def partition_defined(self):
        """ which partition defines the imager or None """
        num_imager_inputs = len(self.imager_inputs)
        li = dict2d.num_items_by_type(self.etendue_inputs,
                                      fld_ape_set, obj_img_set)
        num_field_inputs = li['field']
        num_aperture_inputs = li['aperture']
        partitions = {'imager': num_imager_inputs,
                      'field': num_field_inputs,
                      'aperture': num_aperture_inputs}
        self.partitions = partitions
        max_partition = max(partitions, key=partitions.get)
        max_num_inputs = partitions[max_partition]
        return max_partition if max_num_inputs == 2 else None, max_num_inputs

    def generate_from_inputs(self, imgr_inputs, etendue_inputs):
        """ compute imager and etendue values given input dicts """
        max_partition, max_num_inputs = self.partition_defined()
        num_imager_inputs = self.partitions['imager']
        num_field_inputs = self.partitions['field']
        num_aperture_inputs = self.partitions['aperture']

        conj_type = self.conjugate_type

        imager_inputs = {}
        if max_num_inputs <= 1:
            # fill in imager_inputs with any previous calculations for m or f
            if conj_type == 'finite':
                if num_imager_inputs < 2 and self.imager.m is not None:
                    imager_inputs['m'] = self.imager.m
            else:
                if num_imager_inputs < 2 and self.imager.f is not None:
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

        for fa_key, fa_value in etendue_inputs.items():
            for oi_key, oi_value in fa_value.items():
                if conj_type == 'finite':
                    conj = 'finite'
                elif conj_type == 'afocal':
                    conj = 'infinite'
                elif conj_type == 'infinite':
                    conj = 'finite' if oi_key == 'image' else conj_type
                etendue.fill_in_etendue_data(conj, imager, fa_key,
                                             etendue_inputs[fa_key][oi_key],
                                             etendue_values[fa_key][oi_key])

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
