#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
"""
Created on Fri Jul  6 07:24:40 2018

@author: Michael J. Hayford
"""

import numpy as np

from optical.model_constants import ht, slp, aoi
from util.misc_math import transpose


def compute_third_order(seq_model):
    """ Compute Seidel aberration coefficents. """
    wl = seq_model.optical_spec.spectral_region.reference_wvl
    n_before = seq_model.rndx[0][wl]
    ax_ray, pr_ray, fod = seq_model.optical_spec.parax_data
    opt_inv = fod.opt_inv
    opt_inv_sqr = opt_inv*opt_inv

    SI = 0.
    SII = 0.
    SIII = 0.
    SIV = 0.
    SV = 0.
    third_order = []
    third_order_asp = {}

    # Transfer from object
    p = 0

    for c in range(1, len(ax_ray)-1):
        n_after = seq_model.rndx[c][wl]
        cv = seq_model.ifcs[c].profile_cv()

        A = n_after * ax_ray[c][aoi]
        Abar = n_after * pr_ray[c][aoi]
        P = cv*(1./n_after - 1./n_before)
        delta_slp = ax_ray[c][slp]/n_after - ax_ray[p][slp]/n_before
        SIi = -A**2 * ax_ray[c][ht] * delta_slp
        SIIi = -A*Abar * ax_ray[c][ht] * delta_slp
        SIIIi = -Abar**2 * ax_ray[c][ht] * delta_slp
        SIVi = -opt_inv_sqr * P
        delta_n_sqr = 1./n_after**2 - 1./n_before**2
        SVi = -Abar*(Abar * Abar * delta_n_sqr * ax_ray[c][ht] -
                     (opt_inv + Abar * ax_ray[c][ht])*pr_ray[c][ht]*P)

        # handle case of aspheric profile
        if hasattr(seq_model.ifcs[c], 'profile'):
            to_asp = aspheric_seidel_contribution(seq_model, c)
            if to_asp:
                third_order_asp[c] = to_asp

        SI += SIi
        SII += SIIi
        SIII += SIIIi
        SIV += SIVi
        SV += SVi

        third_order.append([SIi, SIIi, SIIIi, SIVi, SVi])

        p = c
        n_before = n_after

    third_order_total = [sum(v) for v in transpose(third_order)]

    if len(third_order_asp) > 0:
        asp_vals = list(third_order_asp.values())
        sum_asp = [sum(v) for v in transpose(asp_vals)]
        third_order_total = [third_order_total[i]+sum_asp[i]
                             for i in range(len(sum_asp))]
#        third_order_total = np.add(third_order_total, sum_asp)

    third_order.append(third_order_total)
    return (third_order, third_order_total, third_order_asp)


def calc_4th_order_aspheric_term(p):
    G = 0.
    if type(p).__name__ == 'Conic':
        cv = p.cv
        cc = p.cc
        G = cc*cv**3/8.0
    elif type(p).__name__ == 'EvenPolynomial':
        cv = p.cv
        cc = p.cc
        G = cc*cv**3/8.0 + p.coef4

    return G


def aspheric_seidel_contribution(seq_model, i):
    def delta_E(z, y, u, n):
        return -z/(n*y*(y + z*u))
    ax_ray, pr_ray, fod = seq_model.optical_spec.parax_data
    z = -pr_ray[i][ht]/pr_ray[i][slp]
    e = fod.opt_inv*delta_E(z, ax_ray[i][ht], ax_ray[i][slp],
                            seq_model.central_rndx(i))
    G = calc_4th_order_aspheric_term(seq_model.ifcs[i].profile)
    if G == 0.0:
        return None
    delta_n = seq_model.central_rndx(i) - seq_model.central_rndx(i-1)
    SI_star = 8.0*G*delta_n*ax_ray[i][ht]**4
    SII_star = SI_star*e
    SIII_star = SI_star*e**2
    SV_star = SI_star*e**3
    return [SI_star, SII_star, SIII_star, 0., SV_star]


def seidel_to_wavefront(seidel, central_wvl):
    """ Convert Seidel coefficients to wavefront aberrations """
    SI, SII, SIII, SIV, SV = seidel
    W040 = 0.125*SI/central_wvl
    W131 = 0.5*SII/central_wvl
    W222 = 0.5*SIII/central_wvl
    W220 = 0.25*(SIV + SIII)/central_wvl
    W311 = 0.5*SV/central_wvl
    return [W040, W131, W222, W220, W311]


def seidel_to_transverse_aberration(seidel, ref_index, slope):
    """ Convert Seidel coefficients to transverse ray aberrations """
    SI, SII, SIII, SIV, SV = seidel
    cnvrt = 1.0/(2.0*ref_index*slope)
    # TSA = transverse spherical aberration
    TSA = cnvrt*SI
    # TCO = tangential coma
    TCO = cnvrt*3.0*SII
    # TAS = tangential astigmatism
    TAS = cnvrt*(3.0*SIII + SIV)
    # SAS = sagittal astigmatism
    SAS = cnvrt*(SIII + SIV)
    # PTB = Petzval blur
    PTB = cnvrt*SIV
    # DST = distortion
    DST = cnvrt*SV
    return [TSA, TCO, TAS, SAS, PTB, DST]


def seidel_to_field_curv(seidel, ref_index, opt_inv):
    """ Convert Seidel coefficients to astigmatic and Petzval curvatures """
    SI, SII, SIII, SIV, SV = seidel
    cnvrt = ref_index/opt_inv**2
    # TCV = curvature of the tangential image surface
    TCV = cnvrt*(3.0*SIII + SIV)
    # SCV = curvature of the sagittal image surface
    SCV = cnvrt*(SIII + SIV)
    # PCV = curvature of the Petzval surface
    PCV = cnvrt*SIV
    return [TCV, SCV, PCV]


def transform_3rd_order(to_pkg, action):
    """ apply action to 3rd order package, to_pkg.
        to_pkg = to3, to3_total, to3_asp
        action = fct, (args)
        returns a 2d list of the transformed inputs
    """
    to3, to3_total, to3_asp = to_pkg
    fct, (args) = action
    has_asp = False if len(to3_asp) is 0 else True
    transformed_3rd_order = []
    for i, to3i in enumerate(to3[:-1]):
        ta = fct(to3i, *args)
        transformed_3rd_order.append(ta)
        if has_asp:
            try:
                to_asp = to3_asp[i+1]
            except KeyError:
                continue
            else:
                ta = fct(to_asp, *args)
                transformed_3rd_order.append(ta)
    transformed_3rd_order.append(fct(to3_total, *args))
    return transformed_3rd_order


def list_seidel_coefficients(to_pkg):
    to3, to3_total, to3_asp = to_pkg
    has_asp = False if len(to3_asp) is 0 else True
    print("\n                  SI          SII          SIII          SIV"
          "          SV ")
    for i, to3i in enumerate(to3[:-1]):
        print("{}         {:12.4g} {:12.4g} {:12.4g} {:12.4g} {:12.4g}"
              .format(i+1, to3i[0], to3i[1], to3i[2], to3i[3], to3i[4]))
        if has_asp:
            try:
                to_asp = to3_asp[i+1]
            except KeyError:
                continue
            else:
                print("asp       {:12.4g} {:12.4g} {:12.4g} {:12.4g} {:12.4g}"
                      .format(to_asp[0], to_asp[1], to_asp[2], to_asp[3],
                              to_asp[4]))
    print("Total     {:12.4g} {:12.4g} {:12.4g} {:12.4g} {:12.4g}"
          .format(*to3_total))


def list_3rd_order_trans_abr(to_pkg, n_last, u_last):
    to3, to3_total, to3_asp = to_pkg
    has_asp = False if len(to3_asp) is 0 else True
    print("\n                 TSA          TCO          TAS          SAS"
          "          PTB          DST")
    for i, to3i in enumerate(to3[:-1]):
        ta = seidel_to_transverse_aberration(to3i, n_last, u_last)
        print("{}         "
              "{:12.4g} {:12.4g} {:12.4g} {:12.4g} {:12.4g} {:12.4g}"
              .format(i+1, ta[0], ta[1], ta[2], ta[3], ta[4], ta[5]))
        if has_asp:
            try:
                to_asp = to3_asp[i+1]
            except KeyError:
                continue
            else:
                ta = seidel_to_transverse_aberration(to_asp, n_last, u_last)
                print("asp       {:12.4g} {:12.4g} {:12.4g} {:12.4g} {:12.4g}"
                      .format(ta[0], ta[1], ta[2], ta[3], ta[4]))
    tab_total = seidel_to_transverse_aberration(to3_total, n_last, u_last)
    print("Total     {:12.4g} {:12.4g} {:12.4g} {:12.4g} {:12.4g} {:12.4g}"
          .format(*tab_total))


def list_3rd_order_wvfrnt(to_pkg, central_wv):
    to3, to3_total, to3_asp = to_pkg
    has_asp = False if len(to3_asp) is 0 else True
    print("\nWavefront sums    W040         W131         W222         "
          "W220         W311")
    for i, to3i in enumerate(to3[:-1]):
        wv = seidel_to_wavefront(to3i, central_wv)
        print("{}         {:12.4g} {:12.4g} {:12.4g} {:12.4g} {:12.4g}"
              .format(i+1, wv[0], wv[1], wv[2], wv[3], wv[4]))
        if has_asp:
            try:
                to_asp = to3_asp[i+1]
            except KeyError:
                continue
            else:
                wv = seidel_to_wavefront(to_asp, central_wv)
                print("asp       {:12.4g} {:12.4g} {:12.4g} {:12.4g} {:12.4g}"
                      .format(wv[0], wv[1], wv[2], wv[3], wv[4]))
    wv_total = seidel_to_wavefront(to3_total, central_wv)
    print("Total     {:12.4g} {:12.4g} {:12.4g} {:12.4g} {:12.4g}".
          format(*wv_total))


def list_3rd_order_field_cv(to_pkg, n_last, opt_inv):
    to3, to3_total, to3_asp = to_pkg
    TCV, SCV, PCV = seidel_to_field_curv(to3_total, n_last, opt_inv)
    print("\nField Curvatures:   Tangential    Sagittal      Petzval")
    print("Totals           {:12.4g} {:12.4g} {:12.4g}".format(TCV, SCV, PCV))
