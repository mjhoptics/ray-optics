#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
"""
Created on Fri Jul  6 07:24:40 2018

@author: Michael J. Hayford
"""

from optical.model_constants import ax, pr, lns, inv
from optical.model_constants import ht, slp, aoi
from optical.model_constants import pwr, tau, indx, rmd


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

        SI += SIi
        SII += SIIi
        SIII += SIIIi
        SIV += SIVi
        SV += SVi

        third_order.append([SIi, SIIi, SIIIi, SIVi, SVi])

        p = c
        n_before = n_after

    return (third_order, [SI, SII, SIII, SIV, SV])


def seidel_to_wavefront(seidel, central_wvl):
    """ Convert Seidel coefficients to wavefront aberrations """
    SI, SII, SIII, SIV, SV = seidel
    W040 = 0.125*SI/central_wvl
    W131 = 0.5*SII/central_wvl
    W222 = 0.5*SIII/central_wvl
    W220 = 0.25*(SIV + SIII)/central_wvl
    W311 = 0.5*SV/central_wvl
    return [W040, W131, W222, W220, W311]


def seidel_to_transverse_aberration(seidel, indx, slope):
    """ Convert Seidel coefficients to transverse ray aberrations """
    SI, SII, SIII, SIV, SV = seidel
    cnvrt = 1.0/2.0*indx*slope
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


def seidel_to_field_curv(seidel, indx, opt_inv):
    """ Convert Seidel coefficients to astigmatic and Petzval curvatures """
    SI, SII, SIII, SIV, SV = seidel
    cnvrt = indx/opt_inv*opt_inv
    # TCV = curvature of the tangential image surface
    TCV = cnvrt*(3.0*SIII + SIV)
    # SCV = curvature of the sagittal image surface
    SCV = cnvrt*(SIII + SIV)
    # PCV = curvature of the Petzval surface
    PCV = cnvrt*SIV
    return [TCV, SCV, PCV]


def list_seidel_coefficients(to_pckg):
    to3, to3_total = to_pckg
    print("\n                  SI          SII          SIII          SIV"
          "          SV ")
    for i, to3i in enumerate(to3):
        print("{}         {:12.4g} {:12.4g} {:12.4g} {:12.4g} {:12.4g}"
              .format(i+1, to3i[0], to3i[1], to3i[2], to3i[3], to3i[4]))
    print("Total     {:12.4g} {:12.4g} {:12.4g} {:12.4g} {:12.4g}"
          .format(*to3_total))


def list_3rd_order_trans_abr(to_pckg, n_last, u_last):
    to3, to3_total = to_pckg
    print("\n                 TSA          TCO          TAS          SAS"
          "          PTB          DST")
    for i, to3i in enumerate(to3):
        ta = seidel_to_transverse_aberration(to3i, n_last, u_last)
        print("{}         "
              "{:12.4g} {:12.4g} {:12.4g} {:12.4g} {:12.4g} {:12.4g}"
              .format(i+1, ta[0], ta[1], ta[2], ta[3], ta[4], ta[5]))
    tab_total = seidel_to_transverse_aberration(to3_total, n_last, u_last)
    print("Total     {:12.4g} {:12.4g} {:12.4g} {:12.4g} {:12.4g} {:12.4g}"
          .format(*tab_total))


def list_3rd_order_wvfrnt(to_pckg, central_wv):
    to3, to3_total = to_pckg
    print("\nWavefront sums    W040         W131         W222         "
          "W220         W311")
    for i, to3i in enumerate(to3):
        wv = seidel_to_wavefront(to3i, central_wv)
        print("{}         {:12.4g} {:12.4g} {:12.4g} {:12.4g} {:12.4g}"
              .format(i+1, wv[0], wv[1], wv[2], wv[3], wv[4]))
    wv_total = seidel_to_wavefront(to3_total, central_wv)
    print("Total     {:12.4g} {:12.4g} {:12.4g} {:12.4g} {:12.4g}".
          format(*wv_total))


def list_3rd_order_field_cv(to_pckg, n_last, opt_inv):
    to3, to3_total = to_pckg
    TCV, SCV, PCV = seidel_to_field_curv(to3_total, n_last, opt_inv)
    print("\nField Curvatures:   Tangential    Sagittal      Petzval")
    print("Totals           {:12.4g} {:12.4g} {:12.4g}".format(TCV, SCV, PCV))
