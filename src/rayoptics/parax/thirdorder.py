#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
"""thirder order aberration calculation

.. Created on Fri Jul  6 07:24:40 2018

.. codeauthor: Michael J. Hayford
"""

import numpy as np
import pandas as pd

from rayoptics.optical.model_constants import ht, slp, aoi


def compute_third_order(opt_model):
    """ Compute Seidel aberration coefficents. """
    seq_model = opt_model.seq_model
    n_before = seq_model.central_rndx(0)
    parax_data = opt_model['analysis_results']['parax_data']
    ax_ray, pr_ray, fod = parax_data
    opt_inv = fod.opt_inv
    opt_inv_sqr = opt_inv*opt_inv

    third_order = {}

    # Transfer from object
    p = 0
    pd_index = ['S-I', 'S-II', 'S-III', 'S-IV', 'S-V']
    for c in range(1, len(ax_ray)-1):
        n_after = seq_model.central_rndx(c)
        n_after = n_after if seq_model.z_dir[c] > 0 else -n_after
        cv = seq_model.ifcs[c].profile_cv

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

        scoef = pd.Series([SIi, SIIi, SIIIi, SIVi, SVi], index=pd_index)
        col = str(c)
        third_order[col] = scoef

        # handle case of aspheric profile
        if hasattr(seq_model.ifcs[c], 'profile'):
            to_asp = aspheric_seidel_contribution(seq_model, parax_data, c,
                                                  n_before, n_after)
            if to_asp:
                ascoef = pd.Series(to_asp, index=pd_index)
                third_order[col+'.asp'] = ascoef

        p = c
        n_before = n_after

    third_order_df = pd.DataFrame(third_order, index=pd_index)
    third_order_df['sum'] = third_order_df.sum(axis='columns')
    return third_order_df.T


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


def aspheric_seidel_contribution(seq_model, parax_data, i, n_before, n_after):
    def delta_E(z, y, u, n):
        return -z/(n*y*(y + z*u))
    ax_ray, pr_ray, fod = parax_data

    if pr_ray[i][slp] == 0:
        e = pr_ray[i][ht]/ax_ray[i][ht]
    else:
        z = -pr_ray[i][ht]/pr_ray[i][slp]
        e = fod.opt_inv*delta_E(z, ax_ray[i][ht], ax_ray[i][slp], n_after)
    G = calc_4th_order_aspheric_term(seq_model.ifcs[i].profile)
    if G == 0.0:
        return None
    delta_n = n_after - n_before
    SI_star = 8.0*G*delta_n*ax_ray[i][ht]**4
    SII_star = SI_star*e
    SIII_star = SI_star*e**2
    SV_star = SI_star*e**3
    return [SI_star, SII_star, SIII_star, 0., SV_star]


def seidel_to_wavefront(seidel, central_wvl):
    """ Convert Seidel coefficients to wavefront aberrations """
    pd_index = ['W040', 'W131', 'W222', 'W220', 'W311']
    SI, SII, SIII, SIV, SV = seidel
    W040 = 0.125*SI/central_wvl
    W131 = 0.5*SII/central_wvl
    W222 = 0.5*SIII/central_wvl
    W220 = 0.25*(SIV + SIII)/central_wvl
    W311 = 0.5*SV/central_wvl
    return pd.Series([W040, W131, W222, W220, W311], index=pd_index)


def seidel_to_transverse_aberration(seidel, ref_index, slope):
    """ Convert Seidel coefficients to transverse ray aberrations """
    pd_index = ['TSA', 'TCO', 'TAS', 'SAS', 'PTB', 'DST']
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
    return pd.Series([TSA, TCO, TAS, SAS, PTB, DST], index=pd_index)


def seidel_to_field_curv(seidel, ref_index, opt_inv):
    """ Convert Seidel coefficients to astigmatic and Petzval curvatures """
    pd_index = ['TCV', 'SCV', 'PCV']
    SI, SII, SIII, SIV, SV = seidel
    cnvrt = ref_index/opt_inv**2
    # TCV = curvature of the tangential image surface
    TCV = cnvrt*(3.0*SIII + SIV)
    # SCV = curvature of the sagittal image surface
    SCV = cnvrt*(SIII + SIV)
    # PCV = curvature of the Petzval surface
    PCV = cnvrt*SIV
    return pd.Series([TCV, SCV, PCV], index=pd_index)
