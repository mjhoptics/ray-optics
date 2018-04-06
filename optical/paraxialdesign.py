#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" First order paraxial design space

Created on Sat Mar 31 21:14:42 2018

@author: Michael J. Hayford
"""

ax, pr, ln = range(3)
ht, slp, aoi = range(3)
pwr, tau, indx, rmd = range(4)


def build_lens(seq_model):
    sys = seq_model.seq_model_to_paraxial_lens()
    ax_ray, pr_ray, fod = seq_model.optical_spec.parax_data

    ax = []
    ax.append([r[ht] for r in ax_ray])
    ax.append([n*r[slp] for r, n in zip(ax_ray, sys[indx])])
    ax.append([n*r[aoi] for r, n in zip(ax_ray, sys[indx])])

    pr = []
    pr.append([r[ht] for r in pr_ray])
    pr.append([n*r[slp] for r, n in zip(pr_ray, sys[indx])])
    pr.append([n*r[aoi] for r, n in zip(pr_ray, sys[indx])])

    lens = []
    lens.append(ax)
    lens.append(pr)
    lens.append(sys)

    return lens


def ht_to_slope(lens, surf, opt_inv=1.0):
    """ This routine calculates all data dependent on the input
        height coordinates (y,ybar) at surface surf.
    """
    sys = lens[ln]
    ax_ray = lens[ax]
    pr_ray = lens[pr]

    nsm1 = len(sys[pwr]) - 1
    if surf == 0:
        surf += 1

    p = surf - 1
    c = surf

    sys[tau][p] = ((ax_ray[ht][p]*pr_ray[ht][c] - ax_ray[ht][c]*pr_ray[ht][p])
                   / opt_inv)
    ax_ray[slp][p] = (ax_ray[ht][c] - ax_ray[ht][p])/sys[tau][p]
    pr_ray[slp][p] = (pr_ray[ht][c] - pr_ray[ht][p])/sys[tau][p]

    if (surf > 1):
        p2 = surf - 2
        sys[pwr][p] = ((ax_ray[slp][p2]*pr_ray[slp][p] -
                        ax_ray[slp][p]*pr_ray[slp][p2])
                       / opt_inv)

    if (surf < nsm1):
        s = surf + 1

        sys[tau][c] = (ax_ray[ht][c]*pr_ray[ht][s] -
                       ax_ray[ht][s]*pr_ray[ht][c])/opt_inv
        ax_ray[slp][c] = (ax_ray[ht][s] - ax_ray[ht][c])/sys[tau][c]
        pr_ray[slp][c] = (pr_ray[ht][s] - pr_ray[ht][c])/sys[tau][c]
        sys[pwr][c] = (ax_ray[slp][p]*pr_ray[slp][c] -
                       ax_ray[slp][c]*pr_ray[slp][p])/opt_inv

        sys[pwr][s] = (ax_ray[slp][c]*pr_ray[slp][s] -
                       ax_ray[slp][s]*pr_ray[slp][c])/opt_inv

    else:
        ax_ray[slp][c] = ax_ray[slp][p]
        pr_ray[slp][c] = pr_ray[slp][p]
        sys[pwr][c] = 0
        sys[tau][c] = 0


def slope_to_ht(lens, surf, opt_inv=1.0):
    """ This routine calculates all data dependent on the input
        slope coordinates (nu,nubar) at surface surf.
    """
    sys = lens[ln]
    ax_ray = lens[ax]
    pr_ray = lens[pr]

    nsm1 = len(sys[pwr]) - 1
    if nsm1 == 0:
        p = 0
        c = 1
        ax_ray[ht][c] = ax_ray[slp][p]*sys[tau][p] + ax_ray[ht][p]
        pr_ray[ht][c] = pr_ray[slp][p]*sys[tau][p] + pr_ray[ht][p]

    else:
        if (surf == 0):
            surf += 1

        p = surf - 1
        c = surf

        sys[pwr][c] = (ax_ray[slp][p]*pr_ray[slp][c] -
                       ax_ray[slp][c]*pr_ray[slp][p])/opt_inv
        ax_ray[ht][c] = (ax_ray[slp][p] - ax_ray[slp][c])/sys[pwr][c]
        pr_ray[ht][c] = (pr_ray[slp][p] - pr_ray[slp][c])/sys[pwr][c]
        sys[tau][p] = (ax_ray[ht][p]*pr_ray[ht][c] -
                       ax_ray[ht][c]*pr_ray[ht][p])/opt_inv

        if (surf < nsm1):
            s = surf + 1
            ax_ray[ht][s] = ax_ray[slp][c]*sys[tau][c] + ax_ray[ht][c]
            pr_ray[ht][s] = pr_ray[slp][c]*sys[tau][c] + pr_ray[ht][c]


# ParaxTrace() - This routine performs a paraxial raytrace from object
#                (surface 0) to image.  The last operation is a
#                transfer to the image surface.
def paraxial_trace(lens):

    sys = lens[ln]
    ax_ray = lens[ax]
    pr_ray = lens[pr]

    nsm1 = len(sys[pwr]) - 1

    # Transfer from object
    c = 0
    s = 1
    ax_ray[ht][s] = ax_ray[ht][c] + sys[tau][c]*ax_ray[slp][c]
    pr_ray[ht][s] = pr_ray[ht][c] + sys[tau][c]*pr_ray[slp][c]

    for i in range(1, len(sys[pwr])):

        p = c
        c = s

        # Refraction
        ax_ray[slp][c] = ax_ray[slp][p] - ax_ray[ht][c]*sys[pwr][c]
        pr_ray[slp][c] = pr_ray[slp][p] - pr_ray[ht][c]*sys[pwr][c]

        # Transfer
        if (i < nsm1):
            s += 1
            ax_ray[ht][s] = ax_ray[ht][c] + sys[tau][c]*ax_ray[slp][c]
            pr_ray[ht][s] = pr_ray[ht][c] + sys[tau][c]*pr_ray[slp][c]
