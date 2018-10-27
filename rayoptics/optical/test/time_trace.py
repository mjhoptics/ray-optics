#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 12:16:36 2018

@author: Mike
"""
import math
import numpy as np
from pathlib import Path
import timeit

# ray-optics
import rayoptics as ro
from rayoptics.util.misc_math import normalize


def setup(filename):
    root_pth = Path(ro.__file__).resolve().parent
    opm = ro.open_model(root_pth/filename)
    sm = opm.seq_model
    osp = sm.optical_spec
    fld, wvl, foc = osp.lookup_fld_wvl_focus(1)
    vig_pupil = fld.apply_vignetting([0.5, 0.5])
    fod = osp.parax_data.fod
    eprad = fod.enp_radius
    pt1 = np.array([eprad*vig_pupil[0], eprad*vig_pupil[1],
                    fod.obj_dist+fod.enp_dist])
    pt0 = osp.obj_coords(fld)
    dir0 = normalize(pt1 - pt0)
    return (sm, pt0, dir0, wvl)


def run_test(tst_name, trials, setup_stmt, file=None):
    output = '{:s} - {} trials: min time {:.2f}, ' \
             'max time {:.2f}, spread {:.2f}%'
    t = timeit.repeat('ray=rt.trace(*trace_args)', setup=setup_stmt,
                      number=trials, globals=globals())
    pcnt_sprd = 100*(max(t) - min(t))/min(t)
    print(output.format(tst_name, trials, min(t), max(t), pcnt_sprd),
          file=file)
    return [tst_name, trials, min(t), max(t), pcnt_sprd, t]


if __name__ == '__main__':
    trials10 = 10000
    trials15 = 15000
    trials20 = 20000
    trials25 = 25000
    trials40 = 40000
    trials50 = 50000
    setup_str = 'trace_args=setup("{:s}")'

    results = []

    root_pth = Path(ro.__file__).resolve().parent
    with open(root_pth/'optical/test/trace_results.txt', mode='w') as f:

        tst_name = 'singlet'
        trials = trials40
        setup_stmt = setup_str.format("codev/test/singlet.seq")
        results.append(run_test(tst_name, trials, setup_stmt, file=f))

        tst_name = 'landscape lens'
        trials = trials40
        setup_stmt = setup_str.format("codev/test/landscape_lens.seq")
        results.append(run_test(tst_name, trials, setup_stmt, file=f))

        tst_name = 'paraboloid'
        trials = trials40
        setup_stmt = setup_str.format("codev/test/paraboloid.seq")
        results.append(run_test(tst_name, trials, setup_stmt, file=f))

        tst_name = 'double gauss'
        trials = trials15
        setup_stmt = setup_str.format("codev/test/ag_dblgauss.seq")
        results.append(run_test(tst_name, trials, setup_stmt, file=f))

        tst_name = 'cell phone camera'
        trials = trials10
        setup_stmt = setup_str.format("optical/test/cell_phone_camera.roa")
        results.append(run_test(tst_name, trials, setup_stmt, file=f))

        tst_name = 'Sasian Triplet'
        trials = trials20
        setup_stmt = setup_str.format("../test/Sasian Triplet.roa")
        results.append(run_test(tst_name, trials, setup_stmt, file=f))

        tst_name = 'Cassegrain'
        trials = trials40
        setup_stmt = setup_str.format("../test/Cassegrain.roa")
        results.append(run_test(tst_name, trials, setup_stmt, file=f))

    with open(root_pth/'optical/test/trace_data.txt', mode='w') as f:
        print(results, file=f)
