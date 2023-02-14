#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 12:16:36 2018

@author: Mike
"""
import math
import numpy as np
from pathlib import Path
import csv
import timeit

# ray-optics
import rayoptics as ro
from rayoptics.gui.appcmds import open_model
import rayoptics.raytr.raytrace as rt
from rayoptics.util.misc_math import normalize


def setup(filename):
    root_pth = Path(ro.__file__).resolve().parent
    opm = open_model(root_pth/filename)
    sm = opm.seq_model
    osp = opm.optical_spec
    fld, wvl, foc = osp.lookup_fld_wvl_focus(1)
    vig_pupil = fld.apply_vignetting([0.5, 0.5])
    fod = opm['analysis_results']['parax_data'].fod
    eprad = fod.enp_radius
    pt1 = np.array([eprad*vig_pupil[0], eprad*vig_pupil[1],
                    fod.obj_dist+fod.enp_dist])
    pt0, d0 = osp.obj_coords(fld)
    dir0 = normalize(pt1 - pt0)
    return (sm, pt0, dir0, wvl)


def run_test(tst_name, setup_stmt, trials, repeat=5, file=None):
    output = '{:s} - {} trials: min time {:.2f}, ' \
             'max time {:.2f}, spread {:.2f}%'
    t = timeit.repeat('ray=rt.trace(*trace_args)', setup=setup_stmt,
                      number=trials, repeat=repeat, globals=globals())
    pcnt_sprd = 100*(max(t) - min(t))/min(t)
    print(output.format(tst_name, trials, min(t), max(t), pcnt_sprd),
          file=file)
    return [tst_name, trials, min(t), max(t), pcnt_sprd, t]


if __name__ == '__main__':
    trials10 = 10000
    trials15 = 15000
    trials20 = 20000
    trials25 = 25000
    trials30 = 30000
    trials40 = 40000
    trials50 = 50000
    trials80 = 80000
    trials100 = 100000
    setup_str = 'trace_args=setup("{:s}")'

    results = []
    rpt = 5

    root_pth = Path(ro.__file__).resolve().parent
    with open(root_pth/'raytr/tests/trace_results.txt', mode='w') as f:

        tst_name = 'singlet'
        trials = trials100
        setup_stmt = setup_str.format("codev/tests/singlet.seq")
        results.append(run_test(tst_name, setup_stmt, trials, file=f,
                                repeat=rpt))

        tst_name = 'landscape lens'
        trials = trials80
        setup_stmt = setup_str.format("codev/tests/landscape_lens.seq")
        results.append(run_test(tst_name, setup_stmt, trials, file=f,
                                repeat=rpt))

        tst_name = 'Sasian triplet'
        trials = trials50
        setup_stmt = setup_str.format("models/Sasian Triplet.roa")
        results.append(run_test(tst_name, setup_stmt, trials, file=f,
                                repeat=rpt))

        tst_name = 'double gauss'
        trials = trials30
        setup_stmt = setup_str.format("codev/tests/ag_dblgauss.seq")
        results.append(run_test(tst_name, setup_stmt, trials, file=f,
                                repeat=rpt))

        tst_name = '2 spherical mirrors (spheres)'
        trials = trials100
        setup_stmt = setup_str.format("models/TwoSphericalMirror.roa")
        results.append(run_test(tst_name, setup_stmt, trials, file=f,
                                repeat=rpt))

        tst_name = '2 spherical mirrors (conics)'
        trials = trials100
        setup_stmt = setup_str.format("models/TwoMirror.roa")
        results.append(run_test(tst_name, setup_stmt, trials, file=f,
                                repeat=rpt))

        tst_name = 'paraboloid'
        trials = trials100
        setup_stmt = setup_str.format("codev/tests/paraboloid.seq")
        results.append(run_test(tst_name, setup_stmt, trials, file=f,
                                repeat=rpt))

        tst_name = 'Cassegrain'
        trials = trials100
        setup_stmt = setup_str.format("models/Cassegrain.roa")
        results.append(run_test(tst_name, setup_stmt, trials, file=f,
                                repeat=rpt))

        tst_name = 'Ritchey-Chretien'
        trials = trials100
        setup_stmt = setup_str.format("models/Ritchey_Chretien.roa")
        results.append(run_test(tst_name, setup_stmt, trials, file=f,
                                repeat=rpt))

        tst_name = 'cell phone camera'
        trials = trials10
        setup_stmt = setup_str.format("optical/tests/cell_phone_camera.roa")
        results.append(run_test(tst_name, setup_stmt, trials, file=f,
                                repeat=rpt))

    with open(root_pth/'raytr/tests/trace_data.csv', mode='w') as f:
        w = csv.writer(f, delimiter=',', quotechar='"',
                       quoting=csv.QUOTE_MINIMAL)
        for tst in results:
            w.writerow(tst[:5])

    with open(root_pth/'raytr/tests/trace_data.txt', mode='w') as f:
        print(results, file=f)
