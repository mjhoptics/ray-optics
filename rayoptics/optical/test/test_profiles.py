#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 10:59:29 2017

@author: Mike
"""


import unittest
import rayoptics.optical.profile as p
from rayoptics.util.misc_math import normalize
import numpy as np
import numpy.testing as npt
from math import sqrt


class SphericalProfileTestCase(unittest.TestCase):
    def setUp(self):
        self.dir0 = np.array([0., 0., 1.])
        self.p0 = np.array([0., 0., -1.])
        self.p1 = np.array([0., 1., -1.])

    def test_normalize(self):
        v0 = normalize(np.array([0., 0., 0.]))
        npt.assert_array_equal(v0, np.array([0., 0., 0.]))

        v1 = normalize(np.array([1., 1., 1.]))
        sqrRt3 = sqrt(3)/3
        npt.assert_allclose(v1, np.array([sqrRt3, sqrRt3, sqrRt3]),
                            rtol=1e-14)

    def test_planar_sphere(self):
        s1 = p.Spherical(0.0)

        p0s1 = s1.intersect(self.p0, self.dir0)
        self.assertEqual(p0s1[0], 1.0)
        npt.assert_allclose(p0s1[1], np.array([0., 0., 0.]), rtol=1e-14)

        p1s1 = s1.intersect(self.p1, self.dir0)
        self.assertEqual(p1s1[0], 1.0)
        npt.assert_allclose(p1s1[1], np.array([0., 1., 0.]), rtol=1e-14)

    def test_convex_sphere(self):
        r2 = 10
        c2 = 1/r2
        s2 = p.Spherical(c2)
        sag2 = r2 - sqrt(r2*r2 - 1.0)

        p0s2 = s2.intersect(self.p0, self.dir0)
        self.assertEqual(p0s2[0], 1.0)
        npt.assert_allclose(p0s2[1], np.array([0., 0., 0.]), rtol=1e-14)

        p1s2 = s2.intersect(self.p1, self.dir0)
        npt.assert_allclose(p1s2[0], 1+sag2, rtol=1e-14)
        p1s2_truth = np.array([0., 1., sag2])
        npt.assert_allclose(p1s2[1], p1s2_truth, rtol=1e-14)

        dir_p1s2 = s2.normal(p1s2[1])
        dir_p1s2_truth = -normalize(p1s2_truth - np.array([0, 0, r2]))
        npt.assert_allclose(dir_p1s2, dir_p1s2_truth, rtol=1e-14)

    def test_concave_sphere(self):
        r3 = 10
        c3 = 1/r3
        s3 = p.Spherical(-c3)
        sag3 = r3 - sqrt(r3*r3 - 1.0)

        p0s3 = s3.intersect(self.p0, self.dir0)
        self.assertEqual(p0s3[0], 1.0)
        npt.assert_allclose(p0s3[1], np.array([0., 0., 0.]), rtol=1e-14)

        p1s3 = s3.intersect(self.p1, self.dir0)
        npt.assert_allclose(p1s3[0], 1-sag3, rtol=1e-14)
        p1s3_truth = np.array([0., 1., -sag3])
        npt.assert_allclose(p1s3[1], p1s3_truth, rtol=1e-14)

        dir_p1s3 = s3.normal(p1s3[1])
        dir_p1s3_truth = normalize(p1s3_truth - np.array([0, 0, -r3]))
        npt.assert_allclose(dir_p1s3, dir_p1s3_truth, rtol=1e-14)

    def test_dbgauss_s1(self):
        r1 = 56.20238
        c1 = 1/r1
        y0 = 25.0
        s1 = p.Spherical(c1)
        p1 = np.array([0., y0, 0.])
        sag1 = r1 - sqrt(r1*r1 - y0*y0)

        p0s1 = s1.intersect(self.p0, self.dir0)
        self.assertEqual(p0s1[0], 1.0)
        npt.assert_allclose(p0s1[1], np.array([0., 0., 0.]),
                            rtol=1e-14, atol=1e-14)

        p1s1 = s1.intersect(p1, self.dir0)
        npt.assert_allclose(p1s1[0], 5.866433424372758, rtol=1e-14)
        p1s1_truth = np.array([0., y0, sag1])
        npt.assert_allclose(p1s1[1], p1s1_truth, rtol=1e-14)

        dir_p1s1 = s1.normal(p1s1[1])
        dir_p1s1_truth = -normalize(p1s1_truth - np.array([0, 0, r1]))
        npt.assert_allclose(dir_p1s1, dir_p1s1_truth, rtol=1e-14)


if __name__ == '__main__':
    unittest.main(verbosity=3)
