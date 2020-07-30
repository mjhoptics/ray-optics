#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 10:59:29 2017

@author: Mike
"""


import unittest
from pytest import approx
from rayoptics.elem.profiles import Spherical, EvenPolynomial
from rayoptics.util.misc_math import normalize
import numpy as np
import numpy.testing as npt
from math import sqrt


class SphericalProfileTestCase(unittest.TestCase):
    def setUp(self):
        self.dir0 = np.array([0., 0., 1.])
        self.p0 = np.array([0., 0., -1.])
        self.p1 = np.array([0., 1., -1.])
        self.eps = 1.0e-12
        self.z_dir = 1.0

    def test_normalize(self):
        v0 = normalize(np.array([0., 0., 0.]))
        npt.assert_array_equal(v0, np.array([0., 0., 0.]))

        v1 = normalize(np.array([1., 1., 1.]))
        sqrRt3 = sqrt(3)/3
        npt.assert_allclose(v1, np.array([sqrRt3, sqrRt3, sqrRt3]),
                            rtol=1e-14)

    def test_planar_sphere(self):
        s1 = Spherical(c=0.0)
        p_truth = 1.0, np.array([0., 0., 0.])

        # test for p0 ray, (0, 0, 0), (0, 0, 1)
        p0s1 = s1.intersect(self.p0, self.dir0, self.eps, self.z_dir)
        assert (p0s1[0], p0s1[1].all()) == (p_truth[0], p_truth[1].all())

        # test for p1 ray, (0, 1, 0), (0, 0, 1)
        p1s1 = s1.intersect(self.p1, self.dir0, self.eps, self.z_dir)
        assert ((p1s1[0], p1s1[1].all()) ==
                (p_truth[0], np.array([0., 1., 0.]).all()))

    def test_convex_sphere(self):
        r2 = 10
        c2 = 1/r2
        s2 = Spherical(c=c2)

        # test for p0 ray, (0, 0, 0), (0, 0, 1)
        p0_truth = 1.0, np.array([0., 0., 0.])

        p0s2_dir = s2.intersect(self.p0, self.dir0, self.eps, self.z_dir)
        p0s2_tnp = s2.intersect_tangent_plane(self.p0, self.dir0,
                                              self.eps, self.z_dir)
        assert ((p0s2_dir[0], p0s2_dir[1].all()) ==
                (p0s2_tnp[0], p0s2_tnp[1].all()))

        p0s2 = p0s2_dir
        assert ((p0s2[0], p0s2[1].all()) == (p0_truth[0], p0_truth[1].all()))

        # A spherical EvenPolynomial will use iteration to find the
        #  intersection. Use this case to further check results
        sa2 = EvenPolynomial(c=c2)
        p0sa2 = sa2.intersect(self.p0, self.dir0, self.eps, self.z_dir)
        assert ((p0sa2[0], p0sa2[1].all()) ==
                (p0_truth[0], p0_truth[1].all()))

        assert (p0s2_dir[0], p0s2_dir[1].all()) == (p0sa2[0], p0sa2[1].all())

        # test for p1 ray, (0, 1, 0), (0, 0, 1)
        sag2 = r2 - sqrt(r2*r2 - 1.0)
        p1_truth = 1+sag2, np.array([0., 1., sag2])

        p1s2_dir = s2.intersect(self.p1, self.dir0, self.eps, self.z_dir)
        p1s2_tnp = s2.intersect_tangent_plane(self.p1, self.dir0,
                                              self.eps, self.z_dir)
        assert ((p1s2_dir[0], p1s2_dir[1].all()) ==
                (p1s2_tnp[0], p1s2_tnp[1].all()))

        p1s2 = p1s2_dir
        assert ((p1s2[0], p1s2[1].all()) ==
                (p1_truth[0], p1_truth[1].all()))

        dir_p1s2 = s2.normal(p1s2[1])
        dir_p1s2_truth = -normalize(p1_truth[1] - np.array([0, 0, r2]))
        assert dir_p1s2.all() == dir_p1s2_truth.all()

    def test_concave_sphere(self):
        r3 = 10
        s3 = Spherical(r=-r3)

        # test for p0 ray, (0, 0, 0), (0, 0, 1)
        p0_truth = 1.0, np.array([0., 0., 0.])

        p0s3_dir = s3.intersect(self.p0, self.dir0, self.eps, self.z_dir)
        p0s3_tnp = s3.intersect_tangent_plane(self.p0, self.dir0,
                                              self.eps, self.z_dir)
        assert ((p0s3_dir[0], p0s3_dir[1].all()) ==
                (p0s3_tnp[0], p0s3_tnp[1].all()))

        p0s3 = p0s3_dir
        assert ((p0s3[0], p0s3[1].all()) == (p0_truth[0], p0_truth[1].all()))

        # test for p1 ray, (0, 1, 0), (0, 0, 1)
        sag3 = r3 - sqrt(r3*r3 - 1.0)
        p1_truth = 1-sag3, np.array([0., 1., -sag3])

        p1s3_tnp = s3.intersect_tangent_plane(self.p1, self.dir0,
                                              self.eps, self.z_dir)
        p1s3_dir = s3.intersect(self.p1, self.dir0, self.eps, self.z_dir)
        assert ((p1s3_dir[0], p1s3_dir[1].all()) ==
                (p1s3_tnp[0], p1s3_tnp[1].all()))

        p1s3 = p1s3_dir
        assert p1s3[0] == approx(p1_truth[0], rel=1e-14, abs=1e-14)
        assert p1s3[1].all() == p1_truth[1].all()

        dir_p1s3 = s3.normal(p1s3[1])
        dir_p1s3_truth = normalize(p1_truth[1] - np.array([0, 0, -r3]))
        assert dir_p1s3.all() == dir_p1s3_truth.all()

    def test_dbgauss_s1(self):
        r1 = 56.20238
        c1 = 1/r1
        y0 = 25.0
        s1 = Spherical(c=c1)
        p1 = np.array([0., y0, 0.])
        sag1 = r1 - sqrt(r1*r1 - y0*y0)

        # test for p0 ray, (0, 0, 0), (0, 0, 1)
        p0_truth = 1.0, np.array([0., 0., 0.])

        p0s1 = s1.intersect(self.p0, self.dir0, self.eps, self.z_dir)
        assert ((p0s1[0], p0s1[1].all()) ==
                (p0_truth[0], p0_truth[1].all()))
        self.assertEqual(p0s1[0], 1.0)
        npt.assert_allclose(p0s1[1], np.array([0., 0., 0.]),
                            rtol=1e-14, atol=1e-14)

        # test for p1 ray, (0, 1, 0), (0, 0, 1)
        p1_truth = 5.866433424372758, np.array([0., y0, sag1])
        p1s1 = s1.intersect(p1, self.dir0, self.eps, self.z_dir)
        assert p1s1[0] == approx(p1_truth[0], rel=1e-14, abs=1e-14)
        npt.assert_allclose(p1s1[1], p1_truth[1], rtol=1e-14)

        dir_p1s1 = s1.normal(p1s1[1])
        dir_p1s1_truth = -normalize(p1_truth[1] - np.array([0, 0, r1]))
        assert dir_p1s1.all() == dir_p1s1_truth.all()
        npt.assert_allclose(dir_p1s1, dir_p1s1_truth, rtol=1e-14)


if __name__ == '__main__':
    unittest.main(verbosity=3)
