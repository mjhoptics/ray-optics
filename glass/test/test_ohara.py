#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 09:59:50 2017

@author: Mike
"""

import unittest
import glass.ohara as o


class OharaTestCase(unittest.TestCase):
    catalog = o.OharaCatalog()

    def compare_indices(self, glass, tol=5e-6):
        nC = glass.rindex(656.27)
        nd = glass.rindex(587.56)
        ne = glass.rindex(546.07)
        nF = glass.rindex(486.13)
        ng = glass.rindex(435.84)
        nh = glass.rindex(404.66)
        # nI = glass.rindex(365.01)
        indxC = glass.glass_data()[self.catalog.data_index('nC')]
        indxd = glass.glass_data()[self.catalog.data_index('nd')]
        indxe = glass.glass_data()[self.catalog.data_index('ne')]
        indxF = glass.glass_data()[self.catalog.data_index('nF')]
        indxg = glass.glass_data()[self.catalog.data_index('ng')]
        indxh = glass.glass_data()[self.catalog.data_index('nh')]
        # indxI = glass.glass_data()[self.catalog.data_index('ni')]
        self.assertAlmostEqual(nC, indxC, delta=tol)
        self.assertAlmostEqual(nd, indxd, delta=tol)
        self.assertAlmostEqual(ne, indxe, delta=tol)
        self.assertAlmostEqual(nF, indxF, delta=tol)
        self.assertAlmostEqual(ng, indxg, delta=tol)
        self.assertAlmostEqual(nh, indxh, delta=tol)
        # self.assertAlmostEqual(nI, indxI, delta=tol)

    def test_ohara_catalog_glass_index(self):
        stim2 = self.catalog.glass_index('S-TIM 2')
        self.assertEqual(stim2, 48)
        sbsl7 = self.catalog.glass_index('S-BSL 7')
        self.assertEqual(sbsl7, 6)
        snph1 = self.catalog.glass_index('S-NPH 1')
        self.assertEqual(snph1, 132)

    def test_ohara_catalog_data_index(self):
        nd = self.catalog.data_index('nd')
        self.assertEqual(nd, 16)
        vd = self.catalog.data_index('νd')
        self.assertEqual(vd, 24)
        B1 = self.catalog.data_index('B1')
        self.assertEqual(B1, 63)
        glasscode = self.catalog.data_index('Code(d)')
        self.assertEqual(glasscode, 2)
        date = self.catalog.data_index(' D0')
        self.assertEqual(date, 154)

    def test_ohara_glass_stim2(self):
        stim2 = o.OharaGlass('S-TIM 2')
        self.assertIsNotNone(stim2.gindex)
        self.assertEqual(stim2.name(), 'S-TIM 2')
        self.compare_indices(stim2, tol=6e-6)

    def test_ohara_glass_sbsl7(self):
        sbsl7 = o.OharaGlass('S-BSL 7')
        self.assertIsNotNone(sbsl7.gindex)
        self.assertEqual(sbsl7.name(), 'S-BSL 7')
        self.compare_indices(sbsl7)

    def test_ohara_glass_snph1(self):
        snph1 = o.OharaGlass('S-NPH 1')
        self.assertIsNotNone(snph1.gindex)
        self.assertEqual(snph1.name(), 'S-NPH 1')
        self.compare_indices(snph1)


if __name__ == '__main__':
    unittest.main(verbosity=2)
