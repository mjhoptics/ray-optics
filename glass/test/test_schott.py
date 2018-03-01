#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 09:59:50 2017

@author: Mike
"""

import unittest
import glass.schott as s


class SchottTestCase(unittest.TestCase):
    catalog = s.SchottCatalog()

    def compare_indices(self, glass, tol=5e-6):
        nC = glass.rindex(656.27)
        nd = glass.rindex(587.56)
        ne = glass.rindex(546.07)
        nF = glass.rindex(486.13)
        ng = glass.rindex(435.84)
        nh = glass.rindex(404.66)
        nI = glass.rindex(365.01)
        indxC = glass.glass_data()[self.catalog.data_index('  nC')]
        indxd = glass.glass_data()[self.catalog.data_index('  nd')]
        indxe = glass.glass_data()[self.catalog.data_index('  ne')]
        indxF = glass.glass_data()[self.catalog.data_index('  nF')]
        indxg = glass.glass_data()[self.catalog.data_index('  ng')]
        indxh = glass.glass_data()[self.catalog.data_index('  nh')]
        indxI = glass.glass_data()[self.catalog.data_index('  ni')]
        self.assertAlmostEqual(nC, indxC, delta=tol)
        self.assertAlmostEqual(nd, indxd, delta=tol)
        self.assertAlmostEqual(ne, indxe, delta=tol)
        self.assertAlmostEqual(nF, indxF, delta=tol)
        self.assertAlmostEqual(ng, indxg, delta=tol)
        self.assertAlmostEqual(nh, indxh, delta=tol)
        self.assertAlmostEqual(nI, indxI, delta=tol)

    def test_schott_catalog_glass_index(self):
        f2 = self.catalog.glass_index('F2')
        self.assertEqual(f2, 0)
        nbk7 = self.catalog.glass_index('N-BK7')
        self.assertEqual(nbk7, 25)
        sf6ht = self.catalog.glass_index('SF6HT')
        self.assertEqual(sf6ht, 122)

    def test_schott_catalog_data_index(self):
        nd = self.catalog.data_index('nd')
        self.assertEqual(nd, 1)
        vd = self.catalog.data_index('vd')
        self.assertEqual(vd, 3)
        B1 = self.catalog.data_index('B1')
        self.assertEqual(B1, 6)
        glasscode = self.catalog.data_index('Glascode')
        self.assertEqual(glasscode, 158)
        date = self.catalog.data_index('Date')
        self.assertEqual(date, 160)

    def test_schott_glass_f2(self):
        f2 = s.SchottGlass('F2')
        self.assertIsNotNone(f2.gindex)
        self.assertEqual(f2.name(), 'F2')
        self.compare_indices(f2)

    def test_schott_glass_nbk7(self):
        nbk7 = s.SchottGlass('N-BK7')
        self.assertIsNotNone(nbk7.gindex)
        self.assertEqual(nbk7.name(), 'N-BK7')
        self.compare_indices(nbk7)

    def test_schott_glass_sf6ht(self):
        sf6ht = s.SchottGlass('SF6HT')
        self.assertIsNotNone(sf6ht.gindex)
        self.assertEqual(sf6ht.name(), 'SF6HT')
        self.compare_indices(sf6ht)


if __name__ == '__main__':
    unittest.main(verbosity=2)
