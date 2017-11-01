#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 09:59:50 2017

@author: Mike
"""

import unittest
import schott as s


class SchottTestCase(unittest.TestCase):
    def setUp(self):
        self.catalog = s.SchottCatalog()

    def test_schott_catalog_glass_index(self):
        f2 = self.catalog.glass_index('F2')
        self.assertEqual(f2, 4)
        nbk7 = self.catalog.glass_index('N-BK7')
        self.assertEqual(nbk7, 29)
        sf6ht = self.catalog.glass_index('SF6HT')
        self.assertEqual(sf6ht, 126)

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

    def test_schott_glass_nbk7(self):
        nbk7 = s.SchottGlass('N-BK7', self.catalog)
        self.assertIsNotNone(nbk7.gindex)
        self.assertEqual(nbk7.name(), 'N-BK7')


if __name__ == '__main__':
    unittest.main(verbosity=2)
