#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 10:11:18 2017

@author: Mike
"""


from math import sqrt
import xlrd


class SchottCatalog:
    xl_data = None
    data_header = 3
    data_start = 4
    num_glasses = 123
    name_col_offset = 0
    coef_col_offset = 6
    index_col_offset = 117

    def __init__(self):
        # Open the workbook
        dname = '/Users/Mike/Developer/PyProjects/ray-optics/glass/'
        fname = 'schott-optical-glass-overview-excel-table-english-06032017'\
                '.xls'
        xl_workbook = xlrd.open_workbook(dname+fname)
        self.xl_data = xl_workbook.sheet_by_index(0)
        # self.name_col_offset = self.xl_data.row_values(1, 0).index('Glass')
        self.data_header = (self.xl_data.col_values(self.name_col_offset, 0)
                            .index('Glass'))
        self.data_start = self.data_header+1
        gnames = self.xl_data.col_values(self.name_col_offset, self.data_start)
        while gnames and len(gnames[-1]) is 0:
            gnames.pop()
        self.num_glasses = len(gnames)
        colnames = self.xl_data.row_values(self.data_header, 0)
        self.coef_col_offset = colnames.index('B1')
        self.index_col_offset = colnames.index('  n2325.4')

    def glass_index(self, gname):
        gnames = self.xl_data.col_values(self.name_col_offset, self.data_start)
        if gname in gnames:
            gindex = gnames.index(gname)
        else:
            gindex = None

        return gindex

    def data_index(self, dname):
        if dname in self.xl_data.row_values(self.data_header, 0):
            dindex = self.xl_data.row_values(self.data_header, 0).index(dname)
        else:
            dindex = None

        return dindex

    def glass_data(self, row):
        return self.xl_data.row_values(self.data_start+row, 0)

    def catalog_data(self, col):
        return self.xl_data.col_values(col, self.data_start,
                                       self.data_start+self.num_glasses)

    def glass_coefs(self, gindex):
        return (self.xl_data.row_values(self.data_start+gindex,
                                        self.coef_col_offset,
                                        self.coef_col_offset+6))

    def glass_map_data(self, wvl='d'):
        if wvl == 'd':
            nd = self.catalog_data(self.data_index('nd'))
            vd = self.catalog_data(self.data_index('vd'))
            names = self.catalog_data(self.name_col_offset)
            return vd, nd, names
        elif wvl == 'e':
            ne = self.catalog_data(self.data_index('ne'))
            ve = self.catalog_data(self.data_index('ve'))
            names = self.catalog_data(self.name_col_offset)
            return ve, ne, names
        else:
            return None


class SchottGlass:
    catalog = SchottCatalog()

    def __init__(self, gname):
        self.gindex = self.catalog.glass_index(gname)

    def __repr__(self):
        return 'Schott ' + self.name() + ': ' + self.glass_code()

    def glass_code(self):
        nd = self.glass_item('nd')
        vd = self.glass_item('vd')
        return str(1000*round((nd - 1), 3) + round(vd/100, 3))

    def glass_data(self):
        return self.catalog.glass_data(self.gindex)

    def name(self):
        return self.catalog.xl_data.cell_value(
                self.catalog.data_start+self.gindex, 0)

    def glass_item(self, dname):
        dindex = self.catalog.data_index(dname)
        if dindex is None:
            return None
        else:
            return self.glass_data()[dindex]

    def rindex(self, wv_nm):
        wv = 0.001*wv_nm
        wv2 = wv*wv
        coefs = self.catalog.glass_coefs(self.gindex)
        n2 = 1. + coefs[0]*wv2/(wv2 - coefs[3])
        n2 = n2 + coefs[1]*wv2/(wv2 - coefs[4])
        n2 = n2 + coefs[2]*wv2/(wv2 - coefs[5])
        return sqrt(n2)
