#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 21:09:42 2017

@author: Mike
"""

# from __future__ import print_function
import xlrd
from xlrd.sheet import ctype_text


# cd /Users/Mike/Developer/PyProjects/optics/optical

# Open the workbook
dir_name = '/Users/Mike/Developer/PyProjects/optics/optical/'
fname = 'schott-optical-glass-overview-excel-table-english-06032017.xls'
xl_workbook = xlrd.open_workbook(dir_name+fname)
# List sheet names, and pull a sheet by name
#
sheet_names = xl_workbook.sheet_names()
print('Sheet Names', sheet_names)

xl_sheet = xl_workbook.sheet_by_name(sheet_names[0])
print('Sheet name: %s' % xl_sheet.name)

# Pull the first row by index
#  (rows/columns are also zero-indexed)
#
row = xl_sheet.row(0)  # 1st row

# Print 1st row values and types
#

print('(Column #) type:value')
for idx, cell_obj in enumerate(row):
    cell_type_str = ctype_text.get(cell_obj.ctype, 'unknown type')
    print('(%s) %s %s' % (idx, cell_type_str, cell_obj.value))

# Print all values, iterating through rows and columns
#
num_cols = xl_sheet.ncols   # Number of columns

row = xl_sheet.row(1)
for idx, cell_obj in enumerate(row):
    cell_type_str = ctype_text.get(cell_obj.ctype, 'unknown type')
    print('(%s) %s %s' % (idx, cell_type_str, cell_obj.value))

row = xl_sheet.row(2)
for idx, cell_obj in enumerate(row):
    cell_type_str = ctype_text.get(cell_obj.ctype, 'unknown type')
    print('(%s) %s %s' % (idx, cell_type_str, cell_obj.value))

row = xl_sheet.row(3)
for idx, cell_obj in enumerate(row):
    cell_type_str = ctype_text.get(cell_obj.ctype, 'unknown type')
    print('(%s) %s %s' % (idx, cell_type_str, cell_obj.value))

row = xl_sheet.row(4)
for idx, cell_obj in enumerate(row):
    cell_type_str = ctype_text.get(cell_obj.ctype, 'unknown type')
    print('(%s) %s %s' % (idx, cell_type_str, cell_obj.value))

gnames = xl_sheet.col_values(0, 4)
