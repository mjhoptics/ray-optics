#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
""" dict2D constructs a M x N dict using row_keys and col_keys as keys

.. Created on Wed Aug 14 21:08:56 2019

.. codeauthor: Michael J. Hayford
"""


def dict2D(row_keys, col_keys):
    """ returns a 2D dictionary with M row_keys and N col_keys """
    dict_2d = {}
    for r in row_keys:
        row = {}
        dict_2d[r] = row
        for c in col_keys:
            row[c] = {}
    return dict_2d


def row(dict_2d, row_key):
    """ returns a dict of the contents of row **row_key** of dict_2d """
    return dict_2d[row_key]


def col(dict_2d, col_key):
    """ returns a dict of the contents of column **col_key** of dict_2d """
    column = {}
    row_keys = dict_2d.keys()
    for r in row_keys:
        column[r] = dict_2d[r][col_key]
    return column
#    dict([(r, dict_2d[r][col_key]) for r in dict_2d.keys()])


def num_items_by_type(dict_2d, row_keys, col_keys):
    """ return a dict of the number of items in each row/col of dict_2d """
    len_items = dict([(rc, 0) for rc in row_keys+col_keys])
    for r in row_keys:
        for c in dict_2d[r]:
            num_items = len(dict_2d[r][c])
            len_items[r] += num_items
            len_items[c] += num_items
    return len_items


def num_items_by_cell(dict_2d, row_keys, col_keys):
    """ return a list of the number of items in each cell of dict_2d """
    num_items = []
    for r in row_keys:
        for c in dict_2d[r]:
            ni = len(dict_2d[r][c])
            num_items.append((r, c, ni))
    return num_items
