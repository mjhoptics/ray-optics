#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 21 15:00:22 2018

@author: Mike
"""

from . import hoya as h
from . import ohara as o
from . import schott as s

Hoya, Ohara, Schott = range(3)
_cat_names = ["Hoya", "Ohara", "Schott"]


def create_glass(catalog, name):
    cat_name = catalog.capitalize()
    glass_name = name.upper()
    if cat_name == _cat_names[Hoya]:
        return h.HoyaGlass(glass_name)
    elif cat_name == _cat_names[Ohara]:
        return o.OharaGlass(glass_name)
    elif cat_name == _cat_names[Schott]:
        return s.SchottGlass(glass_name)
    else:
        return None


class GlassMapModel():

    def __init__(self):
        self.dataSetList = []
        self.dataSetList.append((h.HoyaCatalog(), _cat_names[Hoya]))
        self.dataSetList.append((o.OharaCatalog(), _cat_names[Ohara]))
        self.dataSetList.append((s.SchottCatalog(), _cat_names[Schott]))

    def get_data_at(self, i):
        return self.dataSetList[i][0].glass_map_data()

    def get_data_set_label_at(self, i):
        return self.dataSetList[i][1]
