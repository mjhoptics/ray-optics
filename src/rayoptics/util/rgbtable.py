#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 20:10:33 2018

@author: Mike
"""


import csv
from pathlib import Path


class RGBTable():
    data_list = []

    def __init__(self, filename='sunset_rgb.csv', data_range=[0., 100.]):
        self.data_range = data_range
        if len(self.data_list) == 0:
            path = Path(__file__).resolve().parent
            with open(path / filename, newline='') as f:
                reader = csv.reader(f)
                for row in reader:
                    self.data_list.append(list(map(int, row)))

    def get_color(self, value):
        fract = ((value - self.data_range[0]) /
                 (self.data_range[1] - self.data_range[0]))
        item_index = int(fract*len(self.data_list))
        return self.data_list[item_index]
