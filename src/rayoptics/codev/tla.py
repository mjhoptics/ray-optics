#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Support for CODE V TLAs

.. codeauthor: Michael J. Hayford
"""
import csv
from pathlib import Path


class MapTLA:
    """ Create and maintain a dictionary of CODE V 3 letter commands """
    _d = {}

    def __init__(self):
        TLA, CmdFct, IndxQuals, DataType, Quals = range(5)
        if len(MapTLA._d) == 0:
            path = Path(__file__).resolve().parent
            with open(path / 'tla_mapping.csv') as f:
                reader = csv.reader(f)
                for row in reader:
                    if row[TLA] != '':
                        if row[Quals] != '':
                            row[Quals] = row[Quals].split(',')
                        MapTLA._d[row[TLA]] = row[CmdFct:]

    def find(self, tla):
        try:
            return MapTLA._d[tla]
        except KeyError:
            return None
