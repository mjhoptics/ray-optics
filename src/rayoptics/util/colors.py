#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2020 Michael J. Hayford
""" manage light and dark interface color schemes

.. Created on Wed Jan 29 20:00:41 2020

.. codeauthor: Michael J. Hayford
"""


solarize = {
    'base03': '#002b36',
    'base02': '#073642',
    'base01': '#586e75',
    'base00': '#657b83',
    'base0': '#839496',
    'base1': '#93a1a1',
    'base2': '#eee8d5',
    'base3': '#fdf6e3',
    'violet': '#6c71c4',
    'blue': '#268bd2',
    'cyan': '#2aa198',
    'green': '#859900',
    'yellow': '#b58900',
    'orange': '#cb4b16',
    'red': '#dc322f',
    'magenta': '#d33682',
    }


def accent_colors(is_dark=True):
    accent = {
        'violet': solarize['violet'],
        'blue': solarize['blue'],
        'cyan': solarize['cyan'],
        'green': solarize['green'],
        'yellow': solarize['yellow'],
        'orange': solarize['orange'],
        'red': solarize['red'],
        'magenta': solarize['magenta'],
        }
    return accent


def foreground_background(is_dark=True):
    if is_dark:
        rgb = {
            'background': solarize['base03'],
            'background1': solarize['base02'],
            'background2': solarize['base01'],
            'foreground': solarize['base0'],
            'foreground1': solarize['base1'],
            'foreground2': solarize['base01'],
            'hilite': solarize['base1']
            }
    else:
        rgb = {
            'background': solarize['base3'],
            'background1': solarize['base2'],
            'background2': solarize['base1'],
            'foreground': solarize['base00'],
            'foreground1': solarize['base01'],
            'foreground2': solarize['base1'],
            'hilite': solarize['base01']
            }
    return rgb
