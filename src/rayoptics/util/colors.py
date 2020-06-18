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


# Replace base2 and base3 with blue tinted hues
solarize_blue = {
    'base03': '#002b36',
    'base02': '#073642',
    'base01': '#586e75',
    'base00': '#657b83',
    'base0': '#839496',
    'base1': '#93a1a1',
    'base2': '#e1e7f2',
    'base3': '#edf3fe',
    # 'base2': '#dae5e7',
    # 'base3': '#e9f3f6',
    # 'base2': '#d4e7ed',
    # 'base3': '#e9f4f6',
    # 'base2': '#d4eaed',
    # 'base3': '#e2f9fd',
    'violet': '#6c71c4',
    'blue': '#268bd2',
    'cyan': '#2aa198',
    'green': '#859900',
    'yellow': '#b58900',
    'orange': '#cb4b16',
    'red': '#dc322f',
    'magenta': '#d33682',
    }


solarize_dict = solarize_blue


def accent_colors(is_dark=True):
    accent = {
        'violet': solarize_dict['violet'],
        'blue': solarize_dict['blue'],
        'cyan': solarize_dict['cyan'],
        'green': solarize_dict['green'],
        'yellow': solarize_dict['yellow'],
        'orange': solarize_dict['orange'],
        'red': solarize_dict['red'],
        'magenta': solarize_dict['magenta'],
        }
    return accent


def foreground_background(is_dark=True):
    if is_dark:
        rgb = {
            'background': solarize_dict['base03'],
            'background1': solarize_dict['base02'],
            'background2': solarize_dict['base01'],
            'foreground': solarize_dict['base0'],
            'foreground1': solarize_dict['base1'],
            'foreground2': solarize_dict['base01'],
            'hilite': solarize_dict['base1']
            }
    else:
        rgb = {
            'background': solarize_dict['base3'],
            'background1': solarize_dict['base2'],
            'background2': solarize_dict['base1'],
            'foreground': solarize_dict['base00'],
            'foreground1': solarize_dict['base01'],
            'foreground2': solarize_dict['base1'],
            'hilite': solarize_dict['base01']
            }
    return rgb
