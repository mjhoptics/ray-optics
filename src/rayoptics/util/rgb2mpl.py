#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" convert RGB data to matplotlib format

.. Created on Mon Jul 15 12:24:33 2019

.. codeauthor: Michael J. Hayford
"""


def rgb2mpl(rgb):
    """ convert 8 bit RGB data to 0 to 1 range for mpl """
    if len(rgb) == 3:
        return [rgb[0]/255., rgb[1]/255., rgb[2]/255.]
    elif len(rgb) == 4:
        return [rgb[0]/255., rgb[1]/255., rgb[2]/255., rgb[3]/255.]


backgrnd_color = rgb2mpl([237, 243, 254])  # light blue
