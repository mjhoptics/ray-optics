#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2020 Michael J. Hayford
""" manage light and dark interface color schemes

.. Created on Mon Feb  3 22:03:16 2020

.. codeauthor: Michael J. Hayford
"""
from pathlib import Path
from shutil import copy2

import matplotlib
from matplotlib import style
from matplotlib.figure import Figure

from rayoptics.util import colors


def copy_styles():
    """Copy rayoptics mpl styles to user's mpl_config dir."""
    pth = Path(__file__).resolve().parent
    styles_dir = Path(pth / 'styles')
    mpl_configdir = Path(matplotlib.get_configdir()) / 'stylelib'
    mpl_configdir.mkdir(exist_ok=True)
    for mpl_style in styles_dir.glob('*.mplstyle'):
        copy2(mpl_style, mpl_configdir)


def apply_style(is_dark):
    """Assign a light or dark style to mpl plots."""
    pth = Path(__file__).resolve().parent
    styles_dir = Path(pth / 'styles')
    if is_dark:
        style_path = styles_dir / 'Solarize_Dark.mplstyle'
        style.use(str(style_path))
    else:
        style_path = styles_dir / 'Solarize_Light_Blue.mplstyle'
        style.use(str(style_path))


class StyledFigure(Figure):
    """Provide a standard implementation for mpl styles."""

    def __init__(self, **kwargs):
        is_dark = kwargs.pop('is_dark', False)
        self._rgb = {**colors.accent_colors(is_dark),
                     **colors.foreground_background(is_dark)}

        super().__init__(**kwargs)

        self.sync_light_or_dark(is_dark, do_refresh=False)

    def sync_light_or_dark(self, is_dark, do_refresh=True):
        self._rgb = {**colors.accent_colors(is_dark),
                     **colors.foreground_background(is_dark)}
        apply_style(is_dark)
        self.set_facecolor(self._rgb['background'])
        if hasattr(self, 'ax'):
            axes = self.ax
            axes.set_facecolor(self._rgb['background1'])
            axes.tick_params(
                colors=self._rgb['foreground'],
                grid_color=self._rgb['background']
                )
        if do_refresh:
            self.refresh()
