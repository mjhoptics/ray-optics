#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2020 Michael J. Hayford
"""Dashboards constructed with ipywidgets.

.. Created on Tue Apr 28 21:15:20 2020

.. codeauthor: Michael J. Hayford
"""

import numpy as np
import ipywidgets as widgets


def create_focus_dashboard(figs, ray_data_items, 
                           foc: float, dfoc_rng: float, shift_rng: float,
                           image_delta, continuous_update: bool=True):
    """Create 3 ipywidgets to modify defocus and image shift interactively.
    
    This function creates 3 ipywidgets (FloatSliders) controlling defoucus, 
    x_shift and y_shift. When the user moves one of the sliders, the callback
    function applies the new slider value to the ray_data_items and causes
    them to 'update'. The ray_data items in the :mod:`~raytr.analyses` package
    all provide a lightweight 'update' implementation that doesn't require
    retracing the model.
    
    Arguments:
        figs: list of Figures to be cleared and replotted
        ray_data_items: list of ray data items to be updated to the new values
        foc: the current defocus
        dfoc_range: the maximum range for the defocus slider
        shift_rng: the maximum range for the image shift sliders
        image_delta: the current image shift, [delta_x: float, delta_y: float]
        continuous_update: if True, model is updated for any slider movement
    """
    if image_delta is None:
        image_delta = [0., 0.]
    defocus = widgets.FloatSlider(min=-dfoc_rng, max=+dfoc_rng,
                                  step=.01*dfoc_rng,
                                  description='defocus', value=foc,
                                  readout_format='.4f',
                                  continuous_update=continuous_update)
    x_shift = widgets.FloatSlider(min=-shift_rng, max=+shift_rng,
                                  step=.01*shift_rng,
                                  description='x shift', value=image_delta[0],
                                  readout_format='.4f',
                                  continuous_update=continuous_update)
    y_shift = widgets.FloatSlider(min=-shift_rng, max=+shift_rng,
                                  step=.01*shift_rng,
                                  description='y shift', value=image_delta[1],
                                  readout_format='.4f',
                                  continuous_update=continuous_update)

    def slider_update(change):
        dfoc_val = defocus.value
        dx = x_shift.value
        dy = y_shift.value

        # apply changes to fans and grids
        for ray_data in ray_data_items:
            ray_data.foc = dfoc_val
            ray_data.image_delta = np.array([dx, dy])
            ray_data.update_data(build='update')

        # update and plot results
        for fig in figs:
            fig.clf()
            fig.plot()

    defocus.observe(slider_update, names='value')
    x_shift.observe(slider_update, names='value')
    y_shift.observe(slider_update, names='value')

    return defocus, x_shift, y_shift


def create_mirror_tilt_dashboard(mirror, app_mgr, figs, ray_data_items,
                                 foc, tilt_rng, shift_rng, oa_ray=None,
                                 continuous_update=True):
    tilt0 = mirror.decenter.euler[0]
    alpha_tilt = widgets.FloatSlider(min=tilt0-tilt_rng, max=tilt0+tilt_rng,
                                     step=.05*tilt_rng,
                                     description='alpha tilt', value=tilt0,
                                     readout_format='.4f',
                                     continuous_update=continuous_update)
    x_shift = widgets.FloatSlider(min=-shift_rng, max=+shift_rng,
                                  step=.05*shift_rng, description='x shift',
                                  value=0., readout_format='.4f',
                                  continuous_update=continuous_update)
    y_shift = widgets.FloatSlider(min=-shift_rng, max=+shift_rng,
                                  step=.05*shift_rng, description='y shift',
                                  value=0., readout_format='.4f',
                                  continuous_update=continuous_update)

    def make_slider_callback(changer):
        def slider_update(change):
            changer.set(change['new'])

            app_mgr.refresh_gui(build='update')

            if oa_ray is not None:
                oa_ray.update_data(build='rebuild')

            # apply changes to fans and grids
            for ray_data in ray_data_items:
                # if oa_ray is not None:
                #     ray_data.image_pt_2d = oa_ray.t_abr
                ray_data.update_data(build='rebuild')

            # update and plot results
            for fig in figs:
                fig.clf()
                fig.plot()

        return slider_update

    alpha_tilt.observe(make_slider_callback(AttrChanger(mirror.decenter,
                                                        'euler', index=0)),
                       names='value')
    x_shift.observe(make_slider_callback(AttrChanger(mirror.decenter,
                                                     'dec', index=0)),
                    names='value')
    y_shift.observe(make_slider_callback(AttrChanger(mirror.decenter,
                                                     'dec', index=1)),
                    names='value')

    return alpha_tilt, x_shift, y_shift


class AttrChanger():
    """Changer built on an object/attribute pair. """
    def __init__(self, obj, attr, index=None):
        self.object = obj
        self.attr = attr
        self.index = index

    def get(self):
        value = getattr(self.object, self.attr, None)
        #print('AttrChanger.get:', self.attr, value, self.index)
        value = value if self.index is None else value[self.index]
        return value

    def set(self, value):
        if self.index is not None:
            seq = getattr(self.object, self.attr, None)
            seq[self.index] = value
            value = seq
        setattr(self.object, self.attr, value)
        #print('AttrChanger.set:', self.attr, value)