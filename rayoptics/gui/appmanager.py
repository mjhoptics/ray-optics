#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Lightweight manager class to connect a model+actions to windows

Created on Fri Aug 17 11:35:01 2018

@author: Michael J. Hayford
"""
import logging

from collections import namedtuple

ModelInfo = namedtuple('ModelInfo', ['model', 'fct', 'arg'])


class AppManager:
    def __init__(self, model):
        self.model = model
        self.wnd_dict = {}

    def add_window(self, window, model_info):
            self.wnd_dict[window] = model_info
            return window

    def delete_window(self, window):
            del self.wnd_dict[window]

    def refresh_gui(self):
        self.model.update_model()
        for model_info in self.wnd_dict.values():
            if model_info.model == self.model:
                model_info.fct(model_info.arg)

    def on_window_activated(self, window):
        if window is not None:
            try:
                model_info = self.wnd_dict[window]
            except KeyError:
                logging.debug('Window "%s" not in wnd_dict',
                              window.windowTitle())
            else:
                model = model_info[0]
                logging.debug("on_window_activated: %s, %s" %
                              (model.name(), window.windowTitle()))
                if model and model != self.model:
                    self.model = model
                    logging.debug("switch model to", model.name())
