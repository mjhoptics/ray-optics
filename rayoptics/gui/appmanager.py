#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Lightweight manager class to connect a model+actions to windows

Created on Fri Aug 17 11:35:01 2018

@author: Michael J. Hayford
"""
import logging

from collections import namedtuple

ModelInfo = namedtuple('ModelInfo', ['model', 'fct', 'args', 'kwargs'])
# model is required; fct, args and kwargs are optional
ModelInfo.__new__.__defaults__ = (None, (), {})


class AppManager:
    def __init__(self, model):
        """
        model is expected to respond to:
            update_model()
            name()

        view is expected to implement:
            windowTitle() - only for debug logging
        """
        self.model = model
        self.view_dict = {}

    def add_view(self, view, model_info):
            self.view_dict[view] = model_info
            return view

    def delete_view(self, view):
            del self.view_dict[view]

    def refresh_gui(self):
        self.model.update_model()
        for mi in self.view_dict.values():
            if mi.model == self.model:
                mi.fct(*mi.args, **mi.kwargs)

    def on_view_activated(self, view):
        if view is not None:
            try:
                model_info = self.view_dict[view]
            except KeyError:
                logging.debug('view "%s" not in view_dict',
                              view.windowTitle())
            else:
                model = model_info[0]
                logging.debug("on_view_activated: %s, %s" %
                              (model.name(), view.windowTitle()))
                if model and model != self.model:
                    self.model = model
                    logging.debug("switch model to", model.name())
