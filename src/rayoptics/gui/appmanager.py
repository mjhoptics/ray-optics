#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Lightweight manager class to connect a model+actions to windows

.. Created on Fri Aug 17 11:35:01 2018

.. codeauthor: Michael J. Hayford
"""
import logging

from collections import namedtuple

ModelInfo = namedtuple('ModelInfo', ['model', 'fct', 'args', 'kwargs'])
""" package of a model and a update function with args and kwargs

    Attributes:
        model: object associated with view update function
        fct: view update function, can be None
        args: list of fct arguments
        kwargs: list of fct keyword arguments

    Note:
        ``model`` is required; ``fct``, ``args`` and ``kwargs`` are optional
"""
ModelInfo.__new__.__defaults__ = (None, (), {})


class AppManager:
    """ Lightweight model/view manager class

    Lightweight manager class to manage connections between a model and
    a ui view.

    The main function of AppManager is refresh_gui(). This is called after
    a user input to the gui to update the model and call a refresh function
    for each ui view of that model.

    Attributes:
        model: the model of the currently active/frontmost ui view

            model is expected to respond to:

                - update_model()
                - name()

        gui_parent: the top level gui manager (optional)

            gui_parent is expected to implement:

                - refresh_app_ui()

        view_dict: keys are ui views, values are ModelInfo tuples

            view is expected to implement:

                - windowTitle() - only for debug logging

    """

    def __init__(self, model, gui_parent=None):
        self.model = self._tag_model(model)
        self.gui_parent = gui_parent
        self.view_dict = {}
        self.figures = []

    def _tag_model(self, model):
        if model is not None:
            if not hasattr(model, 'app_manager'):
                model.app_manager = self

    def set_model(self, model):
        self._tag_model(model)
        self.model = model

    def add_view(self, view, gui_hook, model_info):
        """ Add a new view and model tuple into dictionary

        Args:
            view: the ui view, used as a key
            gui_hook: instance of the GUI component to be refreshed
            model_info: instance of the ModelInfo tuple

        Returns:
            returns the input view
        """
        self._tag_model(model_info.model)
        self.view_dict[view] = gui_hook, model_info
        return view

    def add_figure(self, fig):
        """ Add a new figure to be updated at refresh_gui.

        Args:
            fig: the ui figure

        Returns:
            returns the input figure
        """
        self.figures.append(fig)
        return fig

    def delete_view(self, view):
        """ removes view from the view dictionary

        Args:
            view: view being closed by user
        """
        logging.debug('AppManager.delete_view: {}'.format(view.windowTitle()))
        del self.view_dict[view]

    def close_model(self, view_close_fct=None):
        """ close all ui views associated with the active model

        Args:
            view_close_fct: optional fct called on each closing view, with
                the view as an argument. This function, if used, should call
                delete_view itself.
        """
        cur_model = self.model
        view_items = list(self.view_dict.items())
        for view, info in view_items:
            mi = info[1]
            if mi.model == cur_model:
                if view_close_fct is None:
                    self.delete_view(view)
                else:
                    view_close_fct(view)
        delattr(cur_model, 'app_manager')
        self.model = None

    def refresh_gui(self, **kwargs):
        """ update the active model and refresh its dependent ui views """
        self.model.update_model()
        self.refresh_views(**kwargs)
        self.refresh_figures(**kwargs)

    def refresh_views(self, **kwargs):
        """ refresh the dependent ui views of the active model """
        if self.gui_parent is not None:
            self.gui_parent.refresh_app_ui()
        # traverse a copy of the view dict. this way we can delete any errant
        #  views we might find. we don't seem to be trapping closing a window
        #  via a close box under pyqt5
        for view, info in dict(self.view_dict).items():
            mi = info[1]
            if mi.model == self.model:
                if mi.fct is not None:
                    try:
                        mi.fct(*mi.args, **mi.kwargs)
                    except RuntimeError:
                        del self.view_dict[view]

    def refresh_figures(self, **kwargs):
        """ refresh the dependent ui views of the active model """
        # traverse a copy of the view dict. this way we can delete any errant
        #  views we might find. we don't seem to be trapping closing a window
        #  via a close box under pyqt5
        for fig in self.figures:
            fig.refresh(**kwargs)

    def on_view_activated(self, view):
        """ Makes the model associated with input view the active model

        Args:
            view: view becoming the active view
        """
        if view is not None:
            try:
                info = self.view_dict[view]
            except KeyError:
                logging.debug('view "%s" not in view_dict',
                              view.windowTitle())
            else:
                mi = info[1]
                model = mi.model
                if model:
                    logging.debug("on_view_activated: %s, %s" %
                                  (model.name(), view.windowTitle()))
                if model and model != self.model:
                    self.model = model
                    self.refresh_views()
                    logging.debug("switch model to", model.name())

    def sync_light_or_dark(self, is_dark):
        """ Tells views to update to a light or dark color scheme """
        for view, _ in dict(self.view_dict).items():
            try:
                view.sync_light_or_dark(is_dark)
            except AttributeError:
                pass
