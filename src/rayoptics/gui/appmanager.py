#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Lightweight manager class to connect a model+actions to windows

.. Created on Fri Aug 17 11:35:01 2018

.. codeauthor: Michael J. Hayford
"""
import logging

from collections import namedtuple

logger = logging.getLogger(__name__)

ModelInfo = namedtuple('ModelInfo', ['model', 'fct', 'args', 'kwargs'])
ModelInfo.__new__.__defaults__ = (None, (), {})
ModelInfo.model.__doc__ = "object associated with view update function"
ModelInfo.fct.__doc__ = "view update function, can be None"
ModelInfo.args.__doc__ = "list of fct arguments"
ModelInfo.kwargs.__doc__ = "list of fct keyword arguments"


""" package of a model and a update function with args and kwargs

    Attributes:
        model: object associated with view update function
        fct: view update function, can be None
        args: list of fct arguments
        kwargs: list of fct keyword arguments

    Note:
        ``model`` is required; ``fct``, ``args`` and ``kwargs`` are optional
"""



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

        model_filenames: keys are models, values are Paths or None

    """

    def __init__(self, model, gui_parent=None, filename=None):
        self.view_dict = {}
        self.figures = []
        self.model_filenames = {}
        self.gui_parent = gui_parent
        self.set_model(model, filename=filename)

    def _tag_model(self, model):
        if model is not None:
            if not hasattr(model, 'app_manager'):
                model.app_manager = self

    def set_model(self, model, filename=None):
        self.model = model
        self._tag_model(model)
        if model is not None:
            self.model_filenames[model] = filename

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
        logger.debug(f"AppManager.delete_view: {view.windowTitle()} ["
                     f"x: {view.x()}, y: {view.y()}, "
                     f"w: {view.width()}, h: {view.height()}]")
        self.view_dict.pop(view, None)

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
        if cur_model is not None:
            delattr(cur_model, 'app_manager')
            self.model = None
            if cur_model in self.model_filenames:
                del self.model_filenames[cur_model]

    def refresh_gui(self, **kwargs):
        """ update the active model and refresh its dependent ui views """
        if self.model is not None:
            self.model.update_model(**kwargs)
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
                        mi.fct(*mi.args, **{**kwargs, **mi.kwargs})
                    except RuntimeError:
                        del self.view_dict[view]

    def refresh_figures(self, **kwargs):
        """ refresh the dependent ui views of the active model """
        # traverse a copy of the view dict. this way we can delete any errant
        #  views we might find. we don't seem to be trapping closing a window
        #  via a close box under pyqt5
        for fig in self.figures:
            fig.refresh(**kwargs)

    def listobj_str(self):
        o_str = f"active model_id={id(self.model)}:\n"

        for view, info in dict(self.view_dict).items():
            widget, mi = info
            if mi.model is not None:
                lens_title = "model_id=" + str(id(mi.model))
            else:
                lens_title = "No model"
            win_title = view.windowTitle()
            o_str += f'window: "{win_title}", {lens_title}\n'
                
        for i, fig in enumerate(self.figures):
             o_str += f"fig {i}: {type(fig).__name__}\n"

        return o_str        

    def on_view_activated(self, view):
        """ Makes the model associated with input view the active model

        Args:
            view: view becoming the active view
        """
        if view is not None:
            try:
                info = self.view_dict[view]
            except KeyError:
                logger.debug('view "%s" not in view_dict',
                             view.windowTitle())
            else:
                mi = info[1]
                model = mi.model
                if model:
                    logger.debug("on_view_activated: %s, %s" %
                                 (model.name(), view.windowTitle()))
                if model and model != self.model:
                    self.model = model
                    self.refresh_views()
                    logger.debug("switch model to", model.name())

    def sync_light_or_dark(self, is_dark):
        """ Tells views to update to a light or dark color scheme """
        for view, _ in dict(self.view_dict).items():
            try:
                view.sync_light_or_dark(is_dark)
            except AttributeError:
                pass
