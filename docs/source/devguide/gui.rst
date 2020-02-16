.. currentmodule:: rayoptics.gui

Generic GUI Overview
====================

   The :class:`appmanager.AppManager` manages connections between a model and a UI view. It is purposely lightweight, relying on window manager support for most windowing functions. The :func:`~appmanager.AppManager.refresh_gui` function provides a means of synchronizing model changes across all UI views associated with that model.

   The :class:`~appmanager.AppManager` supports multiple open models at the same time. The :attr:`~appmanager.AppManager.model` attribute corresponds to the model for the active UI view.

Matplotlib Support
==================

   Many forms of graphics can be supported *portably* by using the matplotlib package.

      - 2D lens layout graphics
      - Aberration plots for transverse ray aberrations and wavefront aberrations
      - Wavefront maps
      - Interactive |ybar| and |nubar| diagrams
