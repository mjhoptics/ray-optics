.. currentmodule:: rayoptics.qtgui

.. _qt5_app:

Qt5 app version of rayoptics
============================

    The :mod:`rayoptics.qtgui` subpackage provides a desktop app, :obj:`rayopticsapp`, that runs under Anaconda. It also provides a series of higher level interfaces used by rayoptics. These include:

        - an interface that hosts matplotlib graphics
        - a table grid for numeric model displays (template-based)
        - docking panel support for python objects
        - iPython console window (desktop app only)

Running the app
===============

A desktop application is installed as part of ``rayoptics``. It is invoked by running ``rayoptics`` at the command line.

.. code::

   MacBook-Pro:docs Mike$ rayoptics
