.. currentmodule:: rayoptics.qtgui

Qt application
==============

Qt5 app version of rayoptics
----------------------------

    The :mod:`~rayoptics.qtgui` subpackage provides a desktop app, :obj:`rayopticsapp`. It provides an integration of rayoptics with the Qt GUI toolkit. Capabilities include:

        - an interface that hosts matplotlib graphics
        - a table grid for numeric model displays (template-based)
        - docking panel support for python objects
        - iPython console window (desktop app only)

Running the app
---------------

A desktop application is installed as part of ``rayoptics``. It is invoked by running ``rayoptics`` at the command line.

.. code::

   > rayoptics

On a Windows machine, the ``rayoptics`` command will be located in a Scripts directory underneath the install directory. For example, if using a virtual environment named ``optics``, the command would be

.. code::

   > \optics\Scripts\rayoptics
