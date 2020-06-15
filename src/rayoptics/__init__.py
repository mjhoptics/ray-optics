# -*- coding: utf-8 -*-
""" The **ray-optics** geometrical ray tracing and optical modeling and
    analysis package

    The basic optical functionality is contained in the :mod:`~rayoptics.optical`
    subpackage. It is supported by the following subpackages:

        - :mod:`~rayoptics.codev`: handles import of CODE V .seq files
        - :mod:`opticalglass`: this package interfaces with glass manufacturer optical data

    The :mod:`~rayoptics.gui` subpackage is a layer that implements the platform
    neutral part of the graphical user interface.

    The :mod:`~rayoptics.mpl` subpackage implements a variety of plotting/charting
    based graphics using the **matplotlib** package. The :mod:`~rayoptics.qtgui`
    subpackage implements a desktop ui style application using **Qt** featuring
    interactive layouts and diagrams.

    The :mod:`~rayoptics.util` subpackage provides a variety of different math
    and other miscellaneous calculations.
"""

import pkg_resources

try:
    __version__ = pkg_resources.get_distribution(__name__).version
except:
    __version__ = 'unknown'

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
