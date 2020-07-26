# -*- coding: utf-8 -*-
""" The **ray-optics** geometrical ray tracing and optical modeling and
    analysis package

    The optical model is contained in the :mod:`~.optical` subpackage. It is
    supported by the following subpackages:

        - :mod:`~.seq`: support for the Sequential model
        - :mod:`~.raytr`: support for ray tracing and analysis
        - :mod:`~.parax`: support for paraxial optical design
        - :mod:`~.oprops`: optical property and actions
        - :mod:`~.elem`: support for the Element model
        - :mod:`~.codev`: handles import of CODE V .seq files
        - :mod:`opticalglass`: this package interfaces with glass manufacturer
          optical data

    The :mod:`~.gui` subpackage is a layer that implements the platform
    neutral part of the graphical user interface.

    The :mod:`~.mpl` subpackage implements a variety of plotting/charting
    based graphics using the :doc:`matplotlib <matplotlib:index>` package. The
    :mod:`~.qtgui` subpackage implements a desktop ui style application using
    **Qt** featuring interactive layouts and diagrams.

    The :mod:`~.util` subpackage provides a variety of different math
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
