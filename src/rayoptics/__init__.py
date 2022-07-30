# -*- coding: utf-8 -*-
""" The **ray-optics** geometrical ray tracing and optical modeling and
    analysis package

    The optical model is contained in the :mod:`~.optical` subpackage. It is
    supported by the following subpackages:

        - :mod:`~.optical`: OpticalModel and OpticalBench import
        - :mod:`~.seq`: support for the Sequential model
        - :mod:`~.elem`: Element and PartTree models, geometry support, e.g. 
          for profiles
        - :mod:`~.oprops`: optical property and actions
        - :mod:`~.parax`: support for paraxial optical design
        - :mod:`~.raytr`: support for ray tracing and analysis

        - :mod:`~.codev`: handles import of CODE V .seq files
        - :mod:`~.zemax`: handles import of Zemax .zmx files

        - :mod:`opticalglass`: this package interfaces with glass manufacturer
          optical data and the RefractiveIndex.Info website

    The :mod:`~.gui` subpackage is a layer that implements the platform
    neutral part of the graphical user interface.

    The :mod:`~.mpl` subpackage implements a variety of plotting/charting
    based graphics using the :doc:`matplotlib <matplotlib:index>` package. The
    :mod:`~.qtgui` subpackage implements a desktop ui style application using
    **Qt** featuring interactive layouts and diagrams.

    The :mod:`~.util` subpackage provides a variety of different math
    and other miscellaneous calculations.
"""

from importlib.metadata import version

try:
    __version__ = version(__name__)
except:
    __version__ = 'unknown'


def listobj(obj):
    """ Print wrapper function for listobj_str() method of `obj`. 
    
    listobj() is designed to be used in scripting environments where detailed, 
    textual output is supported. It is a wrapper to a call of `listobj_str` on
    `obj`. 
    
    Classes may implement the `listobj_str` method that returns a string 
    containing a formatted description of the object. Multi-line strings are 
    allowed; each line should end with a newline character. Examples include 
    :meth:`.DecenterData.listobj_str` and :meth:`.EvenPolynomial.listobj_str`.
    """
    try:
        print(obj.listobj_str())
    except AttributeError:
        print(repr(obj))
