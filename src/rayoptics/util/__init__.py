""" package supplying utility functions for math and numpy support

    The :mod:`~rayoptics.util` subpackage provides miscellaneous functions for
    geometric calculations, color calculations and anything else that doesn't
    have an obvious home. These include:

        - miscellaneous math functions,
          :mod:`~.misc_math`
          :mod:`~.line_intersection`
        - support for color handling, :mod:`~.colors` :mod:`~.colour_system`
        - spectral line conversion with :func:`~.spectral_lines.get_wavelength`
          in :mod:`~.spectral_lines`
        - a 2D dict with M x N keys, :mod:`~.dict2d`
"""

import importlib
import logging

def str_to_class(module_name:str, class_name:str, **kwargs):
    """Return a class instance from a string reference"""
    try:
        module_ = importlib.import_module(module_name)
        try:
            class_obj = getattr(module_, class_name)
            class_ = class_obj(**kwargs)
        # except AttributeError:
        except Exception as err:
            logging.error(f'Class "{class_name}" does not exist')
    except ImportError:
        logging.error(f'Module "{module_name}" does not exist')
    return class_ or None
