#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2020 Michael J. Hayford
"""Base class for Interfaces

.. Created on Sat Jun 13 22:04:27 2020

.. codeauthor: Michael J. Hayford
"""


from enum import Enum, auto


class InteractionMode(Enum):
    """ enum for different interact_mode specifications

    Retained to restore old files

    .. deprecated:: 0.4.5
    """
    Transmit = auto()  #: propagate in transmission at this interface
    Reflect = auto()   #: propagate in reflection at this interface


class Interface:
    """Basic part of a sequential model

    The :class:`~sequential.SequentialModel` is a sequence of Interfaces and
    Gaps. The Interface class is a boundary between two adjacent Gaps and
    their associated media. It specifies several methods that must be
    implemented to model the optical behavior of the interface.

    The Interface class addresses the following use cases:

        - support for ray intersection calculation during ray tracing
            - interfaces can be tilted and decentered wrt the adjacent gaps
        - support for getting and setting the optical power of the interface
        - support for various optical properties, i.e. does it reflect or
          transmit
        - supports a basic idea of size, the max_aperture

    Attributes:
        interact_mode: 'transmit' | 'reflect'
        delta_n: refractive index difference across the interface
        decenter: DecenterData for the interface, if specified
        max_aperture: the maximum aperture radius on the interface
    """
    def __init__(self, interact_mode='transmit', delta_n=0.0,
                 max_ap=1.0, decenter=None, phase_element=None, **kwargs):
        self.interact_mode = interact_mode
        self.delta_n = delta_n
        self.decenter = decenter
        self.max_aperture = max_ap
        if phase_element is not None:
            self.phase_element = phase_element

    def update(self):
        if self.decenter is not None:
            self.decenter.update()

    def interface_type(self):
        return type(self).__name__

    def sync_to_restore(self, opt_model):
        if not hasattr(self, 'max_aperture'):
            self.max_aperture = 1.0
        if hasattr(self, 'interact_mode'):
            # don't know why I need to test for the InteractionMode
            #  enum like this, or have to compare enum values, but
            #  that's what works...
            if isinstance(self.interact_mode, Enum):
                imode = self.interact_mode.value
                if imode == InteractionMode.Reflect.value:
                    self.interact_mode = 'reflect'
                elif imode == InteractionMode.Transmit.value:
                    self.interact_mode = 'transmit'
        if hasattr(self, 'refract_mode'):  # really old models
            if self.refract_mode == 'REFL':
                self.interact_mode = 'reflect'
            else:
                self.interact_mode = 'transmit'
            delattr(self, 'refract_mode')

    @property
    def profile_cv(self):
        return 0.0

    def set_optical_power(self, pwr, n_before, n_after):
        pass

    def surface_od(self):
        pass

    def set_max_aperture(self, max_ap):
        """ max_ap is the max aperture radius """
        self.max_aperture = max_ap

    def intersect(self, p0, d, eps=1.0e-12):
        pass

    def normal(self, p):
        pass

    def phase(self, pt, d_in, normal, wl):
        if hasattr(self, 'phase_element'):
            return self.phase_element.phase(pt, d_in, normal, wl=wl)

    def apply_scale_factor(self, scale_factor):
        self.max_aperture *= scale_factor
        if self.decenter:
            self.decenter.apply_scale_factor(scale_factor)
