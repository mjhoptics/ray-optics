#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2020 Michael J. Hayford
"""Base class for Interfaces

.. Created on Sat Jun 13 22:04:27 2020

.. codeauthor: Michael J. Hayford
"""
import numpy as np
from numpy import sqrt
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
        interact_mode: 'transmit' | 'reflect' | 'dummy'
        delta_n: refractive index difference across the interface
        decenter: :class:`~rayoptics.elem.surface.DecenterData` for the interface, if specified
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

    def listobj_str(self):
        o_str = (f"{self.interact_mode}   delta n={self.delta_n}   "
                 f"max aperture={self.max_aperture}\n")
        if hasattr(self, 'phase_element') and self.phase_element is not None:
            o_str += self.phase_element.listobj_str()
        if self.decenter is not None:
            o_str += self.decenter.listobj_str()

        return o_str

    def update(self):
        if self.decenter is not None:
            self.decenter.update()

    def interface_type(self):
        return type(self).__name__

    def ifc_token(self):
        tkn = ''
        if self.interact_mode == 'transmit':
            tkn = 'i'
        elif self.interact_mode == 'reflect':
            tkn = 'r'
        elif self.interact_mode == 'dummy':
            tkn = 'd'
        return tkn

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

    def edge_pt_target(self, rel_dir):
        """ Get a target for ray aiming to aperture boundaries.
        
        The main use case for this function is iterating a ray to the internal 
        edge of a surface. 

        Although `rel_dir` is given as a 2d vector, in practice only the 4 
        quadrant axes are handled in the implementation, a 1D directional 
        search along a coordinate axis.

        Args:
            rel_dir: 2d vector encoding coord axis and direction for edge sample
        
        Returns:
            edge_pt: intersection point of rel_dir with the aperture boundary
        """
        edge_pt = np.array([self.max_aperture*rel_dir[0], 
                            self.max_aperture*rel_dir[1]])
        return edge_pt

    def point_inside(self, x: float, y: float, fuzz: float = 1e-5) -> bool:
        """ Returns True if the point (x, y) is inside the clear aperture. 
        
        Args:
            x: x coodinate of the test point
            y: y coodinate of the test point
            fuzz: tolerance on test pt/aperture comparison, 
                  i.e. pt fuzzy <= surface_od
        """
        return sqrt(x*x + y*y) <= self.max_aperture + fuzz

    def set_max_aperture(self, max_ap):
        """ max_ap is the max aperture radius """
        self.max_aperture = max_ap

    def get_y_aperture_extent(self):
        """ default behavior is returning +/-max_aperture """
        od = [-self.max_aperture, self.max_aperture]
        return od

    def intersect(self, p0, d, eps=1.0e-12, z_dir=1):
        ''' Intersect an :class:`~.Interface`, starting from an arbitrary point.

        Args:
            p0:  start point of the ray in the interface's coordinate system
            d:  direction cosine of the ray in the interface's coordinate system
            z_dir: +1 if propagation positive direction, -1 if otherwise
            eps: numeric tolerance for convergence of any iterative procedure

        Returns:
            tuple: distance to intersection point *s1*, intersection point *p*

        Raises:
            :exc:`~rayoptics.raytr.traceerror.TraceMissedSurfaceError`
        '''
        pass

    def normal(self, p):
        """Returns the unit normal of the interface at point *p*. """
        pass

    def phase(self, pt, in_dir, srf_nrml, ifc_cntxt):
        z_dir, wvl, n_in, n_out, interact_mode = ifc_cntxt
        """Returns a diffracted ray direction and phase increment.

        Args:
            pt: point of incidence in :class:`~.Interface` coordinates
            in_dir: direction cosine of incident ray
            srf_nrml: :class:`~.Interface` surface normal at pt
            ifc_cntxt: a tuple containing
            
                z_dir: -1 if after an odd # of reflections, +1 otherwise
                wl: wavelength in nm for ray, defaults to ref_wl
                n_in: refractive index preceding the interface
                n_out: refractive index following the interface
                interact_mode: 'transmit' or 'reflect'

        Returns:
            (**out_dir, dW**)

            - out_dir: direction cosine of the out going ray
            - dW: phase added by diffractive interaction
        """
        if hasattr(self, 'phase_element'):
            return self.phase_element.phase(pt, in_dir, srf_nrml, ifc_cntxt)

    def apply_scale_factor(self, scale_factor):
        self.max_aperture *= abs(scale_factor)
        if self.decenter:
            self.decenter.apply_scale_factor(scale_factor)

    def update_following_reflection(self):
        """ Notification that incident light is following reflection. """
        pass
