#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
""" Module for diffractive/holographic optical elements

    Classes that implement diffractive optics capabilities must implement
    the function phase() for use by the ray trace engine.

    The :class:`~.DiffractiveElement` and :class:`~.HolographicElement`
    implementations are patterned after Wang, et al, `Ray tracing and wave
    aberration calculation for diffractive optical elements
    <https://doi.org/10.1117/1.600780>`_

.. Created on Fri Jul  5 11:27:13 2019

.. codeauthor: Michael J. Hayford
"""


from math import sqrt
import numpy as np
import importlib

import rayoptics.raytr.raytrace as rt
from rayoptics.util.misc_math import normalize


def radial_phase_fct(pt, coefficients):
    """Evaluate the phase and slopes at **pt**

    Args:
        pt: 3d point of incidence in :class:`~.Interface` coordinates
        coefficients: list of even power radial phase coefficients,
                      e.g. r**2, r**4, ...

    Returns:
        (**dW, dWdX, dWdY**)

        - dW: phase added by diffractive interaction
        - dWdX: slope in x direction
        - dWdY: slope in y direction
    """
    x, y, z = pt
    r_sqr = x*x + y*y
    dW = 0
    dWdX = 0
    dWdY = 0
    for i, c in enumerate(coefficients):
        dW += c*r_sqr**(i+1)
        r_exp = r_sqr**(i)
        factor = 2*(i+1)
        dWdX += factor*c*x*r_exp
        dWdY += factor*c*y*r_exp
    return dW, dWdX, dWdY


class DiffractiveElement:
    """Container class for a phase fct driven diffractive optical element

    Attributes:
        phase_fct: fct the takes an input pt and returns phase and slope
        coefficients: list of coeficients for phase function
        ref_wl: wavelength in nm for phase measurement
        order: which diffracted order to calculate the phase for
        label: optical labeling for listing
    """

    def __init__(self, label='', coefficients=None, ref_wl=550., order=1,
                 phase_fct=None):
        self.label = label
        if coefficients is None:
            self.coefficients = []
        else:
            self.coefficients = coefficients
        self.ref_wl = ref_wl
        self.order = order
        self.phase_fct = phase_fct
        self.debug_output = False

    def __repr__(self):
        return (type(self).__name__ + '(label=' + repr(self.label) +
                ', coefficients=' + repr(self.coefficients) +
                ', ref_wl=' + repr(self.ref_wl) +
                ', order=' + repr(self.order) +
                ', phase_fct=' + repr(self.phase_fct) + ')')

    def __json_encode__(self):
        attrs = dict(vars(self))
        # Save model name and function name of phase_fct, so that fct can
        #  restored later (hopefully)
        del attrs['debug_output']
        del attrs['phase_fct']
        attrs['phase_fct_module'] = self.phase_fct.__module__
        attrs['phase_fct_name'] = self.phase_fct.__name__
        return attrs

    def __json_decode__(self, **attrs):
        module_name = attrs.pop('phase_fct_module')
        fct_name = attrs.pop('phase_fct_name')
        # try to import module and look up function - then assign to phase_fct
        mod = importlib.import_module(module_name)
        phase_fct = getattr(mod, fct_name)
        self.__init__(phase_fct=phase_fct, **attrs)

    def listobj(self):
        if len(self.label) == 0:
            label = 'doe'
        else:
            label = self.label
        print(f"{label}: {self.phase_fct.__name__}")
        print(f"coefficients: {self.coefficients}")
        print(f"ref wl: {self.ref_wl}nm  order: {self.order}")

    def phase(self, pt, in_dir, srf_nrml, z_dir, wl, n_in, n_out):
        """Returns a diffracted ray and phase increment.

        Args:
            pt: point of incidence in :class:`~.Interface` coordinates
            in_dir: incoming direction cosine of incident ray
            srf_nrml: :class:`~.Interface` surface normal at pt
            z_dir: -1 if after an odd # of reflections, +1 otherwise
            wl: wavelength in nm for ray, defaults to ref_wl
            n_in: refractive index preceding the interface
            n_out: refractive index following the interface

        Returns:
            (**out_dir, dW**)

            - out_dir: direction cosine of the out going ray
            - dW: phase added by diffractive interaction
        """
        order = self.order
        normal = normalize(srf_nrml)
        inc_dir = in_dir
        if n_in != 1.0:
            inc_dir = rt.bend(in_dir, srf_nrml, n_in, 1)
        in_cosI = np.dot(inc_dir, normal)
        mu = 1.0 if wl is None else wl/self.ref_wl
        dW, dWdX, dWdY = self.phase_fct(pt, self.coefficients)
        # print(wl, mu, dW, dWdX, dWdY)
        b = in_cosI + order*mu*(normal[0]*dWdX + normal[1]*dWdY)
        c = mu*(mu*(dWdX**2 + dWdY**2)/2 +
                order*(inc_dir[0]*dWdX + inc_dir[1]*dWdY))
        # pick the root based on z_dir
        Q = -b + z_dir*sqrt(b*b - 2*c)
        if self.debug_output:
            print('inc_dir:', inc_dir)
            scale_dir = in_dir
            scale_dir[2] = n_in
            scale_dir = normalize(scale_dir)
            print('scale_dir:', scale_dir)
            print("   mu        dW          dWdX          dWdY          b"
                  "            c           Q")
            print(f"{mu:6.3f} {dW:12.5g} {dWdX:12.5g} {dWdY:12.5g} {b:12.7g}"
                  f" {c:12.7g} {Q:12.7g}")
        out_dir = inc_dir + order*mu*(np.array([dWdX, dWdY, 0])) + Q*normal
        dW *= mu
        if n_in != 1.0:
            out_dir = rt.bend(out_dir, srf_nrml, 1, n_in)

        return out_dir, dW


class HolographicElement:
    """Two point hologram element. """
    def __init__(self, label=''):
        self.label = label
        self.ref_pt = np.array([0., 0., -1e10])
        self.ref_virtual = False
        self.obj_pt = np.array([0., 0., -1e10])
        self.obj_virtual = False
        self.ref_wl = 550.0

    def listobj(self):
        if len(self.label) == 0:
            label = 'hoe'
        else:
            label = self.label
        print(f"{label}: ref wl: {self.ref_wl}nm")
        print(f"ref_pt: {self.ref_pt[0]:12.5g} {self.ref_pt[1]:12.5g} "
              f"{self.ref_pt[2]:12.5g}   virtual: {self.ref_virtual}")
        print(f"obj_pt: {self.obj_pt[0]:12.5g} {self.obj_pt[1]:12.5g} "
              f"{self.obj_pt[2]:12.5g}   virtual: {self.obj_virtual}")

    def phase(self, pt, in_dir, srf_nrml, z_dir, wl, n_in, n_out):
        normal = normalize(srf_nrml)
        ref_dir = normalize(pt - self.ref_pt)
        if self.ref_virtual:
            ref_dir = -ref_dir
        ref_cosI = np.dot(ref_dir, normal)
        obj_dir = normalize(pt - self.obj_pt)
        if self.obj_virtual:
            obj_dir = -obj_dir
        obj_cosI = np.dot(obj_dir, normal)
        in_cosI = np.dot(in_dir, normal)
        mu = 1.0 if wl is None else wl/self.ref_wl
        b = in_cosI + mu*(obj_cosI - ref_cosI)
        refp_cosI = np.dot(ref_dir, in_dir)
        objp_cosI = np.dot(obj_dir, in_dir)
        ro_cosI = np.dot(ref_dir, obj_dir)
        c = mu*(mu*(1.0 - ro_cosI) + (objp_cosI - refp_cosI))
        # pick the root based on z_dir
        Q = -b + z_dir*sqrt(b*b - 2*c)
        out_dir = in_dir + mu*(obj_dir - ref_dir) + Q*normal
        dW = 0.
        return out_dir, dW
