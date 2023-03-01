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


class DiffractionGrating:
    """ Linear (ruled) diffraction grating.
    
    The phase calculation is patterned after Ludwig, 
    `Generalized grating ray-tracing equations <https://doi.org/10.1364/JOSA.63.001105>`_.
    
    Attributes:
        
        grating_normal: grating generation surface normal (**G**)
        grating_lpmm: the grating frequency in lines/mm
        grating_freq_um: grating frequency in lines per micrometer
        order: integer diffraction order used in phase calculation
        interact_mode: 'transmit'|'reflect'
        label: string description
    
    """
    def __init__(self, label='', order=1, grating_normal=None,
                 grating_freq_um=1.0, grating_lpmm=None,
                 interact_mode='transmit'):
        self.label = label
        if grating_normal is None:
            self.grating_normal = np.array([0., 1., 0.])
        else:
            self.grating_normal = grating_normal

        if grating_lpmm is not None:
            self.grating_lpmm = grating_lpmm
        else:
            self.grating_lpmm = 1/(grating_freq_um * 1000)
        self.order = order

        self.interact_mode = interact_mode
        self.debug_output = False

    @property
    def grating_lpmm(self):
        """ the grating spacing in lines/mm. """
        return self._grating_lpmm

    @grating_lpmm.setter
    def grating_lpmm(self, grating_lpmm):
        self._grating_lpmm = grating_lpmm
        self._grating_spacing_nm = 1e6/grating_lpmm

    @property
    def grating_freq_um(self):
        return self._grating_lpmm/1000

    @grating_freq_um.setter
    def grating_freq_um(self, grating_freq_um):
        self.grating_lpmm = grating_freq_um * 1000

    def listobj_str(self):
        if len(self.label) == 0:
            label = 'grating'
        else:
            label = self.label
        o_str = (f"{label}: {self.interact_mode} order: {self.order}, "
                 f"grating_lpmm: {self.grating_lpmm}\n"
                 f"grating_normal: {self.grating_normal}\n")
        return o_str

    def phase(self, pt, in_dir, srf_nrml, ifc_cntxt):
        return self.phase_ludwig(pt, in_dir, srf_nrml, ifc_cntxt)

    def phase_ludwig(self, pt, in_dir, srf_nrml, ifc_cntxt):
        z_dir, wvl, n_in, n_out, interact_mode = ifc_cntxt
        refl = -1 if interact_mode == 'reflect' else 1
        normal = z_dir * normalize(srf_nrml)          # = R

        # grating ruling vector, P = G x R
        P = np.cross(self.grating_normal, normal)
        # ruling separation vector, D = R x P
        D = normalize(np.cross(normal, P))

        mu = n_in / n_out
        T = refl*(wvl * self.order)/(self._grating_spacing_nm * n_out)

        in_cosI = np.dot(in_dir, normal)
        V = mu * in_cosI
        W = mu**2 - 1 + T**2 - 2*mu*T*(np.dot(D, in_dir))
        
        result = np.sqrt(V**2 - W)
        Q1 = result - V
        Q2 = -result - V
        if interact_mode == 'transmit':
            Q = max(Q1, Q2)
        elif interact_mode == 'reflect':
            Q = min(Q1, Q2)

        out_dir = mu*in_dir - T*D + Q*normal
        # `out_dir` is the unit vector in the diffracted ray direction.
        # The `l` and `m` components are correct in the unit circle.
        # The `n` component needs to be adjusted to fall on the unit hemisphere.
        # The sign of the original z-component is applied to the result.
        out_dir[2] = np.copysign(sqrt(1 - out_dir[0]**2 - out_dir[1]**2),
                                 out_dir[2])

        # calculate path difference in wavelengths introduced by grating. 
        in_sinI = sqrt(1 - in_cosI**2)
        out_cosI = np.dot(out_dir, normal)
        out_sinI = sqrt(1 - out_cosI**2)
        dW = (self._grating_spacing_nm/wvl) * (n_in*in_sinI + refl*n_out*out_sinI)

        if self.debug_output:
            from numpy.linalg import norm
            print(f"{interact_mode}: z_dir={z_dir}, {n_in:4.2f}, {n_out:4.2f}")
            print(f"in_dir: {in_dir}")
            print(f"R: {normal},  G: {self.grating_normal}")
            print(f"P = G x R: {P},  D = R x P: {D}")
            print(f"T={T:8.4f}, V={V:8.4f}, W={W:8.4f}")
            print(f"Q={Q:8.4f}, Q1={Q1:8.4f}, Q2={Q2:8.4f}")
            print(f"out_dir: {out_dir}, len={norm(out_dir):8.6f}")

        return out_dir, dW

    def phase_welford(self, pt, in_dir, srf_nrml, ifc_cntxt):
        z_dir, wvl, n_in, n_out, interact_mode = ifc_cntxt
        refl = -1 if interact_mode == 'reflect' else 1

        normal = z_dir * normalize(srf_nrml)          # = R
        cosI = np.dot(in_dir, normal)

        T = refl*(wvl * self.order)/(self._grating_spacing_nm * n_out)

        out_dir = in_dir.copy()
        out_dir[0] = in_dir[0]
        out_dir[1] = in_dir[1] - T
        out_dir[2] = (in_dir[2] - cosI 
                      + sqrt(cosI**2 + 2*in_dir[1]*T - T**2))

        # The `l` and `m` components are correct in the unit circle.
        # The `n` component needs to be adjusted to fall on the unit hemisphere.
        # The sign of the original z-component is applied to the result.
        out_dir[2] = np.copysign(sqrt(1 - out_dir[0]**2 - out_dir[1]**2),
                                 refl)

        # calculate path difference in wavelengths introduced by grating. 
        in_sinI = sqrt(1 - cosI**2)
        out_cosI = np.dot(out_dir, normal)
        out_sinI = sqrt(1 - out_cosI**2)
        dW = (self._grating_spacing_nm/wvl) * (n_in*in_sinI + refl*n_out*out_sinI)

        if self.debug_output:
            from numpy.linalg import norm
            print(f"{interact_mode}: z_dir={z_dir}, {n_in:4.2f}, {n_out:4.2f}")
            print(f"in_dir: {in_dir}")
            print(f"R: {normal},  G: {self.grating_normal}")
            print(f"T={T:8.4f}")
            print(f"out_dir: {out_dir}, len={norm(out_dir):8.6f}")

        return out_dir, dW


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

    def listobj_str(self):
        if len(self.label) == 0:
            label = 'doe'
        else:
            label = self.label
        o_str = f"{label}: {self.phase_fct.__name__}\n"
        o_str += f"coefficients: {self.coefficients}\n"
        o_str += f"ref wl: {self.ref_wl}nm  order: {self.order}\n"
        return o_str

    def phase(self, pt, in_dir, srf_nrml, ifc_cntxt):
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
        z_dir, wvl, n_in, n_out, interact_mode = ifc_cntxt
        order = self.order
        normal = normalize(srf_nrml)
        inc_dir = in_dir
        if n_in != 1.0:
            inc_dir = rt.bend(in_dir, srf_nrml, n_in, 1)
        in_cosI = np.dot(inc_dir, normal)
        mu = 1.0 if wvl is None else wvl/self.ref_wl
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
            out_dir = rt.bend(out_dir, srf_nrml, 1, n_out)

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

    def listobj_str(self):
        if len(self.label) == 0:
            label = 'hoe'
        else:
            label = self.label
        o_str = f"{label}: ref wl: {self.ref_wl}nm\n"
        o_str += (f"ref_pt: {self.ref_pt[0]:12.5g} {self.ref_pt[1]:12.5g} "
                  f"{self.ref_pt[2]:12.5g}   virtual: {self.ref_virtual}\n")
        o_str += (f"obj_pt: {self.obj_pt[0]:12.5g} {self.obj_pt[1]:12.5g} "
                  f"{self.obj_pt[2]:12.5g}   virtual: {self.obj_virtual}\n")
        return o_str

    def phase(self, pt, in_dir, srf_nrml, ifc_cntxt):
        z_dir, wvl, n_in, n_out, interact_mode = ifc_cntxt
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
        mu = 1.0 if wvl is None else wvl/self.ref_wl
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
