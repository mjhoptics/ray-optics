#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2017 Michael J. Hayford
""" Module for different surface profile shapes

Created on Tue Aug  1 13:18:57 2017

@author: Michael J. Hayford
"""

import numpy as np
from math import sqrt, copysign, sin, atan2
from util.misc_math import normalize


class SurfaceProfile:
    def __init__(self):
        self.type = ''

    def __repr__(self):
        return "Profile(%r)" % self.type

    def update(self):
        pass

    def f(self, p):
        pass

    def normal(self, p):
        pass

    def sag(self, x, y):
        pass

    def profile(self, sd, dir=1, steps=6):
        pass

    def intersect(self, p0, d, eps=1.0e-12):
        p = p1 = p0
        s1 = -self.f(p1)/np.dot(d, self.normal(p1))
        delta = abs(s1)
        while delta > eps:
            p = p1 + s1*d
            s2 = s1 - self.f(p)/np.dot(d, self.normal(p))
            delta = abs(s2 - s1)
            # print(s1, s2, delta)
            s1 = s2
        return s1, p


class Spherical(SurfaceProfile):
    """ Spherical surface profile parameterized by curvature. """
    def __init__(self, c=0.0):
        self.type = 'Sphere'
        self.cv = c

    def __str__(self):
        return self.type + " " + str(self.cv)

    def __repr__(self):
        return "Profile(Spherical: c=%r)" % self.cv

    def copyFrom(self, other):
        dispatch[self.__class__, other.__class__](self, other)

    def copyDataFrom(self, other):
        self.cv = other.cv

    def mutate(self, new_profile):
        Spherical.copyDataFrom(new_profile, self)

    def normal(self, p):
        return normalize(np.array(
                [-self.cv*p[0], -self.cv*p[1], 1.0-self.cv*p[2]]))

    def f(self, p):
        return p[2] - 0.5*self.cv*(np.dot(p, p))

    def sag(self, x, y):
        if self.cv != 0.0:
            r = 1/self.cv
            adj = sqrt(r*r - x*x - y*y)
            return r*(1 - abs(adj/r))
        else:
            return 0

    def profile(self, sd, dir=1, steps=6):
        prf = []
        if self.cv != 0.0:
            r = 1/self.cv
            adj = sqrt(r*r - sd*sd)
            ang = atan2(sd, adj)
            da = dir*copysign(ang/steps, self.cv)
            sa = sin(da)
            ca = sqrt(1.0 - sa*sa)
            sab = sb = -dir*(sd/r)
            cab = cb = abs(adj/r)
            for i in range(2*steps):
                prf.append([r*(1-cab), r*sab])
                # print(i, r*(1-cab), r*sab)
                sab = sa*cb + sb*ca
                cab = ca*cb - sa*sb
                sb = sab
                cb = cab
            prf.append([r*(1-cab), r*sab])
            # print(2*steps, r*(1-cab), r*sab)
        else:
            prf.append([0, -dir*sd])
            prf.append([0, dir*sd])
        return prf


class Conic(SurfaceProfile):
    """ Conic surface profile parameterized by curvature and conic constant.

    Conics produced for conic constant values:
        cc > 0.0: oblate spheroid
        cc = 0.0: sphere
        cc < 0.0 and > -1.0: ellipsoid
        cc = -1.0: paraboloid
        cc < -1.0: hyperboloid
    """
    def __init__(self, c=0.0):
        self.type = 'Conic'
        self.cv = c
        self.cc = 0.0

    def __str__(self):
        return self.type + " " + str(self.cv) + " " + str(self.cc)

    def __repr__(self):
        return "Profile(Conic: c=%r, cc=%r)" % (self.cv, self.cc)

    def copyFrom(self, other):
        dispatch[self.__class__, other.__class__](self, other)

    def copyDataFrom(self, other):
        self.cv = other.cv
        self.cc = other.cc

    def normal(self, p):
        return normalize(np.array(
                [-self.cv*p[0],
                 -self.cv*p[1],
                 1.0-(self.cc+1.0)*self.cv*p[2]]))

    def f(self, p):
        return p[2] - 0.5*self.cv*(p[0]*p[0] +
                                   p[1]*p[1] +
                                   (self.cc+1.0)*p[2]*p[2])

    def sag(self, x, y):
        r2 = x*x + y*y
        z = self.cv*r2/(1. + sqrt(1. - (self.cc+1.0)*self.cv*self.cv*r2))
        return z

    def profile(self, sd, dir=1, steps=6):
        prf = []
        if self.cv != 0.0:
            delta = dir*sd/steps
            y = -dir*sd
            z = self.sag(0., y)
            prf.append([z, y])
            for i in range(2*steps):
                y += delta
                z = self.sag(0., y)
                prf.append([z, y])
                # print(i, z, y)
        else:
            prf.append([0, -dir*sd])
            prf.append([0, dir*sd])
        return prf


class EvenPolynomial(SurfaceProfile):
    """ Even Polynomial asphere up to 20th order, on base conic. """
    def __init__(self, c=0.0):
        self.type = 'EvenPolynomial'
        self.cv = c
        self.cc = 0.0
        self.coefs = []
        self.max_nonzero_coef = 0
        self.coef2 = 0.0
        self.coef4 = 0.0
        self.coef6 = 0.0
        self.coef8 = 0.0
        self.coef10 = 0.0
        self.coef12 = 0.0
        self.coef14 = 0.0
        self.coef16 = 0.0
        self.coef18 = 0.0
        self.coef20 = 0.0

    def copyFrom(self, other):
        dispatch[self.__class__, other.__class__](self, other)

    def copyDataFrom(self, other):
        self.cv = other.cv
        self.cc = other.cc
        self.coef2 = other.coef2
        self.coef4 = other.coef4
        self.coef6 = other.coef6
        self.coef8 = other.coef8
        self.coef10 = other.coef10
        self.coef12 = other.coef12
        self.coef14 = other.coef14
        self.coef16 = other.coef16
        self.coef18 = other.coef18
        self.coef20 = other.coef20

    def gen_coef_list(self):
        self.coefs = []
        self.coefs.append(self.coef2)
        self.coefs.append(self.coef4)
        self.coefs.append(self.coef6)
        self.coefs.append(self.coef8)
        self.coefs.append(self.coef10)
        self.coefs.append(self.coef12)
        self.coefs.append(self.coef14)
        self.coefs.append(self.coef16)
        self.coefs.append(self.coef18)
        self.coefs.append(self.coef20)
        self.max_nonzero_coef = -1
        for i, c in enumerate(self.coefs):
            if c != 0.0:
                self.max_nonzero_coef = i
        self.max_nonzero_coef += 1

    def update(self):
        self.gen_coef_list()

    def normal(self, p):
        # sphere + conic contribution
        r2 = p[0]*p[0] + p[1]*p[1]
        e = self.cv/sqrt(1. - (self.cc+1.0)*self.cv*self.cv*r2)

        # polynomial asphere contribution
        e_asp = 0.0
        r_pow = 1.0
        c_coef = 2.0
        for i in range(self.max_nonzero_coef):
            e_asp += c_coef*self.coefs[i]*r_pow
            c_coef += 2.0
            r_pow *= r2

        e_tot = e + e_asp
        return normalize(np.array([-e_tot*p[0], -e_tot*p[1], 1.0]))

    def sag(self, x, y):
        r2 = x*x + y*y
        # sphere + conic contribution
        z = self.cv*r2/(1. + sqrt(1. - (self.cc+1.0)*self.cv*self.cv*r2))

        # polynomial asphere contribution
        z_asp = 0.0
        r_pow = r2
        for i in range(self.max_nonzero_coef):
            z_asp += self.coefs[i]*r_pow
            r_pow *= r2

        z_tot = z + z_asp
        return z_tot

    def f(self, p):
        return p[2] - self.sag(p[0], p[1])

    def profile(self, sd, dir=1, steps=21):
        if steps < 21:
            steps = 21
        prf = []
        if self.max_nonzero_coef > 0 or self.cv != 0.0:
            delta = dir*sd/steps
            y = -dir*sd
            z = self.sag(0., y)
            prf.append([z, y])
            for i in range(2*steps):
                y += delta
                z = self.sag(0., y)
                prf.append([z, y])
                # print(i, z, y)
        else:
            prf.append([0, -dir*sd])
            prf.append([0, dir*sd])
        return prf


dispatch = {
  (Spherical, Spherical): Spherical.copyDataFrom,
  (Spherical, Conic): Spherical.copyDataFrom,
  (Spherical, EvenPolynomial): Spherical.copyDataFrom,
  (Conic, Spherical): Spherical.copyDataFrom,
  (Conic, Conic): Conic.copyDataFrom,
  (Conic, EvenPolynomial): Conic.copyDataFrom,
  (EvenPolynomial, Spherical): Spherical.copyDataFrom,
  (EvenPolynomial, Conic): Conic.copyDataFrom,
  (EvenPolynomial, EvenPolynomial): EvenPolynomial.copyDataFrom,
}


def mutate_profile(old_profile, new_profile_type):
    new_profile = eval(new_profile_type+'()')
    new_profile.copyFrom(old_profile)
    return new_profile


def test():
    s1 = Spherical(0.0)
    s2 = Spherical(0.1)

    dir0 = np.array([0., 0., 1.])
    v1 = normalize(np.array([1., 1., 1.]))
    p0 = np.array([0., 0., -1.])
    p1 = np.array([0., 1., -1.])

    p0s1 = s1.intersect(p0, dir0)
    print(p0s1)
    p1s1 = s1.intersect(p1, dir0)
    print(p1s1)

    p0s2 = s2.intersect(p0, dir0)
    print(p0s2)
    p1s2 = s2.intersect(p1, dir0)
    print(p1s2)
    dir_p1s2 = s2.normal(p1s2[1])
    print(dir_p1s2)

    print("pass")


if __name__ == "__main__":
    test()
