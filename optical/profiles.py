#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 13:18:57 2017

@author: Mike
"""

import numpy as np
from numpy.linalg import norm
from math import sqrt, copysign, sin, atan2


def normalize(v):
    length = norm(v)
    if length == 0.0:
        return v
    else:
        return v/length


class SurfaceProfile:
    def __init__(self):
        self.type = ''

    def __repr__(self):
        return "Profile(%r)" % self.type

    def f(self, p):
        pass

    def normal(self, p):
        pass

    def sag(self, x, y):
        pass

    def profile(self, sd, dir=1, steps=6):
        pass

    def intersect(self, p0, d, eps=1.0e-12):
        # s0 = -p0[2]/d[2]
        # p1 = np.array([p0[0]+s0*d[0], p0[1]+s0*d[1], 0.0])
        # p1 = p0
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


class EvenPolynomial(SurfaceProfile):
    def __init__(self, c=0.0):
        self.type = 'EvenPolynomial'
        self.cv = c
        self.cc = 0.0
        coefs = []
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
        self.coef4 = other.coef4
        self.coef6 = other.coef6
        self.coef8 = other.coef8
        self.coef10 = other.coef10
        self.coef12 = other.coef12
        self.coef14 = other.coef14
        self.coef16 = other.coef16
        self.coef18 = other.coef18
        self.coef20 = other.coef20


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
