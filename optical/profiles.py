#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 13:18:57 2017

@author: Mike
"""

import numpy as np
import numpy.linalg
from math import sqrt


def normalize(v):
    length = numpy.linalg.norm(v)
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

    def bend(self, d_in, normal, n_in, n_out):
        normal_len = numpy.linalg.norm(normal)
        cosI = np.dot(d_in, normal)/normal_len
        sinI_sqr = 1.0 - cosI*cosI
        n_cosIp = sqrt(n_out*n_out - n_in*n_in*sinI_sqr)
        alpha = n_cosIp - n_in*cosI
        d_out = (n_in*d_in + alpha*normal)/n_out
        return d_out


class Spherical(SurfaceProfile):
    def __init__(self, c=0.0):
        self.type = 'Sphere'
        self.cv = c

    def __str__(self):
        return self.type + " " + str(self.cv)

    def __repr__(self):
        return "Profile(Spherical: c=%r)" % self.cv

    def normal(self, p):
        return normalize(np.array(
                [-self.cv*p[0], -self.cv*p[1], 1.0-self.cv*p[2]]))

    def f(self, p):
        return p[2] - 0.5*self.cv*(np.dot(p, p))


class Conic(SurfaceProfile):
    def __init__(self, c=0.0):
        self.type = 'Conic'
        self.cv = c
        self.cc = 0.0

    def __str__(self):
        return self.type + " " + str(self.cv) + " " + str(self.cc)

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
