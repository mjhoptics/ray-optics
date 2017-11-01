#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 13:18:57 2017

@author: Mike
"""

from math import sqrt


def dot(a, b):
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]


def normalize(v):
    length = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])
    if length == 0.0:
        return v
    else:
        return [v[0]/length, v[1]/length, v[2]/length]


class SurfaceProfile:
    def __init__(self):
        self.type = ''

    def f(self, p):
        pass

    def normal(self, p):
        pass

    def intersect(self, p0, d, eps=1.0e-12):
        s0 = -p0[2]/d[2]
        p1 = [p0[0]+s0*d[0], p0[1]+s0*d[1], 0.0]
        p = p1
        s1 = -self.f(p1)/dot(d, self.normal(p1))
        delta = abs(s1)
        while delta > eps:
            p = [p1[0]+s1*d[0], p1[1]+s1*d[1], p1[2]+s1*d[2]]
            s2 = s1 - self.f(p)/dot(d, self.normal(p))
            delta = abs(s2 - s1)
            s1 = s2
        return s0+s1, p


class Spherical(SurfaceProfile):
    def __init__(self, c=0.0):
        self.type = 'Sphere'
        self.cv = c

    def normal(self, p):
        return [-self.cv*p[0], -self.cv*p[1], 1.0-self.cv*p[2]]

    def f(self, p):
        return p[2] - 0.5*self.cv*(p[0]*p[0]+p[1]*p[1])


class Conic(SurfaceProfile):
    def __init__(self, c=0.0):
        self.type = 'Conic'
        self.cv = c
        self.cc = 0.0


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

    dir0 = [0., 0., 1.]
    v1 = normalize([1., 1., 1.])
    p0 = [0., 0., -1.]
    p1 = [0., 1., -1.]

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
