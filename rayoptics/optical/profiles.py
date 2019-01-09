#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2017 Michael J. Hayford
""" Module for different surface profile shapes

.. Created on Tue Aug  1 13:18:57 2017

.. codeauthor: Michael J. Hayford
"""

import numpy as np
from math import sqrt, copysign, sin, atan2
from rayoptics.util.misc_math import normalize
from .traceerror import TraceMissedSurfaceError


def resize_list(lst, new_length, null_item=None):
    return lst + [null_item for item in range(new_length - len(lst))]


def intersect_parabola(cv, p, d, z_dir=1.0):
    ''' Intersect a parabolid, starting from an arbitrary point.

    Args:
        cv: vertex curvature
        p:  start point of the ray in the profile's coordinate system
        d:  direction cosine of the ray in the profile's coordinate system
        z_dir: 
    '''
    # Intersection with a conic, starting from an arbitrary point.
    #
    # For quadratic equation ax**2 + bx + c = 0:
    #  ax2 = 2a
    #  cx2 = 2c
    ax2 = cv*(1. - d[2]*d[2])
    cx2 = cv*(p[0]*p[0] + p[1]*p[1]) - 2.0*p[2]
    b = cv*(d[0]*p[0] + d[1]*p[1]) - d[2]
    # Use z_dir to pick correct root
    s = cx2/(z_dir*sqrt(b*b - ax2*cx2) - b)

    p1 = p + s*d
    return s, p1


class SurfaceProfile:
    def __repr__(self):
        return "{!s}()".format(type(self).__name__)

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

    def intersect(self, p0, d, eps, z_dir):
        ''' Intersect a profile, starting from an arbitrary point.

        Args:
            cv: vertex curvature
            p0:  start point of the ray in the profile's coordinate system
            d:  direction cosine of the ray in the profile's coordinate system
            z_dir: +1 if propagation positive direction, -1 if otherwise
            eps: numeric tolerance for convergence of any iterative procedure

        Returns:
            s1: distance to intersection point
            p: intersection point
        '''
        p = p1 = p0
        s1 = -self.f(p1)/np.dot(d, self.normal(p1))
        delta = abs(s1)
#        print("intersect", s1)
        while delta > eps:
            p = p1 + s1*d
            s2 = s1 - self.f(p)/np.dot(d, self.normal(p))
            delta = abs(s2 - s1)
#            print("intersect", s1, s2, delta)
            s1 = s2
        return s1, p


class Spherical(SurfaceProfile):
    """ Spherical surface profile parameterized by curvature. """
    def __init__(self, c=0.0, r=None):
        """ initialize a Spherical profile.

        Args:
            c: curvature
            r: radius of curvature. If zero, taken as planar. If r is
                specified, it overrides any input for c (curvature).
        """
        if r:
            self.r = r
        else:
            self.cv = c

    @property
    def r(self):
        if self.cv != 0.0:
            return 1.0/self.cv
        else:
            return 0.0

    @r.setter
    def r(self, radius):
        if radius != 0.0:
            self.cv = 1.0/radius
        else:
            self.cv = 0.0

    def __str__(self):
        return type(self).__name__ + " " + str(self.cv)

    def __repr__(self):
        return "{!s}(c={})".format(type(self).__name__, self.cv)

    def copyFrom(self, other):
        dispatch[type(self), type(other)](self, other)

    def copyDataFrom(self, other):
        self.cv = other.cv

    def mutate(self, new_profile):
        Spherical.copyDataFrom(new_profile, self)

    def normal(self, p):
        return normalize(np.array(
                [-self.cv*p[0], -self.cv*p[1], 1.0-self.cv*p[2]]))

    def intersect_tangent_plane(self, p, d, eps, z_dir):
        # Welford's intersection with a sphere, starting from the tangent plane
        # transfer p to tangent plane of surface
        s0 = -p[2]/d[2]
        p0 = np.array([p[0] + s0*d[0], p[1] + s0*d[1], 0.])

        # Welford's 4.8, 4.9 and 4.12
        F = self.cv*(p0[0]*p0[0] + p0[1]*p0[1])
        G = d[2] - self.cv*(d[0]*p0[0] + d[1]*p0[1])
        s1 = F/(G + sqrt(G*G - F*self.cv))

        s = s0 + s1

        p1 = p + s*d
        return s, p1

    def intersect(self, p, d, eps, z_dir):
        ''' Intersection with a sphere, starting from an arbitrary point. '''

        # Substitute expressions equivalent to Welford's 4.8 and 4.9
        # For quadratic equation ax**2 + bx + c = 0:
        #  ax2 = 2a
        #  cx2 = 2c
        ax2 = self.cv
        cx2 = self.cv * p.dot(p) - 2.0*p[2]
        b = 2.0*self.cv * d.dot(p) - d[2]
        # Use z_dir to pick correct root
        s = cx2/(z_dir*sqrt(b*b - ax2*cx2) - b)

        p1 = p + s*d
        return s, p1

    def f(self, p):
        return p[2] - 0.5*self.cv*(np.dot(p, p))

    def sag(self, x, y):
        if self.cv != 0.0:
            r = 1/self.cv
            try:
                adj = sqrt(r*r - x*x - y*y)
            except ValueError:
                raise TraceMissedSurfaceError(self, (x, y))
            finally:
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
    def __init__(self, c=0.0, cc=0.0, r=None, ec=None):
        """ initialize a Conic profile.

        Args:
            c: curvature
            r: radius of curvature. If zero, taken as planar. If r is
                specified, it overrides any input for c (curvature).
            cc: conic constant
            ec: conic asphere (= cc + 1). If ec is specified, it overrides any
                input for the conic constant (cc).
        """
        if r:
            self.r = r
        else:
            self.cv = c

        if ec:
            self.ec = ec
        else:
            self.cc = cc

    @property
    def r(self):
        if self.cv != 0.0:
            return 1.0/self.cv
        else:
            return 0.0

    @r.setter
    def r(self, radius):
        if radius != 0.0:
            self.cv = 1.0/radius
        else:
            self.cv = 0.0

    @property
    def ec(self):
        return self.cc + 1.0

    @ec.setter
    def ec(self, ec):
        self.cc = ec - 1.0

    def __str__(self):
        return type(self).__name__ + " " + str(self.cv) + " " + \
                                           str(self.cc)

    def __repr__(self):
        return "{!s}(c={}, cc={})".format(type(self).__name__,
                                          self.cv, self.cc)

    def copyFrom(self, other):
        dispatch[type(self), type(other)](self, other)

    def copyDataFrom(self, other):
        self.cv = other.cv
        self.cc = other.cc

    def normal(self, p):
        return normalize(np.array(
                [-self.cv*p[0],
                 -self.cv*p[1],
                 1.0-(self.cc+1.0)*self.cv*p[2]]))

    def intersect_tangent_plane(self, p, d, eps, z_dir):
        # Welford's intersection with a conic, starting from the tangent plane
        # transfer p to tangent plane of surface
        s0 = -p[2]/d[2]
        p0 = np.array([p[0] + s0*d[0], p[1] + s0*d[1], 0.])

        # Welford's 4.8, 4.9 and 4.28
        ax2 = self.cv*(1. + self.cc*d[2]*d[2])
        F = self.cv*(p0[0]*p0[0] + p0[1]*p0[1])
        G = d[2] - self.cv*(d[0]*p0[0] + d[1]*p0[1])
        s1 = F/(G + z_dir*sqrt(G*G - F*ax2))

        s = s0 + s1
        p1 = p + s*d
        return s, p1

    def intersect(self, p, d, eps, z_dir):
        ''' Intersection with a conic, starting from an arbitrary point.'''

        # For quadratic equation ax**2 + bx + c = 0:
        #  ax2 = 2a
        #  cx2 = 2c
        ax2 = self.cv*(1. + self.cc*d[2]*d[2])
        cx2 = self.cv*(p[0]*p[0] + p[1]*p[1] + self.ec*p[2]*p[2]) - 2.0*p[2]
        b = self.cv*(d[0]*p[0] + d[1]*p[1] + self.ec*d[2]*p[2]) - d[2]
        # Use z_dir to pick correct root
        s = cx2/(z_dir*sqrt(b*b - ax2*cx2) - b)

        p1 = p + s*d
        return s, p1

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
    def __init__(self, c=0.0, cc=0.0, r=None, ec=None, coefs=None):
        """ initialize a EvenPolynomial profile.

        Args:
            c: curvature
            r: radius of curvature. If zero, taken as planar. If r is
                specified, it overrides any input for c (curvature).
            cc: conic constant
            ec: conic asphere (= cc + 1). If ec is specified, it overrides any
                input for the conic constant (cc).
            coefs: a list of even power coefficents, starting with the
                quadratic term, and not exceeding the 20th order term.
        """
        if r:
            self.r = r
        else:
            self.cv = c

        if ec:
            self.ec = ec
        else:
            self.cc = cc

        if coefs:
            self.coefs = coefs
        else:
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

    @property
    def r(self):
        if self.cv != 0.0:
            return 1.0/self.cv
        else:
            return 0.0

    @r.setter
    def r(self, radius):
        if radius != 0.0:
            self.cv = 1.0/radius
        else:
            self.cv = 0.0

    @property
    def ec(self):
        return self.cc + 1.0

    @ec.setter
    def ec(self, ec):
        self.cc = ec - 1.0

    def __str__(self):
        return type(self).__name__ + " " + str(self.cv) + " " + \
                                           str(self.cc)

    def __repr__(self):
        return "Profile({}: c={}, cc={}".format(type(self).__name__,
                                                self.cv, self.cc)

    def copyFrom(self, other):
        dispatch[type(self), type(other)](self, other)

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
        if len(self.coefs) == 0:
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


class RadialPolynomial(SurfaceProfile):
    """ Radial Polynomial asphere, on base conic.

    Conics produced for conic asphere values:
        ec > 1.0: oblate spheroid
        ec = 1.0: sphere
        ec > 0.0 and < 1.0: ellipsoid
        ec = 0.0: paraboloid
        ec < 0.0: hyperboloid
    """
    initial_size = 10

    def __init__(self, c=0.0, cc=None, r=None, ec=1.0, coefs=None):
        """ initialize a RadialPolynomial profile.

        Args:
            c: curvature
            r: radius of curvature. If zero, taken as planar. If r is
                specified, it overrides any input for c (curvature).
            ec: conic asphere.
            cc: conic constant (= ec - 1). If cc is specified, it overrides any
                input for the conic asphere (ec).
            coefs: a list of radial coefficents, starting with the
                constant term, (and not exceeding the 10th order term).
        """
        if r:
            self.r = r
        else:
            self.cv = c

        if cc:
            self.cc = cc
        else:
            self.ec = ec

        if coefs:
            self.coefs = coefs
        else:
            self.coefs = [0. for i in range(self.initial_size)]
        self.max_nonzero_coef = 0
        self.coef1 = 0.0
        self.coef2 = 0.0
        self.coef3 = 0.0
        self.coef4 = 0.0
        self.coef5 = 0.0
        self.coef6 = 0.0
        self.coef7 = 0.0
        self.coef8 = 0.0
        self.coef9 = 0.0
        self.coef10 = 0.0

    @property
    def r(self):
        if self.cv != 0.0:
            return 1.0/self.cv
        else:
            return 0.0

    @r.setter
    def r(self, radius):
        if radius != 0.0:
            self.cv = 1.0/radius
        else:
            self.cv = 0.0

    @property
    def cc(self):
        return self.ec - 1.0

    @cc.setter
    def cc(self, cc):
        self.ec = cc + 1.0

#    @property
    def get_coef(self, exp):
        return self.coefs[exp]

#    @coef.setter
    def set_coef(self, exp, value):
        try:
            self.coefs[exp] = value
        except IndexError:
            self.coefs = resize_list(self.coefs, exp, null_item=0.0)
            self.coefs[exp] = value

    def __str__(self):
        return type(self).__name__ + " " + str(self.cv) + " " + \
                                           str(self.ec)

    def __repr__(self):
        return "Profile({}: c={}, ec={}".format(type(self).__name__,
                                                self.cv, self.ec)

    def copyFrom(self, other):
        dispatch[type(self), type(other)](self, other)

    def copyDataFrom(self, other):
        self.cv = other.cv
        self.ec = other.ec
        self.coefs = other.coefs.copy()
        self.coefs = other.coefs
        self.coef1 = other.coef1
        self.coef2 = other.coef2
        self.coef3 = other.coef3
        self.coef4 = other.coef4
        self.coef5 = other.coef5
        self.coef6 = other.coef6
        self.coef7 = other.coef7
        self.coef8 = other.coef8
        self.coef9 = other.coef9
        self.coef10 = other.coef10

    def gen_coef_list(self):
        if len(self.coefs) == 0:
            self.coefs = []
            self.coefs.append(self.coef1)
            self.coefs.append(self.coef2)
            self.coefs.append(self.coef3)
            self.coefs.append(self.coef4)
            self.coefs.append(self.coef5)
            self.coefs.append(self.coef6)
            self.coefs.append(self.coef7)
            self.coefs.append(self.coef8)
            self.coefs.append(self.coef9)
            self.coefs.append(self.coef10)
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
        r = sqrt(r2)
        e = self.cv/sqrt(1. - self.ec*self.cv*self.cv*r2)

        # polynomial asphere contribution - compute using Horner's Rule
        e_asp = 0.0
        if r == 0.0:
            r_pow = 1.0
        else:
            # Initialize to 1/r because we multiply by r's components p[0] and
            # p[1] at the final normalization step.
            r_pow = 1.0/r
        c_coef = 1.0
        for coef in self.coefs[1:self.max_nonzero_coef]:
            e_asp += c_coef*coef*r_pow
            c_coef += 1.0
            r_pow *= r

        e_tot = e + e_asp
        return normalize(np.array([-e_tot*p[0], -e_tot*p[1], 1.0]))

    def sag(self, x, y):
        r2 = x*x + y*y
        r = sqrt(r2)
        # sphere + conic contribution
        z = self.cv*r2/(1. + sqrt(1. - self.ec*self.cv*self.cv*r2))

        # polynomial asphere contribution
        z_asp = 0.0
        r_pow = 1.0
        for coef in self.coefs[:self.max_nonzero_coef]:
            z_asp += coef*r_pow
            r_pow *= r

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
  (Spherical, RadialPolynomial): Spherical.copyDataFrom,
  (Conic, Spherical): Spherical.copyDataFrom,
  (Conic, Conic): Conic.copyDataFrom,
  (Conic, EvenPolynomial): Conic.copyDataFrom,
  (Conic, RadialPolynomial): Conic.copyDataFrom,
  (EvenPolynomial, Spherical): Spherical.copyDataFrom,
  (EvenPolynomial, Conic): Conic.copyDataFrom,
  (EvenPolynomial, EvenPolynomial): EvenPolynomial.copyDataFrom,
  (EvenPolynomial, RadialPolynomial): EvenPolynomial.copyDataFrom,
  (RadialPolynomial, Spherical): Spherical.copyDataFrom,
  (RadialPolynomial, Conic): Conic.copyDataFrom,
  (RadialPolynomial, EvenPolynomial): EvenPolynomial.copyDataFrom,
  (RadialPolynomial, RadialPolynomial): RadialPolynomial.copyDataFrom,
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
