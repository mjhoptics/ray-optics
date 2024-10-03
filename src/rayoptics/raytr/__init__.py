""" Package for optical ray tracing and calculations

    The :mod:`~.raytr` subpackage provides core classes and functions
    for optical ray tracing and analyses. These include:

        - Base level ray tracing, :mod:`~.raytrace`
        - Calculation of wavefront aberration, :mod:`~.waveabr`
        - Specification of aperture, field, wavelength and defocus,
          :mod:`~.opticalspec`
        - Higher level ray tracing, in terms of aperture, field and wavelength,
          :mod:`~.trace`
        - Functions setting vignetting and clear apertures and support for 
          pupil exploration, :mod:`~.vigcalc`
        - Tracing of fans, lists and grids of rays, including refocusing of OPD
          values, :mod:`~.analyses`
        - Exception classes for reporting ray trace errors, :mod:`~.traceerror`
        - Sample generation for ray grids, :mod:`~.sampler`

    The overall optical model is managed by the :class:`~.OpticalModel` class
"""

from collections import namedtuple

RayResult = namedtuple('RayResult', ['pkg', 'err'])
RayResult.__doc__ = "A RayPkg and either a RayError or None"
RayResult.pkg.__doc__ = "a RayPkg"
RayResult.err.__doc__ = "a RayError or None, if success"

RayPkg = namedtuple('RayPkg', ['ray', 'op', 'wvl'])
RayPkg.__doc__ = "Ray and optical path length, plus wavelength"
RayPkg.ray.__doc__ = "list of RaySegs"
RayPkg.op.__doc__ = "optical path length between pupils"
RayPkg.wvl.__doc__ = "wavelength (in nm) that the ray was traced in"

RaySeg = namedtuple('RaySeg', ['p', 'd', 'dst', 'nrml'])
RaySeg.__doc__ = "ray intersection and transfer data"
RaySeg.p.__doc__ = "the point of incidence"
RaySeg.d.__doc__ = "ray direction cosine following the interface"
RaySeg.dst.__doc__ = "geometric distance to next point of incidence"
RaySeg.nrml.__doc__ = "surface normal vector at the point of incidence"
