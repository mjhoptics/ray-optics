""" Package for optical ray tracing and calculations

    The :mod:`~.raytr` subpackage provides core classes and functions
    for optical ray tracing and analyses. These include:

        - Primitive and higher level ray tracing, :mod:`~.raytrace`,
          :mod:`~.trace`
        - Specification of aperture, field, wavelength and defocus,
          :mod:`~.opticalspec`
        - Calculation of wavefront aberration, :mod:`~.waveabr`
        - Tracing of fans, lists and grids of rays, including refocusing of OPD
          values, :mod:`~.analyses`
        - Sample generation for ray grids, :mod:`~.sampler`

    The overall optical model is managed by the :class:`~.OpticalModel` class
"""

from collections import namedtuple

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
