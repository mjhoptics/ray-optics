""" Package for optical ray tracing and calculations

    The :mod:`~.raytr` subpackage provides core classes and functions
    for optical ray tracing and analyses. These include:

        - Primitive and higher level ray tracing, :mod:`~.raytrace`,
          :mod:`~.trace`
        - Specification of aperture, field, wavelength and defocus,
          :mod:`~.opticalspec`
        - Tracing of fans, lists and grids of rays, including refocusing of OPD
          values, :mod:`~.analyses`
        - Sample generation for ray grids, :mod:`~.sampler`

    The overall optical model is managed by the :class:`~.OpticalModel` class
"""
