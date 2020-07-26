""" Package for paraxial optical design and analysis

    The :mod:`~.parax` subpackage provides core classes and functions
    for paraxial/first order optical design. These include:

        - Definition of the first order quantities defining an optical system,
          :mod:`~.specsheet`, :mod:`~.idealimager`, and :mod:`~.etendue`
        - Paraxial and third order calculations, :mod:`~.firstorder`,
          :mod:`~.thirdorder`
        - Support for |ybar| and |nubar| diagrams, :mod:`~.paraxialdesign` and
          :mod:`~.diagram`

    The paraxial model is managed by the :class:`~.ParaxialModel` class
"""
