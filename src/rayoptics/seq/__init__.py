""" Package for sequential modeling of optical systems

    The :mod:`~.seq` subpackage provides core classes and functions for
    sequential optical modeling. A sequential optical model is a sequence of
    surfaces and gaps.

    The sequential model has this structure
    ::

        IfcObj  Ifc1  Ifc2  Ifc3 ... Ifci-1   IfcImg
             \  /  \  /  \  /             \   /
             GObj   G1    G2              Gi-1

    where

        - Ifc is a :class:`~rayoptics.seq.interface.Interface` instance
        - G   is a :class:`~rayoptics.seq.gap.Gap` instance

    There are N interfaces and N-1 gaps. The initial configuration has an
    object and image Surface and an object gap.

    The Interface API supports implementation of an optical action, such as
    refraction, reflection, scatter, diffraction, etc. The Interface may be
    realized as a physical profile separating the adjacent gaps or an idealized
    object, such as a thin lens or 2 point HOE.

    The Gap class maintains a simple separation (z translation) and the medium
    filling the gap. More complex coordinate transformations are handled
    through the Interface API.

    Modules comprising the sequential optical modeling capability include:

        - Sequential optical model: :mod:`~.sequential`
        - Interface and gap modules:
          :mod:`~.interface`
          :mod:`~.gap`
          :mod:`~.medium`
        - Modules for specialized optical behavior:
          :mod:`~.twoconicmirrors`

"""
