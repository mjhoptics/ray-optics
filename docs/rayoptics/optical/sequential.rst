.. currentmodule:: rayoptics.optical

.. _seq_model:

****************
Sequential Model
****************

Overview
========

    The Sequential model abstracts the essential parts of a physical optical system when using ray tracing to compute imaging performance. For many use cases, this only requires a minimum of information (curvature, thickness, refractive index) to usefully describe an optical system. Similarly, the needs of ray tracing for performance evaluation is satisfied by having the sequence of optical surfaces involved in image formation.
    
    A sequential optical model is a sequence of interfaces and gaps. The first interface is the object surface and the last is the image surface. Light propagates sequentially throught the interfaaces and gaps.

    The sequential model has this structure
    ::

        IfcObj  Ifc1  Ifc2  Ifc3 ... Ifci-1   IfcImg
             \  /  \  /  \  /             \   /
             GObj   G1    G2              Gi-1

    where

        - Ifc is a :class:`~surface.Interface` instance
        - G   is a :class:`~gap.Gap` instance

    There are N interfaces and N-1 gaps. The initial configuration has an object and image :class:`surface.Surface` and an object gap.

    The :class:`surface.Interface` API supports implementation of an optical action, such as refraction, reflection, scatter, diffraction, etc. The :class:`~surface.Interface` may be realized as a physical profile separating the adjacent gaps or an idealized object, such as a thin lens or 2 point HOE.

    The :class:`gap.Gap` class maintains a simple separation (z translation) and the medium filling the gap. More complex coordinate transformations are handled through the Interface API.