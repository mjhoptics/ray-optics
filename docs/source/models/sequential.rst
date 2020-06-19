.. currentmodule:: rayoptics.optical
.. |opm| replace:: :class:`~opticalmodel.OpticalModel`
.. |sm| replace:: :class:`~sequential.SequentialModel`
.. |ifc| replace:: :class:`~interface.Interface`
.. |surf| replace:: :class:`~surface.Surface`
.. |gap| replace:: :class:`~gap.Gap`

****************
Sequential Model
****************

Overview
========

    The |sm| abstracts the essential parts of a physical optical system when using ray tracing to compute imaging performance. For many use cases, this only requires a minimum of information (curvature, thickness, refractive index) to usefully describe an optical system. Similarly, the needs of ray tracing for performance evaluation is satisfied by having the sequence of optical surfaces involved in image formation.
    
    A sequential optical model is a sequence of interfaces and gaps. The first interface is the object surface and the last is the image surface. Light propagates sequentially throught the interfaaces and gaps.

    The sequential model has this structure
    ::

        IfcObj  Ifc1  Ifc2  Ifc3 ... Ifci-1   IfcImg
             \  /  \  /  \  /             \   /
             GObj   G1    G2              Gi-1

    where

        - Ifc is a |ifc| instance
        - G   is a |gap| instance

    There are N interfaces and N-1 gaps. The initial configuration has an object and image |surf| and an object gap.

    The |ifc| API supports implementation of an optical action, such as refraction, reflection, scatter, diffraction, etc. The |ifc| may be realized as a physical profile separating the adjacent gaps or an idealized object, such as a thin lens or 2 point HOE.

    The |gap| class maintains a simple separation (z translation) and the medium filling the gap. More complex coordinate transformations are handled through the Interface API. The medium in the |gap| is provided via classes that implement :meth:`~medium.Medium.rindex`. For simple modeling, the :mod:`medium` provides constant refractive index classes for :class:`~medium.Air` and :class:`~medium.Glass`. Commercial optical glasses are supported by the `opticalglass <https://opticalglass.readthedocs.io>`_ package.

    The combination of Interfaces and Gaps completely define the path of a ray through the optical model. For convenience and efficiency when doing optical calculations, it is useful to extract a subset of data that is required for ray tracing. The method :meth:`~sequential.SequentialModel.path` in |sm| returns a Python generator that is used to drive the ray tracing process.
