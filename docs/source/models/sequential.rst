.. currentmodule:: rayoptics


****************
Sequential Model
****************

Overview
========

    The :class:`~seq.sequential.SequentialModel` abstracts the essential parts of a physical optical system when using ray tracing to compute imaging performance. For rotationally symmetric optical systems, this requires a minimum of information: a list of curvature, thickness, and refractive index. Similarly, the needs of ray tracing for performance evaluation is satisfied by the sequence of optical surfaces involved in image formation.
    
    A sequential optical model is a sequence of interfaces and gaps. The first interface is the object surface and the last is the image surface. Light propagates sequentially through the interfaces and gaps.

Structure
=========

    The sequential model has this structure::

        IfcObj  Ifc1  Ifc2  Ifc3 ... IfcN-1   IfcN  IfcImg
             \  /  \  /  \  /             \   /  \  /
             GObj   G1    G2              GN-1   GImg

    where

        - Ifc is a :class:`~seq.interface.Interface` instance
        - G   is a :class:`~seq.gap.Gap` instance

    A complete :class:`~seq.sequential.SequentialModel` mates a top level assembly with an Object space and Image space. An advantage of this structure is that it facilitates the substitution of one set of interfaces for another. A top level view of a :class:`~seq.sequential.SequentialModel` looks like this::

        IfcObj  System  IfcImg
             \  /    \  /  
             GObj    GImg

    The Object space consists of an object interface and object gap, with Image space similarly constituted. Note that the Image gap is associated with the Image surface, not the final System interface. This makes editing operations such as inserting or replacing elements straightforward. The user view (e.g. :meth:`~seq.sequential.SequentialModel.list_model`) still groups the image thickness with the final interface, with no loss of generality.

    Sub-assemblies are defined by N interfaces and N-1 gaps. The low level view of System looks like this::

        Ifc1  Ifc2  Ifc3 ... IfcN-1   IfcN
           \  /  \  /             \  /
            G1    G2              GN-1

    The low level view might be easier to manipulate at an intermediate level of detail. For example, a Telescope might be represented as::

        IfcObj  Telescope  IfcImg
             \  /       \  /  
             GObj       GImg

    where the Telescope is::

        Objective  Eyepiece
                \  /  
                GSep

    and GSep is the gap separating the Objective and the Eyepiece.

Constituents
============

    The :class:`~seq.interface.Interface` API supports implementation of an optical action, such as refraction, reflection, scatter, diffraction, etc. The :class:`~seq.interface.Interface` may be realized as a physical profile separating the adjacent gaps or an idealized object, such as a thin lens or 2 point HOE.

    The :class:`~seq.gap.Gap` class maintains a simple separation (z translation) and the medium filling the gap. More complex coordinate transformations are handled through the Interface API.

    The medium in the :class:`~seq.gap.Gap` is provided via classes that implement :meth:`~.medium.Medium.rindex`. For simple modeling, the :mod:`~.medium` provides constant refractive index classes for :class:`~seq.medium.Air` and :class:`~seq.medium.Glass`. 

    If you have a list of refractive indices and their corresponding wavelengths (in nm), you can use the :class:`~seq.medium.InterpolatedGlass` class. This uses the :mod:`scipy` routine :class:`scipy.interpolate.interp1d` to interpolate the data.

    Commercial optical glasses are supported by the `opticalglass <https://opticalglass.readthedocs.io>`_ package.

Sequential Paths
================

    The combination of Interfaces and Gaps completely define the path of a ray through the optical model. For convenience and efficiency when doing optical calculations, it is useful to extract a subset of data that is required for ray tracing. The method :meth:`~seq.sequential.SequentialModel.path` in :class:`~seq.sequential.SequentialModel` returns a Python generator that is used to drive the ray tracing process.
