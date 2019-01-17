.. currentmodule:: rayoptics.optical

.. _opt_model:

*************
Optical Model
*************

Overview
========

   The :class:`~opticalmodel.OpticalModel` serves as a top level container of model properties. Key aspects are built-in element and surface based repesentations of the optical surfaces.

   The :class:`opticalmodel.OpticalModel` contains several different representations of aspects of an optical system. There is the physical model of the components, contained in the :class:`elements.ElementModel` class. A major requirement is the elements be able to render themselves in a 2D drawing view.

   A second representation is the :class:`sequential.SequentialModel`, that describes a sequence of Interfaces and Gaps that light propagates through to form an image. The basic purpose of this representation is for geometric ray tracing. The :func:`raytrace.trace_raw` function is the building block of the ray trace system.

   The :class:`opticalspec.OpticalSpecs` class holds the optical usage definition of the model. Aperture, field of view, wavelength, and focal position are all aspects of the :class:`~opticalspec.OpticalSpecs`. The first order properties are calculated and maintained by :class:`~opticalspec.OpticalSpecs` in the ``parax_data`` variable.

   Paraxial layout operations are handled via the :class:`paraxialdesign.ParaxialModel`. This enables model manipulation via the :math:`y-\overline{y}` or :math:`\omega-\overline{\omega}` diagrams.

   Finally, system units and other information is maintained in the :class:`opticalmodel.SystemSpec` class.
