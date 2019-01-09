.. currentmodule:: rayoptics.optical

***************************
Optical Usage Specification
***************************

.. _optical_spec:

Optical Specification Overview
==============================

   The :class:`opticalspec.OpticalSpecs` class holds the optical usage definition of the model. Aperture, field of view, wavelength, and focal position are all aspects of the :class:`~opticalspec.OpticalSpecs`.

   The first order properties are calculated and maintained by :class:`~opticalspec.OpticalSpecs` in the ``parax_data`` variable. These include the paraxial axial and chief rays, and the :class:`~firstorder.FirstOrderData` that contains first order properties.
