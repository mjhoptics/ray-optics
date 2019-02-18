.. currentmodule:: rayoptics.optical

.. |Wvl| replace:: :class:`~opticalspec.WvlSpec`
.. |Pupil| replace:: :class:`~opticalspec.PupilSpec`
.. |FieldSpec| replace:: :class:`~opticalspec.FieldSpec`
.. |Field| replace:: :class:`~opticalspec.Field`
.. |FocusRange| replace:: :class:`~opticalspec.FocusRange`

***************************
Optical Usage Specification
***************************

.. _optical_spec:

Optical Specification Overview
==============================

   The :class:`opticalspec.OpticalSpecs` class holds the optical usage definition of the model. Aperture, field of view, wavelength, and focal position are all aspects of the :class:`~opticalspec.OpticalSpecs`.

   The first order properties are calculated and maintained by :class:`~opticalspec.OpticalSpecs` in the ``parax_data`` variable. These include the paraxial axial and chief rays, and the :class:`~firstorder.FirstOrderData` that contains first order properties.

   The optical configuration is broken into four parts::

      * aperture
      * field of view
      * spectral region
      * focus range

   The |Pupil| class maintains the aperture specification. The type of specification is given by the :class:`model_enums.PupilType` enum. Four types of specification are provided::

      * FNO: image space f/#
      * EPD: entrance pupil diameter
      * NA: image space numerical aperture
      * NAO: object space numerical aperture

   The |Pupil| class allows rays to be specified as fractions of the pupil dimension. A list of pupil_rays and ray_labels define rays to be used to establish clear aperture dimensions on optical elements and rays to be drawn for the lens layout. A default set of pupil rays is provided that is appropriate for circular pupil systems with plane symmetry.

   The |FieldSpec| maintains a list of |Field| instances. Each |Field| contains an absolute field specification of the type specified in |FieldSpec|. It can also have attributes set for the chief ray data at the field (:attr:`~opticalspec.Field.chief_ray`) as well as the definition of the reference sphere for OPD calculations (:attr:`~opticalspec.Field.ref_sphere`). These are calculated by :func:`trace.setup_pupil_coords`.
