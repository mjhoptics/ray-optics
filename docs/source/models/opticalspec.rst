.. currentmodule:: rayoptics.raytr

.. |Wvl| replace:: :class:`~opticalspec.WvlSpec`
.. |Pupil| replace:: :class:`~opticalspec.PupilSpec`
.. |FieldSpec| replace:: :class:`~opticalspec.FieldSpec`
.. |Field| replace:: :class:`~opticalspec.Field`
.. |FocusRange| replace:: :class:`~opticalspec.FocusRange`

***************************
Optical Usage Specification
***************************

Overview
========

   The :class:`~opticalspec.OpticalSpecs` class holds the optical usage definition of the model. Aperture, field of view, wavelength, and focal position are all aspects of the :class:`~opticalspec.OpticalSpecs`.

   The optical configuration is broken into four parts::

      * aperture
      * field of view
      * spectral region
      * focus range

   The pupil and field specifications can be specified in a variety of ways. The ``key`` keyword argument takes a list of 2 strings. The first string indicates whether the specification is in object or image space. The second one indicates which parameter is the defining specification.

Pupil Specification
-------------------

   The |Pupil| class maintains the aperture specification. The `PupilSpec` can be defined in object or image space. The defining parameters can be ``epd``, ``f/#`` or ``NA``, where ``epd`` is the pupil diameter.

   .. code:: ipython3

       osp.pupil = PupilSpec(osp, key=['object', 'epd'], value=12.5)

   The |Pupil| class allows rays to be specified as fractions of the pupil dimension. A list of pupil_rays and ray_labels define rays to be used to establish clear aperture dimensions on optical elements and rays to be drawn for the lens layout. A default set of pupil rays is provided that is appropriate for circular pupil systems with plane symmetry.

Field Specification
-------------------

   The |FieldSpec| can be defined in object or image space. The defining parameters can be ``height`` or ``angle``, where ``angle`` is given in degrees.

   .. code:: ipython3

       osp.field_of_view = FieldSpec(osp, key=['object', 'angle'], flds=[0., 20.0])

   The |FieldSpec| maintains a list of |Field| instances. Each |Field| contains an absolute field specification of the type specified in |FieldSpec|. It can also have attributes set for the chief ray data at the field (:attr:`~opticalspec.Field.chief_ray`) as well as the definition of the reference sphere for OPD calculations (:attr:`~opticalspec.Field.ref_sphere`). These are calculated by :func:`trace.get_chief_ray_pkg` and :func:`trace.setup_pupil_coords` respectively.

Wavelength Specification
------------------------

   The |Wvl| defines the wavelengths and weights to use when evaluating the model. The wavelength values can be given in either nanometers or a spectral line designation.

   .. code:: ipython3

       osp.spectral_region = WvlSpec([('F', 0.5), (587.5618, 1.0), ('C', 0.5)], ref_wl=1)

   The :class:`~.opticalspec.FocusRange` defines the amount of defocus and, optionally, focal range to use when evaluating the model.

   .. code:: ipython3

       osp.defocus = FocusRange(focus_shift=0.01)

