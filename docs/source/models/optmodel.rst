.. currentmodule:: rayoptics
.. |opm| replace:: :class:`~optical.opticalmodel.OpticalModel`
.. |sm| replace:: :class:`~seq.sequential.SequentialModel`
.. |osp| replace:: :class:`~raytr.opticalspec.OpticalSpecs`
.. |specsheet| replace:: :class:`~parax.specsheet.SpecSheet`
.. |pm| replace:: :class:`~parax.paraxialdesign.ParaxialModel`
.. |em| replace:: :class:`~elem.elements.ElementModel`

*************
Optical Model
*************

Overview
========

   The |opm| serves as a top level container of model properties. Key aspects are built-in sequential surface, paraxial and element based repesentations of the optical model. Each representation can be used build and modify the optical model. 

   For design and analysis of image forming systems, the sequential surface model is often the most direct method for specifying an optical system. The |sm| describes a sequence of Interfaces and Gaps that light propagates through to form an image. First order properties and ray tracing can be done directly with this representation. It is straightforward to produce an element based representation from the |sm| in many cases. Changes in the |sm| can be propagated to the peer models in the |opm| by calling :meth:`~optical.opticalmodel.OpticalModel.update_model` on |opm|. 

   A |sm| is a high fidelity representation of the geometry and optical properties of the optical model. In contrast, the |pm| is obtained by assuming all optical ray angles are small so that trigonometric relationships can be replaced by the angles themselves. This enables deriving the constructional parameters from paraxial ray specifications. This can be done using paraxial ray solves and/or |ybar| and |nubar| diagrams.

   Finally, the physical model of the components is managed in the |em| class. A major requirement is the elements be able to render themselves in a 2D drawing view. It is planned to use the grouping of interfaces into elements to allow |ybar| diagram to collapse individual nodes for interfaces into a single node for elements or groups of elements. It should also be possible in the future to ray trace the |em| to generate a |sm| corresponding to a ray's path through the elements.

   The |osp| class holds the optical usage definition of the model. Aperture, field of view, wavelength, and focal position are all aspects of the |osp|. The first order properties are calculated and maintained by |osp| in the :attr:`~raytr.opticalspec.OpticalSpecs.parax_data` variable.  The usage specification of the optical model is given by the |specsheet| class in combination with the |osp| class. 

   Paraxial layout operations are handled via the |pm|. This enables model manipulation via the |ybar| or |nubar| diagrams. The goal is to smoothly integrate the information in the |specsheet| with constraints on movement in the |ybar| diagram or changes to the paraxial rays in the Optical Layout.

   The system units and other information is contained in the :class:`~optical.opticalmodel.SystemSpec` class. It is the  :attr:`~optical.opticalmodel.OpticalModel.system_spec` attribute in |opm|.
