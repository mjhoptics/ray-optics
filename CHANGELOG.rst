.. currentmodule:: rayoptics

=========
Changelog
=========


Version 0.9.3
=============
In response to issue 157, a wide angle raytrace mode was added to the :class:`~.raytr.opticalspec.FieldSpec` model. This can handle lenses over 90 deg half angle, and performs well for fields>45 deg. Chief ray aiming is done by locating the real entrance pupil location in object space rather than aiming at the paraxial entrance pupil, as in the default field specification. Investigation led to the realization that many wide angle models were specified in terms of real image height, which can be vastly different from the paraxial value. The chief ray is traced in reverse from the image point into object space. Both the real entrance pupil and the real image height calculations are in the new :mod:`~.wideangle` module. A different approach towards calculating vignetting based on a binary search to find the limiting aperture was added. This performs better than the default for wide angle systems, where ray failures during the search need to be handled robustly. Finally, the entire ray trace stack from the :mod:`~.analyses` module to :func:`~.raytr.trace.trace_ray` to :func:`~.raytr.raytrace.trace_raw` was refactored to uniformly pass kwargs, and provide ray trace exception handling. The ability to display rays up to the point they were clipped is helpful for model debugging. The code will continue to be streamlined and refactored as needed to accomodate new functionality and get rid of obsolete ideas.

There were improvements to the importers, both to use the wide angle support and real image height capabilities, and routine maintenance.

Thanks to @Tohuvavohu for correcting an error in the :mod:`~.zmxread` module.


Version 0.9.2
=============
Revised the call list for :func:`~.raytr.trace.trace_ray` to make it straightforward to use in the simplest cases. The function returns a tuple of the ray_pkg and an ray error instance; the latter is None if are no errors. The :func:`~.raytr.trace.list_ray` function was modified to handle the trace_ray return data in addition to ray and ray_pkg inputs.

Other changes include:

- Updated the doc for the profile module and added an update step during profile creation (issue 152).
- Fixed several subtle transformation issues and shortcomings when rendering cemented elements (issue 153). Reworked the :func:`~.elem.transform.compute_global_coords` function.
- Changed the paraxial NA to be ref index * slope, in keeping with OSLO and the literature (issue 154).
- Introduce :meth:`~.seq.sequential.SequentialModel.apply_scale_factor_over` in seq_model to permit scaling a range of surfaces. Helps with issue 156.
- Remove 'aperture' and 'field' from PupilSpec :attr:`~.PupilSpec.key` and FieldSpec :attr:`~.FieldSpec.key` key definitions; they were redundant information.

Thanks to @quentGit for correcting an error in the :mod:`~.waveabr` module.


Version 0.9.1
=============
Fix issue #150, need to constrain the version of numpy to < 2.0.0


Version 0.9.0
=============
A goal of this version is specify optical systems either via the sequential model, the |ybar| diagram, or the element/part tree description. The rayoptics app supports system entry by :meth:`~.seq.sequential.SequentialModel.add_surface` or :meth:`~.optical.opticalmodel.OpticalModel.add_lens` functions in the iPython console, by sketching a |ybar| diagram, and by opening existing optical models. This required a major update to functionality in the :class:`~.parax.paraxialdesign.ParaxialModel`, especially the |ybar| diagram functionality. 

To accomodate the different ways models may be built up interactively, a new internal approach of using a simple grammar (:mod:`~.elem.sgz2ele`) to take a sequential model as input and produce an element model/part tree description as output has been implemented. This should eliminate the need to use :meth:`~.optical.opticalmodel.OpticalModel.rebuild_from_seq` to get a correct lens layout. It works well with add_surface, add_lens, add_mirror, etc.

The :mod:`~.qtgui.rayopticsapp` module has been refined to better support use of the iPython console with the graphics windows and menu bar. A button, ``Refresh GUI``, was added to the iPython console to do an update_model() call and reshresh the graphics windows. The ``New`` item on the ``File`` menu opens an iPython console with an empty model and **opm**, **sm**, **osp**, etc. predefined. ``New Diagram`` will open an interactive |ybar| diagram window for sketching a diagram (after clicking the ``Sketch Diagram`` button). Panels for optical specification data are open and docked by default.

Optical systems often can have object or pupil planes at infinity. Even with IEEE floating point arithmetic, **inf** values cause problems in common calculations. Furthermore, it is common in imaging forming software to model an "infinite" object distance as a large number (e.g. 1e10); this convention is embedded in many imported models. The places in the code that dealt with infinite values were identified and guarded against bad outcomes using functions  :func:`~.util.misc_math.infinity_guard` and  :func:`~.util.misc_math.is_kinda_big`.

Effort was put into covering various corner cases that broke the code. Telecentric and afocal systems, object or image spaces with non-air materials, and models with a single surface (e.g. an eye model) were some of the areas that were addressed.

Other changes include:

- Changed :attr:`raytr.opticalspec.PupilSpec.key` `value_key` literal from ``pupil`` to ``epd`` for clarity.
- Added :meth:`~.optical.opticalmodel.OpticalModel.apply_scale_factor` method to OpticalModel, ParaxialModel, ElementModel, and the OpticalSpecs models.
- Fix the surface normal calculation for :class:`~.elem.profiles.XToroid` (issue 147).
- Look specifically for ``www.photonstophotos.net`` in stringified input (issue 145).
- Added wavefront/opd calculation based on an infinite reference sphere (issue 142).

Version 0.8.7
=============
Get dependencies correct for Python 3.8 on conda-forge and ReadTheDocs. Last release for Python 3.8


Version 0.8.6
=============
Implement a :class:`~.oprops.doe.DiffractionGrating` with ray trace and opd calculations. Add ability to draw single rays and control their color. Fix issues 102 (defocus application), 107 (Zemax import), 110 (defocus and image shift), 114 (Zemax import), 115 (ray aiming), 124 (immersed image), 125 (Zemax constant index), 127 (model modifications), 138 (local loggers).
There were a variety of small fixes and improvements, e.g. to :func:`~listobj`.
Thanks to @mpetroff and @dominikonysz for their contributions to **ray-optics**.


Version 0.8.5
=============
Fix crashing issue #101. Adjustments to the :mod:`~.medium` package. Remove many work files that weren't in the repo.


Version 0.8.4
=============
Improve :meth:`~.seq.sequential.SequentialModel.add_surface` handling of materials and semi-diameter.
Maintenance work on ray trace related issues, esp. aperture handling for list and grid traces.
Move functionality in :mod:`~.seq.medium` into the :mod:`opticalglass` package and update code.
Remove jupyterlab and ipympl from the installation. These should be installed as desired by the user.
Fix a ray trace bug in XToroid profile; other miscellaneous bug fixes.
Revise doc to use RTD yaml file; other doc updates.


Version 0.8.3
=============
Add :func:`~.gui.appcmds.set_pupil` function to set the pupil specification based on the aperture stop size. Add a hole feature to :class:`~.elem.elements.Element` and :class:`~.elem.elements.Mirror`. InteractiveLayout will draw holes but they haven't yet been included in the sequential model or the ray trace. The :meth:`~.seq.sequential.SequentialModel.add_surface` method was fixed and enhanced to take additional inputs for materials and allow use of a `sd` keyword argument to specify the semi-diameter of the surface. A variety of bug fixes were made to import and saving of models as well as other fixes; see the update log for specifics. All of the .roa files on the distro were restored and saved with the current version. At the same time, the long deprecated Pupil and Field enums were removed from the code. Finally, the build and packaging of **ray-optics** was updated to use current python technologies.

Version 0.8.2
=============
Fixes to drawing cemented element flats when the inner surfaces intersect the outer flats. Accept URLs from the OpticalBench database as arguments to :func:`~.gui.appcmds.open_model`. Enhance :func:`~listobj` output for :class:`~.seq.sequential.SequentialModel` and :class:`~.seq.gap.Gap`. Bump pyqt5 compatibility to v5.15.

Version 0.8.1
=============
Major review and cleanup of the :mod:`~.raytr` package. Moved wave aberration calculations and support to new module :mod:`~.raytr.waveabr`.  

Integrate :meth:`~.seq.interface.Interface.point_inside` query into the ray trace and add a new ray trace exception, :class:`~.raytr.traceerror.TraceRayBlockedError`. Use this to implement a pupil ray-based vignetting calculation and provide :func:`~.gui.appcmds.set_vignetting` to set the vignetted pupil rays for each field point and the companion function, :func:`~.gui.appcmds.set_apertures`, to set clear apertures based on the vignetted pupils at each field. to clear aperture. These new calculations are in the :mod:`~.raytr.vigcalc` module.

The addition of the :func:`~.gui.appcmds.set_vignetting` function uncovered a bug in the ray plots; they didn't handle the case of needing to overfill the paraxial entrance pupil in wide angle lenses. The ray plots, :class:`~.mpl.axisarrayfigure.RayFanFigure` and :class:`~.mpl.analysisfigure.RayFanPlot`, were updated to correctly display overfilled pupils, and an example lens was added to demonstrate the issue.

A bug affecting the display of flats on concave surfaces was fixed. The cell phone lens example was updated to explain some of the controls added for greater control over lens element display.

Numerous updates to the documentation.

Version 0.8.0
=============
Add :meth:`~.optical.opticalmodel.OpticalModel.flip` function to :class:`~.optical.opticalmodel.OpticalModel`. Can flip a range of surfaces, flip an element or an assembly. 

Added :class:`~.elem.elements.Part` as an ABC for elements, et al, and added an :class:`~.elem.elements.Assembly` class allow groups of `Parts`. Parts and Profiles are now saved and restored correctly. 

Revise :meth:`~.optical.opticalmodel.OpticalModel.update_model` operation, including splitting out derived property calculations into a separate :meth:`~.optical.opticalmodel.OpticalModel.update_optical_properties` function. Added `src_model` keyword to update_model so that submodels can filter on who originated the update.

Upgrade the |nubar| diagram to support the layer view control used in the |ybar| diagram. An `assembly` layer was added to the diagrams. The object/image conjugate shift on |ybar| diagram was fixed for the case of a virtual object.

Add `analysis_results` dictionary to submodels maintained by the `OpticalModel`. Move first order data there under keyword `parax_data`. Going forward, this is the natural place to put results from the different evaluators in the :mod:`~.raytr.analyses` module.

Removed convoluted ray offset calculation for `InteractiveLayout` and let the view bounds handle it.
Included the effect of defocus in transverse aberration plots (fixes bug in implementation).

Version 0.7.5
=============
Implemented a formatted object output protocol, the listobj_str() method, that returns a formatted, descriptive string about the object's contents. This supersedes the use of listobj as a per-class method (see 0.7.1 release notes). This protocol was implemented in particular for profiles and parts of the optical specification. Added a console print function :func:`~listobj` that will work with objects that implement a listobj_str() method; the fallback is to invoke repr(obj). Thanks to @asebian and @rlabs-oss for their issue reports and @dibyendumajumdar for continued discussions.

Version 0.7.4
=============
Pass 2 on the |ybar| diagram layer capability. Handle thick elements and include a top level 'sys' layer. Fix insertion of system from file. Add support for models from the `OpticalBenchHub <https://www.photonstophotos.net/GeneralTopics/Lenses/OpticalBench/OpticalBenchHub.htm>`_ portion of Bill Claff's `PhotonsToPhotos <https://www.photonstophotos.net/>`_ website. Support odd polynomial surfaces in Zemax import. Added additional control over the use of flats when drawing lens elements, see ray-optics notebook `Cell Phone lens <https://github.com/mjhoptics/ray-optics-notebooks/blob/master/Cell%20Phone%20lens.ipynb>`_ for an example. Thanks also to @wuffi for contributing 2 fixes to make the interactive ray-optics app more robust.

Version 0.7.3
=============
Miscellaneous bug fixes, see check-in comments.

Version 0.7.2
=============
Add RayFans to interactive layout. Add a multiple layer |ybar| diagram capability. Works well for thin lens systems. Systems with thick lenses (e.g. Double Gauss) don't work well. Fixes in the :mod:`~.raytr.analyses` module include getting the sign right for defocused transverse aberrations and using the image gap distance instead of parax img_dist for reference sphere definition. Miscellaneous fixes.

Version 0.7.1
=============
Switch to use of strings for :class:`~.optical.model_enums.DecenterType` and :class:`~.optical.model_enums.DimensionType`. Use of the types is now deprecated. Add listobj() methods to all :mod:`~.elem.profiles` classes that lists all of the data for each profile type. Completed the fix of `opticalglass` issue #5 by using `difflib` to find better matches; seems to almost eliminate need for .smx files. Provided an alternative replacement glass spec for .smx files that requires just the glass and catalog names. Miscellaneous other fixes.

Version 0.7.0
=============

Add :class:`~.elem.parttree.PartTree` to :class:`~.optical.opticalmodel.OpticalModel`. Many changes to get SequentialModel, ElementModel and PartTree to work in sync. Add :class:`~.elem.elements.CementedElement` element to :mod:`~.elem.elements` module. Fix paraxial image distance calculation. Fix DOE implementation for DOE on exiting side of element. Add fractional field height specification to :class:`~.raytr.opticalspec.FieldSpec`. Add a mapping interface to :class:`~.optical.opticalmodel.OpticalModel` and :class:`~.raytr.opticalspec.OpticalSpecs`; allows submodel access like opm['sm'] and opm['osp]['fov']. Add :meth:`~.seq.sequential.SequentialModel.reverse_path` method to facilitate reverse ray tracing through the model. Add semi-diameter to list of items in :meth:`~.seq.sequential.SequentialModel.add_surface` method. Add exception handling for ray trace errors in :mod:`~.raytr.analyses` module. Many miscellaneous fixes and improvements.

Version 0.6.5
=============

Fixes for ray normal calculation in Even and Radial polynomial profiles, the ray trace was incorrect for these surfaces. Add do_aperture flag to :class:`~.seq.sequential.SequentialModel` to control aperture setting by :meth:`~.seq.sequential.SequentialModel.update_model`. Miscellaneous other fixes and enhancements.

Version 0.6.4
==============

Rework :mod:`~.mpl.analysisfigure` module to better separate axis vs figure functions. Simplify default layout options, but support more complicated layouts as needed. Correctly update the model when dropping a glass onto the |ybar| diagram. Fix issue #18 and update documentation.

Version 0.6.3
==============

Fixes for aspheric ray intersection calculation, implemented by reworking the :mod:`~.elem.profiles` module. Generalize the :func:`~.elem.elements.create_mirror` fct to accept a :class:`~.elem.profiles.SurfaceProfile` instance in addition.

Version 0.6.2
==============

Fixes for .zmx import. Refactor GlassHandler to allow CODEV import to use as well.

Version 0.6.1
==============

Interpret Zemax .zmx files as having UTF-16 encoding, with a fallback to UTF-8. Push classes in :mod:`~.seq.medium` towards a common set of queries with :mod:`opticalglass` glass and catalog classes. Add :func:`~.raytr.analyses.trace_list_of_rays` specified by starting point, direction cosines and wavelength.

Version 0.6.0
==============

Add import of Zemax .zmx files. Demonstrate access to lens-designs.com collection of zemax files and use that to test import capability. Build a strong tie with :mod:`opticalglass` package, including drag and drop from a glassmap view to the layout, yybar and lens table views. Add x, y-toroids. Add a paraxial-based set vig calculation to Tools menu. Set up a protocol for aperture checking in the ray trace. Add conda installation via conda-forge.

Version 0.5.0
==============

Refactor/repackage the :mod:`~.optical` subpackage into more manageable chunks. New subpackages include :mod:`~.seq`, :mod:`~.raytr`, :mod:`~.parax`, :mod:`~.elem` and :mod:`~.oprops`. The :mod:`~.layout` and :mod:`~.diagram` modules were moved from the :mod:`~.gui` subpackage to the :mod:`~.elem` and :mod:`~.parax` subpackages, respectively. The refactoring broke reading of existing .roa files; a preprocessing step was added in the :mod:`~.roafile` module to handle mapping the old structure to the new one. Documentation and Jupyter notebooks were updated to reflect the changes. 

Version 0.4.12
==============

Misc bug fixes and improvements: draw mirror substrates consistently, old .roa files should be resaved. Consolidate InteractiveFigure color and lw control in dicts, and retrofit layout usage to match diagrams. 

Version 0.4.11
==============

Add single field monochromatic PSF calc and display. Add dashboards for jupyter notebook usage for refocusing, etc. Revisions to app manager protocols for wider usage, esp. including matplotlib Figures. Use brute force OPD calc instead of equally inclined chords until eic behavior with decentered systems is fixed.

Version 0.4.10
==============

Add add_from_file() method to OpticalModel to enable importing pieces of optical
models to a master model. Include a Jupyter notebook demonstration.

Version 0.4.9
=============

Add single panel refactoring of ray fan, spot diagram and wavefront map plots
that are designed for use in a scripting environment.

Version 0.4.8
=============

Bug fixes and doc improvements

Version 0.4.7
=============

UI improvements. For Qt app, add light and dark mode, add zoom, and pan
capabilities to InteractiveFigures, rework main menu bar. For Diagram, add
object and stop shifts, lens thickness and bending, replace node with system
option. Reorganize doc structure.

Version 0.4.6
=============

Add SpecSheet capability. V2 of |ybar| diagram with more editing capability

Version 0.4.0
=============

Add initial version of Interactive Layout (Live Layout). Work in Progress...

Version 0.3.5
=============

Update aperture import and layout rendering for y-z plane tilts and decenters

Version 0.3.0
=============

Rework packaging using PyScaffold and Versioneer. Deploy documentation to ReadTheDocs

Version 0.2.0
=============

first version of documentation support for rayoptics. continue refactoring

Version 0.1.5
=============

separate qt gui files from mpl and generic ones. add rayoptics script to start qt version

Version 0.1.0
=============

major update that reads CV files, supports Hoya, Ohara and Schott glass catalogs, computes paraxial and real rays, and has Qt UI for lens data table, 2D lens picture and glass map view.

