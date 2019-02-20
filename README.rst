==========
ray-optics
==========

Installation
------------

rayoptics makes use of PyQt5 for the **desktop** version GUI. This software is automatically installed and available when using `Anaconda <https://www.anaconda.com/>`_. Alternatively, PyQt5 may be downloaded using ``pip3``

.. code::

    > pip3 install pyqt5

rayoptics itself is installed using ``pip3`` as well.

.. code::

    > pip3 install rayoptics

Documentation
-------------

The documentation for **ray-optics** is hosted at `Read the Docs <https://rayoptics.readthedocs.io>`_

Tools for image forming optical design and analysis
---------------------------------------------------

The **ray-optics** project has several goals, both in the optical domain and
in the software development domain

* Rethink how image forming optical calculations are done absent historical
  constraints on computer speed and memory
* Serve as a reference implementation of basic ray tracing and wavefront
  analysis algorithms
* Leverage Python libraries to avoid reinventing the wheel
* Look for architectures that can be deployed on the desktop (e.g. using Qt),
  using (Jupyter) notebooks, and in mobile environments

Image forming optical design and analysis was one of the first domain areas to
be implemented on the first electronic computing machines. The calculations
were so voluminous and tedious that automation was a tremendous benefit. The
computational workflow initially followed that used by the human computers of
the day. The capabilities of electronic computers were also extremely limited.

Computers are vastly more powerful now than they were when the venerable
CODE V program was initially developed. Zemax and Oslo date from the early
IBM PC days, when both speed and memory were limiting factors. The software
tools available then were limited as well. In order to gain acceptable
performance, compiled languages such as C and FORTRAN were required. Graphical
user interfaces were also expensive to develop and were often considered
secondary in importance by vendors.

Optical calculation technology can be considered a mature field. There is a
long history in the literature of investigations tying optical theory to
calculation approaches that maintain accuracy, handle corner cases, etc. Much
of this knowledge is embedded in production optical design software; having it
available in an open source implementation brings this knowledge out of the
literature and makes it accessible to students and researchers.

The advent of scripting environments like Python make a fresh look at optical
design calculations worthwhile. Python has many advantages for scientific and
engineering applications, including libraries for numerical calculations, a
variety of graphics options including a very full featured plotting library,
and support for data science applications that look promising for managing
data generated during extensive analysis of optical systems. There is also
good support for unit testing of Python software, as well as debugging and
performance profiling.

Finally, computational notebooks, such as those of the Jupyter project,
provide the ability to document **reproducibly** the implementation of an
algorithm and results produced by the algorithm.

Note
----

This project has been set up using PyScaffold 3.1. For details and usage information on PyScaffold see https://pyscaffold.org/.
