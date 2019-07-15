.. currentmodule:: rayoptics.optical

***********
Ray Tracing
***********

Ray Tracing Overview
====================

   Real Ray Tracing is a fundamental operation in optical design and analysis. Ray tracing is two operations, repeated until reaching the image. They are:
   ::
   
      - find the next ray-interface intersection
      - scatter the ray based on the interface properties

   Sequential models of optical systems make the first operation very straightforward. No search for the next closest surface is performed, rather, the sequence of the interfaces is used to determine the next surface interface to calculate. Sequential models also have only one dominant optical behavior per interface. This means there will be a one to one relationship between the incident and exiting rays from an interface.

   Sequential ray tracing in **ray-optics** by default is configured to ray trace the seq_model of the OpticalModel. The :func:`raytrace.trace_raw` function is the building block of the ray trace system. The first argument to this function is an iterator that returns a tuple of interface, gap, refractive index, transform to the next interface and a z direction. Python provides a generator capability that could be used to programmatically generate a sequence of interfaces, for example for a ghost image analysis, without explicitly constructing a list of interfaces.