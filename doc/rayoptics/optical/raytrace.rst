.. currentmodule:: rayoptics.optical

***********
Ray Tracing
***********

.. _raytrace:

Ray Tracing Overview
====================

   Real Ray Tracing is a fundamental operation in optical design and analysis. Ray tracing is two operations, repeated until reaching the image. They are:
   ::
   
      - find the next ray-interface intersection
      - scatter the ray based on the interface properties

   Sequential models of optical systems make the first operation very straightforward. No search for the next closest surface is performed, rather, the sequence of the interfaces is used to determine the next surface interface to calculate. Sequential models also have only one dominant optical behavior per interface. This means there will be a one to one relationship between the incident and exiting rays from an interface.