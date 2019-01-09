******************
Jupyter setup code
******************

.. code::

   # initialization
   import math
   import numpy as np
   import matplotlib as mpl
   import matplotlib.pyplot as plt
   from pathlib import Path
   import pandas as pd

   import ipywidgets as widgets
   from ipywidgets import interact, interactive, fixed, interact_manual
   from IPython.display import display

   # ray-optics
   import rayoptics as ro
   from rayoptics.gui.appmanager import AppManager, ModelInfo

   # optical model
   import rayoptics.optical.opticalmodel as om

   import rayoptics.optical.surface as srf
   from rayoptics.optical.gap import Gap
   from rayoptics.optical.thinlens import ThinLens
   from rayoptics.optical.profiles import Spherical, Conic
   from rayoptics.optical.opticalspec import FieldSpec, PupilSpec

   # ray-optics first and third order
   import rayoptics.optical.firstorder as fo
   from rayoptics.optical.firstorder import compute_first_order
   import rayoptics.optical.thirdorder as to
   from rayoptics.optical.thirdorder import compute_third_order

   # ray tracing
   from rayoptics.optical.trace import list_ray, trace, trace_base, trace_with_opd
   import rayoptics.optical.raytrace as rt

   # paraxial design
   import rayoptics.optical.model_constants as mc
   #from rayoptics.optical.model_constants import ht, slp, aoi
   #from rayoptics.optical.model_constants import pwr, tau, indx, rmd

   import rayoptics.mpl.paraxdgnfigure as pdfig
   from rayoptics.mpl.paraxdgnfigure import ParaxialDesignFigure, Dgm

   # axis array figures
   from rayoptics.mpl.axisarrayfigure import Fit
   from rayoptics.mpl.axisarrayfigure import RayFanFigure
   from rayoptics.mpl.axisarrayfigure import SpotDiagramFigure
   from rayoptics.mpl.axisarrayfigure import WavefrontFigure

   # lens layout
   import rayoptics.mpl.lenslayoutfigure as lay
   from rayoptics.mpl.lenslayoutfigure import LensLayoutFigure
