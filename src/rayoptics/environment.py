#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
""" script file providing an environment for using rayoptics

.. Created on Sun Feb 10 22:42:36 2019

.. codeauthor: Michael J. Hayford
"""

# initialization
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd

# import ipywidgets as widgets
# from ipywidgets import interact, interactive, fixed, interact_manual
# from IPython.display import display

# ray-optics
import rayoptics
from rayoptics.gui.appmanager import AppManager, ModelInfo
from rayoptics.gui.appcmds import create_new_model, open_model

# optical model
from rayoptics.optical.opticalmodel import OpticalModel

import rayoptics.optical.surface as srf
from rayoptics.optical.gap import Gap
from rayoptics.optical import elements
from rayoptics.optical.thinlens import ThinLens
from rayoptics.optical.profiles import (Spherical, Conic, EvenPolynomial,
                                        RadialPolynomial)
from rayoptics.optical.opticalspec import (WvlSpec, FieldSpec, Field,
                                           PupilSpec, FocusRange)
from rayoptics.optical.model_enums import (DimensionType, DecenterType)

# ray-optics first and third order
import rayoptics.optical.firstorder as fo
from rayoptics.optical.firstorder import compute_first_order
import rayoptics.optical.thirdorder as to
from rayoptics.optical.thirdorder import compute_third_order

# ray tracing
from rayoptics.optical.trace import (RayPkg, RaySeg, list_ray,
                                     trace, trace_base, trace_with_opd,
                                     trace_astigmatism)
import rayoptics.optical.raytrace as rt
from rayoptics.optical import analyses

# paraxial design
import rayoptics.optical.model_constants as mc

import rayoptics.mpl.paraxdgnfigure as pdfig
from rayoptics.mpl.paraxdgnfigure import ParaxialDesignFigure
from rayoptics.mpl.interactivediagram import InteractiveDiagram

# axis array figures
from rayoptics.mpl.axisarrayfigure import Fit
from rayoptics.mpl.axisarrayfigure import RayFanFigure
from rayoptics.mpl.axisarrayfigure import SpotDiagramFigure
from rayoptics.mpl.axisarrayfigure import WavefrontFigure

from rayoptics.mpl.analysisfigure import (AnalysisFigure,
                                          RayFanPlot, RayGeoPSF,
                                          Wavefront, DiffractionPSF)

from rayoptics.mpl import analysisplots

# lens layout
import rayoptics.mpl.lenslayoutfigure as lay
from rayoptics.mpl.lenslayoutfigure import LensLayoutFigure
from rayoptics.mpl.interactivelayout import InteractiveLayout

# opticalglass
from opticalglass.glassfactory import create_glass
from rayoptics.util.spectral_lines import get_wavelength