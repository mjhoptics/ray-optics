***************
Scripting Usage
***************

The ``environment.py`` module imports many useful classes and functions. All the symbols defined in the module are intended to be imported into a rayoptics interactive session.

.. code::

   > from rayoptics.environment import *

Create a new OpticalModel

.. code::

   > opm = OpticalModel()

Alternatively, restore a model from a file

.. code::

   > opm = open_model('ritchey_chretien.roa')

CODE V sequence files can be imported

.. code::

   > opm = open_model('ag_dblgauss.seq')

Zemax .zmx files can be imported

.. code::

   > opm = open_model('US05831776-1.zmx')

Please note that the author doesn't have access to either CODE V or Zemax so if you find problems with the import, please create an issue in Git-Hub with the problem file and the expected results.

Whether created new or restored from a file, setting up the following names makes reuse of code snippets much easier

.. code::

   > sm = opm['seq_model']
   > osp = opm['optical_spec']
   > pm = opm['parax_model']
   > em = opm['ele_model']
   > pt = opm['part_tree']
   > ar = opm['analysis_results']

******************
iPython setup code
******************

One can use the flow above in the iPython environment. If you want to import fewer things into your iPython environment, you can import just the modules you need:

.. code::

   import numpy as np
   import matplotlib.pyplot as plt
   from rayoptics.optical.opticalmodel import OpticalModel
   from rayoptics.gui.appcmds import open_model
   
   # open an OpticalModel from a file
   opm = open_model("codev/test/landscape_lens.seq")

   # or create a new OpticalModel
   #opm = OpticalModel()
   sm = opm['seq_model']
   osp = opm['optical_spec']
   pm = opm['parax_model']
   em = opm['ele_model']
   pt = opm['part_tree']
   ar = opm['analysis_results']

******************
Jupyter setup code
******************

Depending on whether you are using Jupyter or JupyterLab, you may want to set a particular matplotlib back-end for your session. If the "out of the box" default isn't working, you can try changing it *at the beginning of your session*:

.. code::

   > %matplotlib widget

