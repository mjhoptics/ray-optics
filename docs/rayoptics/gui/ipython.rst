******************
iPython setup code
******************

.. code::

   import numpy as np
   import matplotlib.pyplot as plt
   import pandas as pd
   from pathlib import Path
   import rayoptics as ro

   # open an OpticalModel from a file
   root_pth = Path(ro.__file__).resolve().parent
   opm = ro.open_model(root_pth/"codev/test/landscape_lens.seq")

   # create a new OpticalModel
   #opm = ro.OpticalModel()
   sm = opm.seq_model
   osp = opm.optical_spec
   pm = opm.parax_model
