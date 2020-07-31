#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2020 Michael J. Hayford
"""Test module name replacement for .roa file reading

.. Created on Thu Jul 30 16:26:00 2020

.. codeauthor: Michael J. Hayford
"""


import unittest
import warnings
from pathlib import Path

import rayoptics as ro
from rayoptics.gui.appcmds import open_model
from rayoptics.optical.opticalmodel import OpticalModel


class RestoreOldFilesTestCase(unittest.TestCase):
    """Test module name replacement for .roa file reading

       This tc suite tests the compatability of .roa files written before the
       v0.5 repackaging.

       These tests ensure that the preprocessing of changed module locations
       works correctly. The failure the preprocessing fixes is tested as
       well by passing an empty mapping dictionary into open_model.
    """

    def test_cell_phone(self):
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        root_pth = Path(ro.__file__).resolve().parent
        opm = open_model(root_pth/'gui/tests/cell_phone_camera_old.roa')
        assert isinstance(opm, OpticalModel)

    def test_Sasian_Triplet(self):
        root_pth = Path(ro.__file__).resolve().parent
        opm = open_model(root_pth/'gui/tests/Sasian Triplet old.roa')
        assert isinstance(opm, OpticalModel)

    def test_cell_phone_fail(self):
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        root_pth = Path(ro.__file__).resolve().parent
        with self.assertRaises(ImportError):
            open_model(root_pth/'gui/tests/cell_phone_camera_old.roa',
                       mapping={})

    def test_Sasian_Triplet_fail(self):
        root_pth = Path(ro.__file__).resolve().parent
        with self.assertRaises(ImportError):
            open_model(root_pth/'gui/tests/Sasian Triplet old.roa',
                       mapping={})


if __name__ == '__main__':
    unittest.main(verbosity=3)
