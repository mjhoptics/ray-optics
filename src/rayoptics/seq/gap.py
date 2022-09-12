#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2017 Michael J. Hayford
""" Gap container class

.. Created on Fri Sep 15 17:06:17 2017

.. codeauthor: Michael J. Hayford
"""


from opticalglass import opticalmedium as om


def gap_medium(g) -> str:
    """ return a formatted string of the gap's medium. """
    m = g.medium
    return f"{m.name():<12} {m.__class__.__name__:<14} {m.__module__:<27}"


class Gap:
    """ Gap container class.

    The gap class represents the space between 2 surfaces. It contains the
    media definition for the space and a (z) displacement between the
    adjacent surfaces.

    The most common use case is an optical system with surfaces centered on a
    common axis. The Gap structure implements this case in the simplest manner.
    More complicated transformations between surfaces are implemented using
    transformations associated with the surfaces themselves.

    Attributes:
        thi: the length (along z) of the gap
        medium: a :class:`~optical.medium.Medium` or a catalog glass instance
    """

    def __init__(self, t=0.0, med=om.Air()):
        self.thi = t
        self.medium = med

    def __repr__(self):
        return "Gap(t=%r, medium=%r)" % (self.thi, self.medium)

    def listobj_str(self):
        o_str = f"t={self.thi}"
        if self.medium is not None:
            o_str += ";    "
            name = self.medium.name()
            catalog = self.medium.catalog_name()
            if len(catalog) > 0:
                o_str += f"medium: {name}, {catalog}\n"
            else:
                o_str += f"medium: {name}\n"
        else:
            o_str += "\n"

        return o_str

    def sync_to_restore(self, seq_model):
        if hasattr(self.medium, 'sync_to_restore'):
            self.medium.sync_to_restore()

    def apply_scale_factor(self, scale_factor):
        self.thi *= scale_factor
