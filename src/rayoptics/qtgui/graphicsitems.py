#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Qt GraphicsItem subclasses for 2D graphics views

.. Created on Tue Mar 20 10:39:27 2018

.. codeauthor: Michael J. Hayford
"""

from PyQt5.QtCore import QPointF
from PyQt5.QtGui import (QColor, QPen, QPolygonF)
from PyQt5.QtWidgets import QGraphicsPolygonItem

import rayoptics.gui.layout as layout


class OpticalElement(QGraphicsPolygonItem):
    """ QGraphicsItem subclass for 2D lens element graphics """
    def __init__(self, element):
        super().__init__()

        self.oe = layout.OpticalElement(element)
        self.update_shape()

        color = self.oe.render_color()
        self.setBrush(QColor(*color))
        pen = QPen()
        pen.setCosmetic(True)
        self.setPen(pen)

    def update(self):
        super().update()
        print("in update")

    def update_shape(self):
        poly, bbox = self.oe.update_shape()
        polygon = QPolygonF()
        for p in poly:
            # Qt canvas coordinates have 0,0 at top left of window, i.e. +y is
            #  pointing down
            polygon.append(QPointF(p[0], -p[1]))

        self.setPolygon(polygon)


class RayBundle(QGraphicsPolygonItem):
    """ QGraphicsItem subclass for ray bundle from a single field point """
    def __init__(self, opt_model, field_num, start_offset):
        super().__init__()
        fld, wvl, foc = opt_model.optical_spec.lookup_fld_wvl_focus(field_num)
        self.rb = layout.RayBundle(opt_model, fld, wvl, start_offset)
        self.update_shape()

        self.setBrush(QColor(254, 197, 254, 64))  # magenta, 25%
        pen = QPen()
        pen.setCosmetic(True)
        self.setPen(pen)

    def update_shape(self):
        poly, bbox = self.rb.update_shape()
        polygon = QPolygonF()
        for p in poly:
            # Qt canvas coordinates have 0,0 at top left of window, i.e. +y is
            #  pointing down
            polygon.append(QPointF(p[0], -p[1]))

        self.setPolygon(polygon)
