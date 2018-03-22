#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Qt GraphicsItem subclasses for 2D graphics views

Created on Tue Mar 20 10:39:27 2018

@author: Michael J. Hayford
"""

from PyQt5.QtCore import QPointF
from PyQt5.QtGui import (QColor, QPen, QPolygonF)
from PyQt5.QtWidgets import QGraphicsPolygonItem


class OpticalElement(QGraphicsPolygonItem):
    """ QGraphicsItem subclass for 2D lens element graphics """
    def __init__(self, element):
        super(OpticalElement, self).__init__()
        self.element = element
        self.update_shape()

        self.setBrush(QColor(*self.element.render_color))
        pen = QPen()
        pen.setCosmetic(True)
        self.setPen(pen)

    def update(self):
        super(OpticalElement, self).update()
        print("in update")

    def update_shape(self):
        poly = self.element.shape()
        polygon = QPolygonF()
        for p in poly:
            polygon.append(QPointF(p[0], p[1]))

        self.setPolygon(polygon)

        t = self.element.tfrm[1]
        self.setPos(QPointF(t[2], -t[1]))


class RayBundle(QGraphicsPolygonItem):
    """ QGraphicsItem subclass for ray bundle from a single field point """
    def __init__(self, seq_model, field_num, start_offset):
        super(RayBundle, self).__init__()
        self.seq_model = seq_model
        self.field_num = field_num
        self.start_offset = start_offset
        self.update_shape()

        self.setBrush(QColor(254, 197, 254, 64))  # magenta, 25%
        pen = QPen()
        pen.setCosmetic(True)
        self.setPen(pen)

    def update_shape(self):
        offset = self.start_offset
        seq_model = self.seq_model
        tfrms = seq_model.transforms
        rayset = seq_model.trace_boundary_rays_at_field(self.field_num)

        # If the object distance (tfrms[0][1][2]) is greater than the
        #  start_offset, then modify rayset start to match start_offset.
        # Remember object transformation for resetting at the end.
        tfrtm0 = tfrms[0]
        try:
            if abs(tfrms[0][1][2]) > self.start_offset:
                r, t = seq_model.setup_shift_of_ray_bundle(offset)
                tfrms[0] = (r, t)
                seq_model.shift_start_of_ray_bundle(rayset, offset, r, t)

            poly1 = []
            for i, r in enumerate(rayset[3][0][0:]):
                rot, trns = tfrms[i]
                p = rot.dot(r[0]) + trns
    #            print(i, r[0], rot, trns, p)
                poly1.append(QPointF(p[2], -p[1]))

            poly2 = []
            for i, r in enumerate(rayset[4][0][0:]):
                rot, trns = tfrms[i]
                p = rot.dot(r[0]) + trns
    #            print(i, r[0], rot, trns, p)
                poly2.append(QPointF(p[2], -p[1]))

            poly2.reverse()
            poly1.extend(poly2)
            polygon = QPolygonF()
            for p in poly1:
                polygon.append(p)

            self.setPolygon(polygon)

        finally:
            tfrms[0] = tfrtm0
