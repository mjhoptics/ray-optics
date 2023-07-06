#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
""" Definitions for user interface actions

    `Action` classes are those that contain a `dict` called ``actions``. The
    ``actions`` dict contains functions that are executed for UI actions 
    corresponding to the keyword.

    Action keywords include:
        'press'
        'drag'
        'release'

    The argument lists for the called functions include ``fig`` and ``event``.
    Duck typing is relied on loosen the binding with a particular GUI
    implementation - or at least that's the theory.

    ``fig`` must support a function refresh_gui()

.. Created on Wed Apr 24 14:26:14 2019

.. codeauthor: Michael J. Hayford
"""

import math

from opticalglass import glassfactory as gfact

kwupdate = {
    'build': 'update',
}

class Action():
    """ Action built on a set/get function pair """
    def __init__(self, getf, setf):
        self.getf = getf
        self.setf = setf
        self.cur_value = None
        self.new_value = None
        self.actions = {}

        def on_select(fig, event):
            self.cur_value = getf()
            return self.cur_value
        self.actions['press'] = on_select

        def on_edit(fig, event, new_value):
            setf(self.cur_value, new_value)
            fig.refresh_gui(**kwupdate)
        self.actions['drag'] = on_edit

        def on_release(fig, event):
            self.new_value = getf()
            fig.refresh_gui(**kwupdate)
        self.actions['release'] = on_release


class AttrAction():
    """ Action built on an object/attribute pair """
    def __init__(self, obj, attr):
        self.object = obj
        self.attr = attr
        self.cur_value = getattr(self.object, self.attr, None)
        self.new_value = None
        self.actions = {}

        def on_select(fig, event):
            self.cur_value = getattr(self.object, self.attr, None)
#            print('AttrAction.on_select:', self.attr, self.cur_value)
            return self.cur_value
        self.actions['press'] = on_select

        def on_edit(fig, event, delta_value):
            setattr(self.object, self.attr, self.cur_value+delta_value)
#            print('AttrAction.on_edit:', self.attr, self.cur_value+delta_value)
            fig.refresh_gui(**kwupdate)
        self.actions['drag'] = on_edit

        def on_release(fig, event):
            self.new_value = getattr(self.object, self.attr, None)
            fig.refresh_gui(**kwupdate)
        self.actions['release'] = on_release


class SagAction():
    """ Action to set the profile_cv via surface sag (related to x,y input) """
    def __init__(self, surf):
        self.surf = surf
        self.cur_value = None
        self.new_value = None
        self.actions = {}

        def on_select(fig, event):
            self.cur_value = self.surf.z_sag((event.x, event.ydata))
#            print('SagAction.on_select:', self.cur_value)
            return self.cur_value
        self.actions['press'] = on_select

        def on_edit(fig, event, value):
            self.surf.set_z_sag(value)
#            cv = self.surf.calc_cv_from_zsag(value)
#            print('SagAction.on_edit (x, y, cv):', value, cv)
            fig.refresh_gui(**kwupdate)
        self.actions['drag'] = on_edit

        def on_release(fig, event):
            self.new_value = self.surf.z_sag((event.x, event.ydata))
            fig.refresh_gui(**kwupdate)
        self.actions['release'] = on_release


class BendAction():
    """ Action to bend a lens element, using x component of user input """
    def __init__(self, e):
        self.ele = e
        self.select_pt_lcl = None
        self.cv_orig = None
        self.cv_new = None
        self.actions = {}

        def sag(cv, y):
            if cv != 0.0:
                r = 1/cv
                adj = math.sqrt(r*r - y*y)
                return r*(1 - abs(adj/r))
            else:
                return 0

        def on_select(fig, event):
            self.select_pt_lcl = event.lcl_pt
            self.cv_orig = self.ele.reference_interface().profile_cv
            self.xsag_orig = sag(self.cv_orig, event.lcl_pt[1])
#            print('on_select: {:.3f} {:.3f} {:.5f}'
#                  .format(event.lcl_pt[0], self.xsag_orig, self.cv_orig))
            return self.select_pt_lcl
        self.actions['press'] = on_select

        def on_edit(fig, event, lcl_pt):
            cv1 = self.ele.s1.profile_cv
            cv2 = self.ele.s2.profile_cv
            x, y = lcl_pt
            xsag_new = self.xsag_orig + (x - self.select_pt_lcl[0])
            new_pt = xsag_new, lcl_pt[1]
            cv1_new = self.ele.reference_interface().calc_cv_from_zsag(new_pt)
#            print('on_edit: {:.3f} {:.3f} {:.5f}'.format(x, xsag_new, cv1_new))
            delta_cv = cv1_new - cv1
            cv2_new = cv2 + delta_cv
            self.ele.s1.profile_cv = cv1_new
            self.ele.s2.profile_cv = cv2_new
            fig.refresh_gui(**kwupdate)
        self.actions['drag'] = on_edit

        def on_release(fig, event):
            self.cv_new = self.ele.reference_interface().profile_cv
            xsag = sag(self.cv_new, event.lcl_pt[1])
#            print('on_release: {:.3f} {:.3f} {:.5f}'
#                  .format(event.lcl_pt[0], xsag, self.cv_new))
            fig.refresh_gui(**kwupdate)
        self.actions['release'] = on_release


class ReplaceGlassAction():
    """ Action for replacing an element's glass from a drag/drop action. """

    def __init__(self, gap, update=True):
        self.gap = gap
        self.update = update
        self.actions = {}

        def null_action(fig, event):
            pass
        self.actions['press'] = null_action
        self.actions['drag'] = null_action
        self.actions['release'] = null_action

    def __call__(self, fig, event):
        mime = event.mimeData()
        # comma separated list
        glass_name, catalog_name = mime.text().split(',')
        mat = gfact.create_glass(glass_name, catalog_name)
        self.gap.medium = mat
        if self.update:
            fig.refresh_gui(**kwupdate)
