#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
""" Qt5 dialog box for ideal imager ui

.. Created on Tue Jun 18 13:09:19 2019

.. codeauthor: Michael J. Hayford
"""

import math

from rayoptics.optical import etendue
from rayoptics.optical.etendue import obj_img_set, fld_ape_set, fld_labels, ap_labels

from rayoptics.optical import specsheet

from rayoptics.util import dict2d
from rayoptics.util.dict2d import dict2D

from PyQt5.QtCore import Qt as qt
from PyQt5.QtWidgets import (QApplication, QDialog, QRadioButton,
                             QStackedWidget, QFormLayout, QGridLayout,
                             QGroupBox, QHBoxLayout, QLabel, QLineEdit,
                             QPushButton, QCheckBox, QVBoxLayout)


def value_to_text(value, fmt_str="{:> #.5f}"):
    if value is None:
        value_text = ''
    elif type(value) is str:
        value_text = value
    else:
        value_text = fmt_str.format(value)
    return value_text


class IdealImagerDialog(QDialog):
    NumGridRows = 3
    NumButtons = 4

    def __init__(self, conjugate_type, imager=None, imager_inputs=None,
                 **kwargs):
        super().__init__(**kwargs)

        self.dlog_attrs = {}

        self.conjugate_type = conjugate_type
        self.conjugate_box = self.createConjugateBox(itype=self.conjugate_type)

        # setup finite conjugate defaults
        self.imager_stack = {}
        figb = ImagerSpecGroupBox()
        self.imager_stack['finite'] = figb

        self.etendue_stack = {}
        fegb = EtendueGroupBox(self, 'finite')
        self.etendue_stack['finite'] = fegb

        # setup infinite conjugate defaults
        enabled_list = [False, False, False, False, True]
        chkbox_enabled_list = [False]*5
        imager_inputs = {'s': -math.inf, 'f': None}
        imager = specsheet.IdealImager(None, -math.inf, None, None, None)
        iigb = ImagerSpecGroupBox(enabled_list=enabled_list, imager=imager,
                                  chkbox_enabled_list=chkbox_enabled_list,
                                  imager_inputs=imager_inputs)
        self.imager_stack['infinite'] = iigb

        iegb = EtendueGroupBox(self, 'infinite')
        self.etendue_stack['infinite'] = iegb

        self.imager_groupbox_stack = QStackedWidget()
        self.imager_groupbox_stack.addWidget(iigb)
        self.imager_groupbox_stack.addWidget(figb)

        self.etendue_groupbox_stack = QStackedWidget()
        self.etendue_groupbox_stack.addWidget(iegb)
        self.etendue_groupbox_stack.addWidget(fegb)

        imager_groupbox = self.imager_stack[self.conjugate_type]
        self.imager_groupbox_stack.setCurrentWidget(imager_groupbox)

        etendue_groupbox = self.etendue_stack[self.conjugate_type]
        self.etendue_groupbox_stack.setCurrentWidget(etendue_groupbox)

        mainLayout = QHBoxLayout()
        mainLayout.addWidget(self.conjugate_box)
        mainLayout.addWidget(self.imager_groupbox_stack)
        mainLayout.addWidget(self.etendue_groupbox_stack)
        self.setLayout(mainLayout)

        self.setWindowTitle("Optical Spec Sheet")

    def get_imager_and_inputs(self):
        imager_groupbox = self.imager_stack[self.conjugate_type]
        return imager_groupbox.imager, imager_groupbox.imager_inputs

    def set_imager_and_inputs(self, imager, inputs):
        imager_groupbox = self.imager_stack[self.conjugate_type]
        imager_groupbox.imager = imager
        imager_groupbox.imager_inputs = inputs
        imager_groupbox.update_values()

    def createConjugateBox(self, itype='infinite'):
        conjugate_box = QGroupBox("Conjugates")
        layout = QVBoxLayout()

        infinite = QRadioButton('Infinite')
        self.dlog_attrs['infinite'] = infinite
        layout.addWidget(infinite)
        infinite.clicked.connect(lambda: self.change_conjugate('infinite'))

        finite = QRadioButton('Finite')
        self.dlog_attrs['finite'] = finite
        layout.addWidget(finite)
        finite.clicked.connect(lambda: self.change_conjugate('finite'))

        conjugate_box.setLayout(layout)

        self.dlog_attrs[itype].setChecked(True)

        return conjugate_box

    def change_conjugate(self, conj_type):
        print("change_conjugate:", conj_type)
        if self.conjugate_type is not conj_type:
            self.update_conjugate(conj_type)

    def update_conjugate(self, conj_type):
#        print(conj_type)
        prev_imager_groupbox = self.imager_groupbox_stack.currentWidget()
        self.conjugate_type = conj_type
        new_imager_groupbox = self.imager_stack[conj_type]
        if prev_imager_groupbox.imager is not None:
            if new_imager_groupbox.imager is None:
                new_imager_groupbox.imager = prev_imager_groupbox.imager
            else:
                prev_enabled_list, _ = prev_imager_groupbox.get_enabled_lists()
                new_enabled_list, _ = new_imager_groupbox.get_enabled_lists()
                new_imager = list(new_imager_groupbox.imager)
                for i, p in enumerate(prev_imager_groupbox.imager):
#                    print(i, p, prev_enabled_list[i], new_enabled_list[i])
                    # only transfer values if both items are enabled
                    if (prev_enabled_list[i] and new_enabled_list[i]):
                        new_imager[i] = p
                        key = prev_imager_groupbox.keys[i]
                        if key in prev_imager_groupbox.imager_inputs:
                            new_imager_groupbox.chkbox_change(qt.Checked, key)
                            new_imager_groupbox.imager_inputs[key] = p
                new_imager_groupbox.imager = specsheet.IdealImager(*new_imager)
        self.imager_groupbox_stack.setCurrentWidget(new_imager_groupbox)
        new_imager_groupbox.update_values(new_imager)

        etendue_groupbox = self.etendue_stack[conj_type]
        self.etendue_groupbox_stack.setCurrentWidget(etendue_groupbox)
        etendue_groupbox.update_values()

    def update_values(self):
        imager_groupbox = self.imager_stack[self.conjugate_type]
        imager = imager_groupbox.imager
        imager_inputs = imager_groupbox.imager_inputs
        conj_type = self.conjugate_type
        if conj_type == 'finite':
            imager_defined = True if imager.m is not None else False
        else:
            imager_defined = True if imager.f is not None else False

        etendue_groupbox = self.etendue_stack[conj_type]

        etendue_inputs = etendue_groupbox.etendue_inputs
        etendue_grid = etendue_groupbox.etendue_grid
        li = dict2d.num_items_by_type(etendue_inputs, fld_ape_set, obj_img_set)
        if imager_defined:
            if li['field'] == 1 and li['aperture'] == 1:
                # we have enough data to calculate all of the etendue grid
                etendue.do_etendue_via_imager(conj_type, imager_inputs, imager,
                                              etendue_inputs, etendue_grid)
            elif li['field'] == 1:
                # we have enough data to calculate all of the etendue grid
                row = dict2d.row(etendue_inputs, 'field')
                obj_img_key = 'object' if len(row['object']) else 'image'
                etendue.do_field_via_imager(conj_type, imager, etendue_inputs,
                                            obj_img_key, etendue_grid)
            elif li['aperture'] == 1:
                # we have enough data to calculate all of the etendue grid
                row = dict2d.row(etendue_inputs, 'aperture')
                obj_img_key = 'object' if len(row['object']) else 'image'
                etendue.do_aperture_via_imager(conj_type, imager, etendue_inputs,
                                               obj_img_key, etendue_grid)
        else:  # imager not specified
            if li['field'] == 2 or li['aperture'] == 2:
                # solve for imager
                ii = etendue.do_etendue_to_imager(etendue_inputs, etendue_grid)
                imager_inputs[ii[0]] = ii[1]
                imager_groupbox.update_values()
                # update etendue grid
                etendue.do_etendue_via_imager(conj_type, imager_inputs,
                                              imager_groupbox.imager,
                                              etendue_inputs, etendue_grid)
            else:  # don't wipe out inputs
                for fld_ape_key, fld_ape_value in etendue_inputs.items():
                    for obj_img_key, cell in fld_ape_value.items():
                        for key in cell:
                            etendue_grid[fld_ape_key][obj_img_key][key] = \
                                etendue_inputs[fld_ape_key][obj_img_key][key]

        etendue_groupbox.update_values()


class ImagerSpecGroupBox(QGroupBox):
    def __init__(self, keys=None, enabled_list=None, imager=None,
                 imager_inputs=None, chkbox_enabled_list=None, labels=None,
                 **kwargs):
        super().__init__(title='Imager specs', **kwargs)
        self.dlog_attrs = {}

        self.keys = keys if keys else specsheet.ideal_imager_keys
        self.labels = labels if labels else specsheet.ideal_imager_labels

        if imager is None:
            imager = specsheet.IdealImager(None, None, None, None, None)
        self.imager = imager

        self.imager_inputs = imager_inputs if imager_inputs else {}

        self.enabled_list = enabled_list
        self.chkbox_enabled_list = chkbox_enabled_list
        enabled_list, chkbox_enabled_list = self.get_enabled_lists()

        layout = QGridLayout()
        for i, key in enumerate(self.keys):
            label = QLabel(self.labels[i]+':')
            lineEdit = QLineEdit()
            checkBox = QCheckBox()
            row = i + 1
            layout.addWidget(label, row, 0)
            layout.addWidget(lineEdit, row, 1)
            layout.addWidget(checkBox, row, 2)
            self.dlog_attrs[key] = (label, lineEdit, checkBox)
            lineEdit.setText(value_to_text(imager[i]))
            lineEdit.setEnabled(enabled_list[i])
            # Use returnPressed signal rather than editingFinished. The latter
            #  fires on change of focus as well as return. We only want a
            #  signal when the lineEdit contents change and are committed by
            #  the user.
            lineEdit.returnPressed.connect(lambda k=key: self.value_change(k))

            if key in self.imager_inputs:
                if not checkBox.isChecked():
                    checkBox.setCheckState(qt.Checked)
            checkBox.setEnabled(chkbox_enabled_list[i])
            checkBox.stateChanged.connect(lambda state, k=key:
                                          self.chkbox_change(state, k))

        self.setLayout(layout)

    def value_change(self, imager_key):
        label, lineEdit, checkBox = self.dlog_attrs[imager_key]
#        print(" value_change: {} {} '{}'".format(str(self.imager_inputs),
#              imager_key, lineEdit.text()))
        try:
            value = float(lineEdit.text())
        except ValueError:
            return
        else:
            if imager_key in self.imager_inputs or len(self.imager_inputs) < 2:
                self.imager_inputs[imager_key] = value
                if not checkBox.isChecked():
                    checkBox.setCheckState(qt.Checked)
            self.update_values()

    def update_values(self, values=None, enabled_list=None):
#        print('update_values:', self.imager_inputs)
        if values is None:
            inputs = {k: v for (k, v) in self.imager_inputs.items()
                      if v is not None}
            imager = specsheet.ideal_imager_setup(**inputs)
            if imager is not None:
                self.imager = imager

            if self.imager:
                value_list = self.imager
            else:
                value_list = [None]*len(self.keys)
        else:
#            print('update_values: values', values)
            value_list = values

        if self.imager is not None:
            # refill gui widgets and set enabled state
            for i, key in enumerate(self.keys):
                value = value_list[i]
                label, lineEdit, checkBox = self.dlog_attrs[key]
                lineEdit.setText(value_to_text(value))

        self.update_checkboxes()

    def chkbox_change(self, state, imager_key):
        label, lineEdit, checkBox = self.dlog_attrs[imager_key]
        checked = state == qt.Checked
#        print("chkbox_change:", self.imager_inputs, imager_key, checked)
        if checked:
            try:
                value = float(lineEdit.text())
            except ValueError:
                value = None
            finally:
                self.imager_inputs[imager_key] = value
                if not checkBox.isChecked():
                    checkBox.setChecked(True)
        else:
            if imager_key in self.imager_inputs:
                del self.imager_inputs[imager_key]
                if checkBox.isChecked():
                    checkBox.setChecked(False)
        self.update_checkboxes()
#        print("  chkbox_exit:", self.imager_inputs, imager_key, checked)

    def get_enabled_lists(self):
        if self.enabled_list is not None:
            enabled_list = self.enabled_list
        else:
            if len(self.imager_inputs) < 2:
                enabled_list = [True]*5
            else:
                enabled_list = [False]*5
            for key in self.imager_inputs:
                enabled_list[self.keys.index(key)] = True

        if self.chkbox_enabled_list is not None:
            chkbox_enabled_list = self.chkbox_enabled_list
        else:
            chkbox_enabled_list = enabled_list

        return enabled_list, chkbox_enabled_list

    def update_checkboxes(self):
        enabled_list, chkbox_enabled_list = self.get_enabled_lists()

        for i, key in enumerate(self.keys):
            label, lineEdit, checkBox = self.dlog_attrs[key]

            lineEdit.setEnabled(enabled_list[i])
            checkBox.setEnabled(chkbox_enabled_list[i])


class EtendueGroupBox(QGroupBox):
    def __init__(self, parent, itype,
                 object_inputs=None, image_inputs=None, **kwargs):
        super().__init__(title='optical invariant', **kwargs)

        self.parent = parent

        self.groupboxes = dict2D(fld_ape_set, obj_img_set)
        self.etendue_inputs = dict2D(fld_ape_set, obj_img_set)
        self.etendue_grid = dict2D(fld_ape_set, obj_img_set)

        if itype is 'infinite':
            gb00 = SpaceGroupBox(parent=self, title='object space',
                                 obj_img='object', fld_ape='field',
                                 keys=(fld_labels[1],))
            self.groupboxes['field']['object'] = gb00
            self.etendue_grid['field']['object'] = dict([
                    (fld_labels[1], None)])
            gb10 = SpaceGroupBox(parent=self, title='',
                                 obj_img='object', fld_ape='aperture',
                                 keys=(ap_labels[0],))
            self.groupboxes['aperture']['object'] = gb10
            self.etendue_grid['aperture']['object'] = dict([
                    (ap_labels[0], None)])
            gb01 = SpaceGroupBox(parent=self, title='image space',
                                 obj_img='image', fld_ape='field',
                                 keys=(fld_labels[0],))
            self.groupboxes['field']['image'] = gb01
            self.etendue_grid['field']['image'] = dict([
                    (fld_labels[0], None)])
            gb11 = SpaceGroupBox(parent=self, title='',
                                 obj_img='image', fld_ape='aperture',
                                 keys=(ap_labels[2], ap_labels[1]))
            self.groupboxes['aperture']['image'] = gb11
            self.etendue_grid['aperture']['image'] = dict([
                    (ap_labels[2], None), (ap_labels[1], None)])

        else:
            gb00 = SpaceGroupBox(parent=self, title='object space',
                                 fld_ape='field', obj_img='object',
                                 keys=(fld_labels[0],))
            self.groupboxes['field']['object'] = gb00
            self.etendue_grid['field']['object'] = dict([
                    (fld_labels[0], None)])
            gb10 = SpaceGroupBox(parent=self, title='',
                                 fld_ape='aperture', obj_img='object',
                                 keys=(ap_labels[2], ap_labels[1]))
            self.groupboxes['aperture']['object'] = gb10
            self.etendue_grid['aperture']['object'] = dict([
                    (ap_labels[2], None), (ap_labels[1], None)])
            gb01 = SpaceGroupBox(parent=self, title='image space',
                                 fld_ape='field', obj_img='image',
                                 keys=(fld_labels[0],))
            self.groupboxes['field']['image'] = gb01
            self.etendue_grid['field']['image'] = dict([
                    (fld_labels[0], None)])
            gb11 = SpaceGroupBox(parent=self, title='',
                                 fld_ape='aperture', obj_img='image',
                                 keys=(ap_labels[2], ap_labels[1]))
            self.groupboxes['aperture']['image'] = gb11
            self.etendue_grid['aperture']['image'] = dict([
                    (ap_labels[2], None), (ap_labels[1], None)])

        layout = QGridLayout()
        layout.addWidget(gb00, 0, 0)
        layout.addWidget(gb10, 1, 0)
        layout.addWidget(gb01, 0, 1)
        layout.addWidget(gb11, 1, 1)

        self.setLayout(layout)

    def get_imager(self):
        imager, inputs = self.parent.get_imager_and_inputs()
        if self.parent.conjugate_type == 'finite':
            return imager if imager.m is not None else None
        else:
            return imager if imager.f is not None else None

    def value_change(self, fld_ape, obj_img, key):
        """ callback routine for item value widget """
        inputs = self.etendue_inputs[fld_ape][obj_img]
        gb = self.groupboxes[fld_ape][obj_img]
        label, lineEdit, checkBox = gb.dlog_attrs[key]
        print("value_change: {} {} {}".format(obj_img, fld_ape, key))
#        print(" value_change: {} {} '{}'".format(str(self.etendue_inputs),
#              key, lineEdit.text()))
        try:
            value = float(lineEdit.text())
        except ValueError:
            return
        else:
            if (key in inputs or len(inputs) < 2):
                inputs[key] = value
                if not checkBox.isChecked():
                    checkBox.setCheckState(qt.Checked)
            self.parent.update_values()

    def chkbox_change(self, state, fld_ape, obj_img, key):
        """ callback routine for item checkbox """
        inputs = self.etendue_inputs[fld_ape][obj_img]
        gb = self.groupboxes[fld_ape][obj_img]
        label, lineEdit, checkBox = gb.dlog_attrs[key]
        checked = state == qt.Checked
#        print("chkbox_change:", self.etendue_inputs, key, checked)
        if checked:
            try:
                value = float(lineEdit.text())
            except ValueError:
                value = None
            finally:
                inputs[key] = value
                if not checkBox.isChecked():
                    checkBox.setChecked(True)
        else:
            if key in inputs:
                del inputs[key]
                if checkBox.isChecked():
                    checkBox.setChecked(False)
        self.update_checkboxes()
#        print("  chkbox_exit:", self.etendue_inputs, key, checked)

    def get_enabled_lists(self):
        imager = self.get_imager()

        if self.enabled_list is not None:
            enabled_list = self.enabled_list
        else:
            if len(self.etendue_inputs) < 2:
                enabled_list = [True]*7
            else:
                enabled_list = [False]*7
            for key in self.etendue_inputs:
                enabled_list[self.keys.index(key)] = True

        if self.chkbox_enabled_list is not None:
            chkbox_enabled_list = self.chkbox_enabled_list
        else:
            chkbox_enabled_list = enabled_list

        return enabled_list, chkbox_enabled_list

    def update_values(self):
        etendue_grid = self.etendue_grid

        for fld_ape_key, fld_ape_value in etendue_grid.items():
            for obj_img_key, cell in fld_ape_value.items():
                sgb = self.groupboxes[fld_ape_key][obj_img_key]
                sgb.update_values(cell)

    def update_checkboxes(self):
        pass
#        enabled_list, chkbox_enabled_list = self.get_enabled_lists()
#
#        for i, key in enumerate(self.keys):
#            label, lineEdit, checkBox = self.dlog_attrs[key]
#
#            lineEdit.setEnabled(enabled_list[i])
#            checkBox.setEnabled(chkbox_enabled_list[i])


class SpaceGroupBox(QGroupBox):
    def __init__(self, parent, title, fld_ape, obj_img, keys,
                 enabled_list=None, labels=None,
                 item_inputs=None, chkbox_enabled_list=None, **kwargs):
        super().__init__(title=title, **kwargs)
        self.dlog_attrs = {}

        self.parent = parent
        self.fld_ape = fld_ape
        self.obj_img = obj_img

        self.keys = keys
        if labels is None:
            self.labels = keys
        else:
            self.labels = labels

        self.item_inputs = item_inputs if item_inputs else {}

        self.enabled_list = enabled_list
        self.chkbox_enabled_list = chkbox_enabled_list
        enabled_list, chkbox_enabled_list = self.get_enabled_lists()

        layout = QGridLayout()

        for i, key in enumerate(self.keys):
            label = QLabel(self.labels[i]+':')
            lineEdit = QLineEdit()
            checkBox = QCheckBox()
            row = i + 1
            layout.addWidget(label, row, 0)
            layout.addWidget(lineEdit, row, 1)
            layout.addWidget(checkBox, row, 2)
            self.dlog_attrs[key] = (label, lineEdit, checkBox)
            
            if key in self.item_inputs:
                lineEdit.setText(value_to_text(self.item_inputs[key]))
            lineEdit.setEnabled(enabled_list[i])
            # Use returnPressed signal rather than editingFinished. The latter
            #  fires on change of focus as well as return. We only want a
            #  signal when the lineEdit contents change and are committed by
            #  the user.
            lineEdit.returnPressed.connect(lambda k=key:
                                           parent.value_change(fld_ape,
                                                               obj_img, k))

            if key in self.item_inputs:
                if not checkBox.isChecked():
                    checkBox.setCheckState(qt.Checked)
            checkBox.setEnabled(chkbox_enabled_list[i])
            checkBox.stateChanged.connect(lambda state, k=key:
                                          parent.chkbox_change(state, fld_ape,
                                                               obj_img, k))

        self.setLayout(layout)

    def value_change(self, key):
        label, lineEdit, checkBox = self.dlog_attrs[key]
#        print(" value_change: {} {} '{}'".format(str(self.item_inputs),
#              key, lineEdit.text()))
        try:
            value = float(lineEdit.text())
        except ValueError:
            return
        else:
            if (key in self.item_inputs or len(self.item_inputs) < 2):
                self.item_inputs[key] = value
                if not checkBox.isChecked():
                    checkBox.setCheckState(qt.Checked)
            self.update_values()

    def update_values(self, cell):
        """ update the display for the etendue cell being updated """
        print('sgb.update_values:', cell)
        for key in self.keys:
            value = cell[key]
            label, lineEdit, checkBox = self.dlog_attrs[key]
            lineEdit.setText(value_to_text(value))

        self.update_checkboxes()

    def chkbox_change(self, state, key):
        label, lineEdit, checkBox = self.dlog_attrs[key]
        checked = state == qt.Checked
#        print("chkbox_change:", self.item_inputs, key, checked)
        if checked:
            try:
                value = float(lineEdit.text())
            except ValueError:
                value = None
            finally:
                self.item_inputs[key] = value
                if not checkBox.isChecked():
                    checkBox.setChecked(True)
        else:
            if key in self.item_inputs:
                del self.item_inputs[key]
                if checkBox.isChecked():
                    checkBox.setChecked(False)
        self.update_checkboxes()
#        print("  chkbox_exit:", self.item_inputs, key, checked)

    def get_enabled_lists(self):
        if self.enabled_list is not None:
            enabled_list = self.enabled_list
        else:
            if len(self.item_inputs) < 2:
                enabled_list = [True]*7
            else:
                enabled_list = [False]*7
            for key in self.item_inputs:
                enabled_list[self.keys.index(key)] = True

        if self.chkbox_enabled_list is not None:
            chkbox_enabled_list = self.chkbox_enabled_list
        else:
            chkbox_enabled_list = enabled_list

        return enabled_list, chkbox_enabled_list

    def update_checkboxes(self):
        pass
#        enabled_list, chkbox_enabled_list = self.get_enabled_lists()
#
#        for i, key in enumerate(self.keys):
#            label, lineEdit, checkBox = self.dlog_attrs[key]
#
#            lineEdit.setEnabled(enabled_list[i])
#            checkBox.setEnabled(chkbox_enabled_list[i])


if __name__ == '__main__':

    import sys

    app = QApplication(sys.argv)
    dialog = IdealImagerDialog('infinite')
    dialog.exec()
