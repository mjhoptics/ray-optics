#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
""" Qt5 dialog box for ideal imager ui

.. Created on Tue Jun 18 13:09:19 2019

.. codeauthor: Michael J. Hayford
"""

import math

from rayoptics.optical import specsheet

from PyQt5.QtCore import Qt as qt
from PyQt5.QtWidgets import (QApplication, QDialog, QRadioButton,
                             QStackedWidget, QFormLayout, QGridLayout,
                             QGroupBox, QHBoxLayout, QLabel, QLineEdit,
                             QPushButton, QCheckBox, QVBoxLayout)


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

        # setup infinite conjugate defaults
        enabled_list = [False, False, False, False, True]
        chkbox_enabled_list = [False]*5
        imager_inputs = {'s': -math.inf, 'f': None}
        imager = specsheet.IdealImager(None, -math.inf, None, None, None)
        iigb = ImagerSpecGroupBox(enabled_list=enabled_list, imager=imager,
                                  chkbox_enabled_list=chkbox_enabled_list,
                                  imager_inputs=imager_inputs)
        self.imager_stack['infinite'] = iigb

        self.imager_groupbox_stack = QStackedWidget()
        self.imager_groupbox_stack.addWidget(iigb)
        self.imager_groupbox_stack.addWidget(figb)

        imager_groupbox = self.imager_stack[self.conjugate_type]
        self.imager_groupbox_stack.setCurrentWidget(imager_groupbox)

        mainLayout = QHBoxLayout()
        mainLayout.addWidget(self.conjugate_box)
        mainLayout.addWidget(self.imager_groupbox_stack)
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
        prev = self.imager_groupbox_stack.currentWidget()
        self.conjugate_type = conj_type
        new = self.imager_stack[conj_type]
        if prev.imager is not None:
            if new.imager is None:
                new.imager = prev.imager
            else:
                prev_enabled_list, _ = prev.get_enabled_lists()
                new_enabled_list, _ = new.get_enabled_lists()
                new_imager = list(new.imager)
                for i, p in enumerate(prev.imager):
#                    print(i, p, prev_enabled_list[i], new_enabled_list[i])
                    # only transfer values if both items are enabled
                    if (prev_enabled_list[i] and new_enabled_list[i]):
                        new_imager[i] = p
                        key = prev.keys[i]
                        if key in prev.imager_inputs:
                            new.chkbox_change(qt.Checked, key)
                            new.imager_inputs[key] = p
                new.imager = specsheet.IdealImager(*new_imager)
        self.imager_groupbox_stack.setCurrentWidget(new)
        new.update_values(new_imager)


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
            lineEdit.setText(ImagerSpecGroupBox.value_to_text(imager[i]))
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
                lineEdit.setText(ImagerSpecGroupBox.value_to_text(value))

        self.update_checkboxes()

    def value_to_text(value):
        if value is None:
            value_text = ''
        elif type(value) is str:
            value_text = value
        else:
            value_text = "{:> #.5f}".format(value)
        return value_text

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


if __name__ == '__main__':

    import sys

    app = QApplication(sys.argv)
    dialog = IdealImagerDialog('infinite')
    dialog.exec()
