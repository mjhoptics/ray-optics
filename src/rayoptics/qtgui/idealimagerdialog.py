#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2019 Michael J. Hayford
""" Qt5 dialog box for ideal imager ui

.. Created on Tue Jun 18 13:09:19 2019

.. codeauthor: Michael J. Hayford
"""

from rayoptics.parax.specsheet import create_specsheets
from rayoptics.parax import idealimager
from rayoptics.parax.etendue import (obj_img_set, fld_ape_set,
                                     fld_labels, ap_labels)

from rayoptics.util.dict2d import dict2D

from PyQt5.QtCore import Qt as qt
from PyQt5.QtWidgets import (QApplication, QDialog, QRadioButton, QWidget,
                             QStackedWidget, QGridLayout, QGroupBox,
                             QHBoxLayout, QLabel, QLineEdit, QCheckBox,
                             QVBoxLayout, QDialogButtonBox, QMessageBox)


def value_to_text(value, fmt_str="{:> #.5f}"):
    if value is None:
        value_text = ''
    elif type(value) is str:
        value_text = value
    else:
        value_text = fmt_str.format(value)
    return value_text


class IdealImagerDialog(QWidget):
    def __init__(self, conjugate_type, specsheets, cmd_fct=None,
                 **kwargs):
        super().__init__(**kwargs)

        self.dlog_attrs = {}

        self.conjugate_type = conjugate_type
        self.conjugate_box = self.createConjugateBox(itype=self.conjugate_type)

        self.specsheet_dict = specsheets

        self.imager_stack = {}
        self.etendue_stack = {}
        self.imager_groupbox_stack = QStackedWidget()
        self.etendue_groupbox_stack = QStackedWidget()

        for key, specsheet in specsheets.items():
            isgb = ImagerSpecGroupBox(self, specsheet)
            self.imager_stack[key] = isgb
            self.imager_groupbox_stack.addWidget(isgb)
            egb = EtendueGroupBox(self, key, specsheet)
            self.etendue_stack[key] = egb
            self.etendue_groupbox_stack.addWidget(egb)

        imager_groupbox = self.imager_stack[self.conjugate_type]
        self.imager_groupbox_stack.setCurrentWidget(imager_groupbox)

        etendue_groupbox = self.etendue_stack[self.conjugate_type]
        self.etendue_groupbox_stack.setCurrentWidget(etendue_groupbox)

        # construct the top level layout
        overallLayout = QVBoxLayout(self)

        mainLayout = QHBoxLayout()
        mainLayout.addWidget(self.conjugate_box)
        mainLayout.addWidget(self.imager_groupbox_stack)
        mainLayout.addWidget(self.etendue_groupbox_stack)

        self.button_box = self.createButtonBox(cmd_fct)
        overallLayout.addLayout(mainLayout)
        overallLayout.addWidget(self.button_box)

        self.setLayout(overallLayout)

        self.setWindowTitle("Optical Spec Sheet")

    def createButtonBox(self, cmd_fct):
        def clicked(button):
            command = button.text()
            specsheet = self.specsheet_dict[self.conjugate_type]
            if cmd_fct:
                try:
                    cmd_fct(self, command, specsheet)
                except Exception as e:
                    print(str(e))
                    QMessageBox.warning(self,
                                        self.tr("Ray-Optics"), 
                                        self.tr("Please provide correct inputs."))
            else:
                print(button.text(), 'button pressed')

        buttonbox = QDialogButtonBox(qt.Horizontal, self)
        buttonbox.addButton('New', QDialogButtonBox.ApplyRole)
        buttonbox.addButton(QDialogButtonBox.Apply)
        buttonbox.addButton('Update', QDialogButtonBox.ApplyRole)
        buttonbox.addButton(QDialogButtonBox.Close)
        for b in buttonbox.buttons():
            b.setAutoDefault(False)
#        buttonbox.setCenterButtons(True)
        buttonbox.clicked.connect(clicked)
        return buttonbox

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
        self.update_conjugate(conj_type)

    def update_conjugate(self, conj_type):
        prev_specsheet = self.specsheet_dict[self.conjugate_type]
        self.conjugate_type = conj_type

        new_specsheet = self.specsheet_dict[conj_type]
        new_partition, max_inputs = new_specsheet.partition_defined()
        if new_partition is None:
            if prev_specsheet.imager.f is not None:
                prev_f = prev_specsheet.imager.f
                new_specsheet.imager_inputs['f'] = prev_f

        new_imager_groupbox = self.imager_stack[conj_type]
        self.imager_groupbox_stack.setCurrentWidget(new_imager_groupbox)

        etendue_groupbox = self.etendue_stack[conj_type]
        self.etendue_groupbox_stack.setCurrentWidget(etendue_groupbox)

        self.update_values()

    def update_values(self):
        """ callback routine for any dialog value change """
        conj_type = self.conjugate_type
        specsheet = self.specsheet_dict[conj_type]

        imager_groupbox = self.imager_stack[conj_type]
        etendue_groupbox = self.etendue_stack[conj_type]

        imager_inputs = imager_groupbox.specsheet.imager_inputs
        etendue_inputs = etendue_groupbox.specsheet.etendue_inputs
        imager, etendue_values = specsheet.generate_from_inputs(imager_inputs,
                                                                etendue_inputs)

        specsheet.imager_inputs = imager_inputs
        specsheet.etendue_inputs = etendue_inputs

        imager_groupbox.update_values()
        etendue_groupbox.update_values()

    def update_checkboxes(self):
        """ callback routine for any dialog checkbox change """
        imager_groupbox = self.imager_stack[self.conjugate_type]
        etendue_groupbox = self.etendue_stack[self.conjugate_type]

        imager_groupbox.update_checkboxes()
        etendue_groupbox.update_checkboxes()


class ImagerSpecGroupBox(QGroupBox):
    def __init__(self, parent, specsheet, keys=None, labels=None, **kwargs):
        super().__init__(title='Imager specs', **kwargs)

        self.parent = parent
        self.specsheet = specsheet

        self.dlog_attrs = {}

        self.keys = keys if keys else idealimager.ideal_imager_keys
        self.labels = labels if labels else idealimager.ideal_imager_labels

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
            lineEdit.setText(value_to_text(specsheet.imager[i]))
            # Use returnPressed signal rather than editingFinished. The latter
            #  fires on change of focus as well as return. We only want a
            #  signal when the lineEdit contents change and are committed by
            #  the user.
            lineEdit.returnPressed.connect(lambda k=key: self.value_change(k))

            if key in specsheet.imager_inputs:
                if not checkBox.isChecked():
                    checkBox.setCheckState(qt.Checked)
            checkBox.stateChanged.connect(lambda state, k=key:
                                          self.chkbox_change(state, k))

        self.setLayout(layout)
        self.update_values()

    def value_change(self, imager_key):
        label, lineEdit, checkBox = self.dlog_attrs[imager_key]

        try:
            value = float(lineEdit.text())
        except ValueError:
            value = None
        finally:
            imager_inputs = self.specsheet.imager_inputs
            if value is not None:
                if imager_key in imager_inputs or len(imager_inputs) < 2:
                    imager_inputs[imager_key] = value
                    if not checkBox.isChecked():
                        checkBox.setCheckState(qt.Checked)
            else:
                if imager_key in imager_inputs:
                    del imager_inputs[imager_key]
                    if checkBox.isChecked():
                        checkBox.setChecked(False)
            self.parent.update_values()

    def chkbox_change(self, state, imager_key):
        label, lineEdit, checkBox = self.dlog_attrs[imager_key]

        checked = state == qt.Checked
        if checked:
            try:
                value = float(lineEdit.text())
            except ValueError:
                value = None
            finally:
                self.specsheet.imager_inputs[imager_key] = value
                if not checkBox.isChecked():
                    checkBox.setChecked(True)
        else:
            if imager_key in self.specsheet.imager_inputs:
                del self.specsheet.imager_inputs[imager_key]
                if checkBox.isChecked():
                    checkBox.setChecked(False)
        self.parent.update_checkboxes()

    def update_values(self):
        if self.specsheet.imager is not None:
            value_list = self.specsheet.imager
            # refill gui widgets and set enabled state
            for i, key in enumerate(self.keys):
                value = value_list[i]
                label, lineEdit, checkBox = self.dlog_attrs[key]
                lineEdit.setText(value_to_text(value))

        self.update_checkboxes()

    def update_checkboxes(self):
        partition, max_inputs = self.specsheet.partition_defined()
        if partition == 'imager':
            # the inputs for imager are sufficient to calculate all; mark
            # inputs editable and disable the others
            enabled_list = [False]*5
            for key in self.specsheet.imager_inputs:
                enabled_list[self.keys.index(key)] = True
        elif bool(partition):
            # an imager is defined by other partitions, this means:
            #  1) the imager key is not editable
            #  2a) if one imager input, the remainder are uneditable
            #  2b) if no imager inputs, all but the imager key is editable
            if len(self.specsheet.imager_inputs):
                enabled_list = [False]*5
                for key in self.specsheet.imager_inputs:
                    enabled_list[self.keys.index(key)] = True
                imager_key = self.specsheet.imager_defined()
                if imager_key:
                    enabled_list[self.keys.index(imager_key)] = False
            else:
                enabled_list = [True]*5
                imager_key = self.specsheet.imager_defined()
                if imager_key:
                    enabled_list[self.keys.index(imager_key)] = False
        else:
            # no imager defined, all items are editable
            enabled_list = [True]*5

        # update the ui from the enabled list and inputs
        for i, key in enumerate(self.keys):
            label, lineEdit, checkBox = self.dlog_attrs[key]

            if key in self.specsheet.imager_inputs:
                if not checkBox.isChecked():
                    checkBox.setCheckState(qt.Checked)
            else:
                if checkBox.isChecked():
                    checkBox.setCheckState(qt.Unchecked)

            if self.specsheet.frozen_imager_inputs[i]:
                enabled_list[i] = False

            lineEdit.setEnabled(enabled_list[i])
            checkBox.setEnabled(enabled_list[i])


class EtendueGroupBox(QGroupBox):
    def __init__(self, parent, itype, specsheet, **kwargs):
        super().__init__(title='etendue definition', **kwargs)

        self.parent = parent
        self.specsheet = specsheet

        self.groupboxes = dict2D(fld_ape_set, obj_img_set)

        if itype == 'infinite':
            gb00 = SpaceGroupBox(parent=self, title='object space',
                                 obj_img='object', fld_ape='field',
                                 keys=(fld_labels[1],))
            self.groupboxes['field']['object'] = gb00
            gb10 = SpaceGroupBox(parent=self, title='',
                                 obj_img='object', fld_ape='aperture',
                                 keys=(ap_labels[0],))
            self.groupboxes['aperture']['object'] = gb10
            gb01 = SpaceGroupBox(parent=self, title='image space',
                                 obj_img='image', fld_ape='field',
                                 keys=(fld_labels[0],))
            self.groupboxes['field']['image'] = gb01
            gb11 = SpaceGroupBox(parent=self, title='',
                                 obj_img='image', fld_ape='aperture',
                                 keys=(ap_labels[2], ap_labels[1]))
            self.groupboxes['aperture']['image'] = gb11

        else:
            gb00 = SpaceGroupBox(parent=self, title='object space',
                                 fld_ape='field', obj_img='object',
                                 keys=(fld_labels[0],))
            self.groupboxes['field']['object'] = gb00
            gb10 = SpaceGroupBox(parent=self, title='',
                                 fld_ape='aperture', obj_img='object',
                                 keys=(ap_labels[2], ap_labels[1]))
            self.groupboxes['aperture']['object'] = gb10
            gb01 = SpaceGroupBox(parent=self, title='image space',
                                 fld_ape='field', obj_img='image',
                                 keys=(fld_labels[0],))
            self.groupboxes['field']['image'] = gb01
            gb11 = SpaceGroupBox(parent=self, title='',
                                 fld_ape='aperture', obj_img='image',
                                 keys=(ap_labels[2], ap_labels[1]))
            self.groupboxes['aperture']['image'] = gb11

        layout = QGridLayout()
        layout.addWidget(gb00, 0, 0)
        layout.addWidget(gb10, 1, 0)
        layout.addWidget(gb01, 0, 1)
        layout.addWidget(gb11, 1, 1)

        self.setLayout(layout)

    def value_change(self, fld_ape, obj_img, key):
        """ callback routine for item value widget """
        inputs = self.specsheet.etendue_inputs[fld_ape][obj_img]
        gb = self.groupboxes[fld_ape][obj_img]
        label, lineEdit, checkBox = gb.dlog_attrs[key]
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
        inputs = self.specsheet.etendue_inputs[fld_ape][obj_img]
        gb = self.groupboxes[fld_ape][obj_img]
        label, lineEdit, checkBox = gb.dlog_attrs[key]
        checked = state == qt.Checked
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
        self.parent.update_checkboxes()

    def update_values(self):
        etendue_inputs = self.specsheet.etendue_inputs
        etendue_values = self.specsheet.etendue_values

        for fld_ape_key, fld_ape_value in etendue_values.items():
            for obj_img_key, values in fld_ape_value.items():
                sgb = self.groupboxes[fld_ape_key][obj_img_key]
                inputs = etendue_inputs[fld_ape_key][obj_img_key]
                sgb.update_values(inputs, values)
        self.update_checkboxes()

    def update_checkboxes(self):
        """ update the enabled and checked state for the etendue groupbox """
        etendue_inputs = self.specsheet.etendue_inputs
        etendue_values = self.specsheet.etendue_values
        partition, max_inputs = self.specsheet.partition_defined()

        for fld_ape_key, fld_ape_value in etendue_values.items():
            num_fld_ape_inputs = self.specsheet.partitions[fld_ape_key]
            for obj_img_key, values in fld_ape_value.items():
                sgb = self.groupboxes[fld_ape_key][obj_img_key]
                inputs = etendue_inputs[fld_ape_key][obj_img_key]
                if num_fld_ape_inputs == 2:
                    partition_defined = True
                elif num_fld_ape_inputs == 1 and len(inputs) == 1:
                    partition_defined = True
                elif num_fld_ape_inputs == 1 and len(inputs) == 0:
                    partition_defined = True if bool(partition) else False
                else:  # num_fld_ape_inputs == 0 and len(inputs) == 0
                    partition_defined = False

                sgb.update_checkboxes(inputs, values,
                                      partition_defined=partition_defined)


class SpaceGroupBox(QGroupBox):
    def __init__(self, parent, title, fld_ape, obj_img, keys,
                 labels=None, **kwargs):
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

            # Use returnPressed signal rather than editingFinished. The latter
            #  fires on change of focus as well as return. We only want a
            #  signal when the lineEdit contents change and are committed by
            #  the user.
            lineEdit.returnPressed.connect(lambda k=key:
                                           parent.value_change(fld_ape,
                                                               obj_img, k))

            checkBox.stateChanged.connect(lambda state, k=key:
                                          parent.chkbox_change(state, fld_ape,
                                                               obj_img, k))

        self.setLayout(layout)

    def update_values(self, inputs, values):
        """ update the display for the etendue cell being updated """
        for key in self.keys:
            value = values[key]
            label, lineEdit, checkBox = self.dlog_attrs[key]
            lineEdit.setText(value_to_text(value))

    def update_checkboxes(self, inputs, values, partition_defined=False):
        """ update the display for the etendue cell being updated

        A partition is an aperture or field pair of object/image inputs.

        If it is defined, this means that all attrs can be supplied.
        In this case, the inputs will have editable values and (checked)
        checkboxes; the remaining attrs will have uneditable values and
        (unchecked) checkboxes.

        If the partition is not defined, all values and checkboxes will be
        editable, the input attrs, if any, will be checked.
        """
        if partition_defined:
            for key in self.keys:
                label, lineEdit, checkBox = self.dlog_attrs[key]

                if key in inputs:
                    lineEdit.setEnabled(True)
                    checkBox.setEnabled(True)
                    if not checkBox.isChecked():
                        checkBox.setCheckState(qt.Checked)
                else:
                    lineEdit.setEnabled(False)
                    checkBox.setEnabled(False)
                    if checkBox.isChecked():
                        checkBox.setCheckState(qt.Unchecked)

        else:
            for key in self.keys:
                label, lineEdit, checkBox = self.dlog_attrs[key]
                lineEdit.setEnabled(True)
                checkBox.setEnabled(True)

                if key in inputs:
                    if not checkBox.isChecked():
                        checkBox.setCheckState(qt.Checked)
                else:
                    if checkBox.isChecked():
                        checkBox.setCheckState(qt.Unchecked)


if __name__ == '__main__':

    import sys

    class Dialog2(QDialog):
        def __init__(self, parent, conjugate_type, specsheets):
            QDialog.__init__(self, parent)
#            self.setModal(0)
            iid = IdealImagerDialog(conjugate_type, specsheets, parent=self)
            iid.update_values()

        def exit(self):
            self.close()

    app = QApplication(sys.argv)
    specsheets = create_specsheets()
    dialog = Dialog2(None, 'infinite', specsheets)
    dialog.show()
#    dialog.exec()
