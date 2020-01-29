#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

.. Created on Tue Nov 27 21:26:08 2018

.. codeauthor: Michael J. Hayford
"""

import logging
from collections import namedtuple

from PyQt5 import QtCore
from PyQt5.QtCore import QDate, QSize, Qt
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import (QCheckBox, QComboBox, QDateTimeEdit, QAction,
                             QDialog, QGroupBox, QHBoxLayout, QLabel,
                             QLineEdit, QListView, QListWidget, QDockWidget,
                             QListWidgetItem, QPushButton, QSpinBox,
                             QStackedWidget, QVBoxLayout, QFormLayout, QWidget)

from rayoptics.optical.model_enums import PupilType
from rayoptics.optical.model_enums import FieldType
from rayoptics.optical.model_enums import DimensionType


PanelInfo = namedtuple('PanelInfo', ['dock', 'panel_widget', 'menu_action'])

panels = {}


def create_dock_windows(gui_app):
    panels['wavelengths'] = create_dock_widget(gui_app,
                                               'wavelengths', 'Wavelengths',
                                               SpectrumWavelengthsPanel, False)
    panels['aperture'] = create_dock_widget(gui_app,
                                            'aperture', 'Aperture',
                                            AperturePanel, False)
    panels['field'] = create_dock_widget(gui_app,
                                         'field', 'Field of View',
                                         FieldOfViewPanel, False)
    panels['system'] = create_dock_widget(gui_app,
                                          'system', 'System Info',
                                          SystemSpecPanel, False)


def create_dock_widget(gui_app, item_key, label, panel, state):
    dock = QDockWidget(item_key, gui_app)
    panel_widget = panel(gui_app, dock)
    dock.setWidget(panel_widget)
    menu_action = create_menu_action(gui_app, item_key, label, state)
    gui_app.addDockWidget(Qt.RightDockWidgetArea, dock)
    dock.setVisible(state)
    return PanelInfo(dock, panel_widget, menu_action)


def update_dock_windows(gui_app):
    for pi in panels.values():
        if pi.menu_action.isChecked():
            pi.panel_widget.update(gui_app.app_manager.model)


def create_menu_action(gui_app, item_key, label, state=False):
    menu_action = QAction(label, gui_app, checkable=True)
    menu_action.setChecked(state)
    menu_action.triggered.connect(lambda state:
                                  togglePanel(gui_app, state, item_key))
    return menu_action


def togglePanel(gui_app, state, item_key):
    panel_info = panels[item_key]
    panel_info.dock.setVisible(state)
    if state:
        panel_info.panel_widget.update(gui_app.app_manager.model)


class UIWidget():
    def __init__(self, gui_app, widgetLabel, widgetClass,
                 rootEvalStr, get_eval_str, set_eval_str=None):
        self.gui_app = gui_app
        self.rootEvalStr = rootEvalStr
        self.get_eval_str = get_eval_str
        if set_eval_str is not None:
            self.set_eval_str = set_eval_str
        else:
            self.set_eval_str = get_eval_str + '={}'
#            self.set_eval_str = get_eval_str + '="{:s}"'
        self.widgetLabel = widgetLabel
        self.widget = widgetClass()

    def get_root_object(self):
        root = self.gui_app.app_manager.model
        if len(self.rootEvalStr) == 0:
            return root
        else:
            root_eval_str = ('root' + self.rootEvalStr)
            try:
                root = eval(root_eval_str)
            except IndexError:
                return root
            else:
                return root

    def data(self):
        root = self.get_root_object()
        eval_str = 'root' + self.get_eval_str
        try:
            val = eval(eval_str)
            return val
        except IndexError:
            return ''
        except TypeError:
            print('Data type error: ', eval_str, val)
            return ''

    def setData(self, value):
        root = self.get_root_object()
        exec_str = ('root' + self.set_eval_str).format(value)
        try:
            exec(exec_str)
            return True
        except IndexError:
            return False
        except SyntaxError:
            logging.info('Syntax error: "%s"', value)
            return False

    def updateData(self):
        """ push widget data to backend """
        pass

    def refresh(self):
        """ push backend data to widget """
        pass


class TextFieldWidget(UIWidget):
    def __init__(self, gui_app, widgetLabel, widgetClass,
                 rootEvalStr, get_eval_str, valueFormat, set_eval_str=None):
        super().__init__(gui_app, widgetLabel, widgetClass, rootEvalStr,
                         get_eval_str, set_eval_str)
        self.valueFormat = valueFormat
        self.widget.editingFinished.connect(self.updateData)
        self.widget.setAlignment(QtCore.Qt.AlignLeft)

    def updateData(self):
        """ push widget data to backend """
        self.setData(self.widget.text())
        self.gui_app.refresh_gui()

    def refresh(self):
        """ push backend data to widget """
        valStr = self.valueFormat.format(self.data())
        self.widget.setText(valStr)


class ChoiceWidget(UIWidget):
    def __init__(self, gui_app, widgetLabel, widgetClass,
                 rootEvalStr, get_eval_str, combo_items, set_eval_str=None):
        super().__init__(gui_app, widgetLabel, widgetClass, rootEvalStr,
                         get_eval_str, set_eval_str)
        for item in combo_items:
            self.widget.addItem(item)
        self.widget.currentIndexChanged.connect(self.updateData)

    def updateData(self):
        """ push widget data to backend """
        self.setData(self.widget.currentIndex())
        self.gui_app.refresh_gui()

    def refresh(self):
        """ push backend data to widget """
        self.widget.setCurrentIndex(self.data())


class SpectrumWavelengthsPanel(QWidget):
    rootEvalStr = '.optical_spec.spectral_region'
    evalStr = '.central_wvl', '.wavelengths[0]', '.wavelengths[-1]'

    def __init__(self, gui_app, parent=None):
        super().__init__(parent)

        self.ctrl_wvl_edit = TextFieldWidget(
            gui_app, 'central', QLineEdit,
            SpectrumWavelengthsPanel.rootEvalStr,
            SpectrumWavelengthsPanel.evalStr[0],
            '{:7.1f}')

        self.red_wvl_edit = TextFieldWidget(
            gui_app, 'red', QLineEdit,
            SpectrumWavelengthsPanel.rootEvalStr,
            SpectrumWavelengthsPanel.evalStr[1],
            '{:7.1f}')

        self.blue_wvl_edit = TextFieldWidget(
            gui_app, 'blue', QLineEdit,
            SpectrumWavelengthsPanel.rootEvalStr,
            SpectrumWavelengthsPanel.evalStr[2],
            '{:7.1f}')

        wavlnsLayout = QFormLayout()
        wavlnsLayout.addRow('central', self.ctrl_wvl_edit.widget)
        wavlnsLayout.addRow('red', self.red_wvl_edit.widget)
        wavlnsLayout.addRow('blue', self.blue_wvl_edit.widget)

        self.achroCheckBox = QCheckBox("Achromatic")
        wavlnsLayout.addRow('', self.achroCheckBox)

        self.setLayout(wavlnsLayout)

    def update(self, opt_model):
        """ push backend data to widgets """
        num_wvls = len(opt_model.optical_spec.spectral_region.wavelengths)
        if num_wvls == 1:
            self.ctrl_wvl_edit.refresh()
            self.red_wvl_edit.widget.setText('')
            self.blue_wvl_edit.widget.setText('')
        elif num_wvls == 2:
            self.ctrl_wvl_edit.widget.setText('')
            self.red_wvl_edit.refresh()
            self.blue_wvl_edit.refresh()
        elif num_wvls > 2:
            self.ctrl_wvl_edit.refresh()
            self.red_wvl_edit.refresh()
            self.blue_wvl_edit.refresh()


class AperturePanel(QWidget):
    rootEvalStr = '.optical_spec.pupil'
    evalStr = '.value', '.get_pupil_type()'
    comboItems = ["Ent Pupil Diam", "Object NA", "F/#", "NA"]
    set_combo_str = '.mutate_pupil_type(PupilType({}))'

    def __init__(self, gui_app, parent=None):
        super().__init__(parent)

        apertureLayout = QFormLayout()

        self.aperture_combo = ChoiceWidget(gui_app, 'pupil_type', QComboBox,
                                           AperturePanel.rootEvalStr,
                                           AperturePanel.evalStr[1],
                                           AperturePanel.comboItems,
                                           set_eval_str=
                                           AperturePanel.set_combo_str)
        apertureLayout.addRow('Type', self.aperture_combo.widget)

        self.aperture_edit = TextFieldWidget(gui_app, 'value', QLineEdit,
                                             AperturePanel.rootEvalStr,
                                             AperturePanel.evalStr[0],
                                             '{:12.5f}')
        apertureLayout.addRow('value', self.aperture_edit.widget)

        self.setLayout(apertureLayout)

    def update(self, opt_model):
        """ push backend data to widgets """
        self.aperture_combo.refresh()
        self.aperture_edit.refresh()


class FieldOfViewPanel(QWidget):
    rootEvalStr = '.optical_spec.field_of_view'
    evalStr = '.value', '.get_field_type()'
    comboItems = ["Object Angle", "Object Height", "Image Height"]
    set_combo_str = '.mutate_field_type(FieldType({}))'

    def __init__(self, gui_app, parent=None):
        super().__init__(parent)

        fieldLayout = QFormLayout()

        self.field_combo = ChoiceWidget(gui_app, 'field_type', QComboBox,
                                        FieldOfViewPanel.rootEvalStr,
                                        FieldOfViewPanel.evalStr[1],
                                        FieldOfViewPanel.comboItems,
                                        set_eval_str=
                                        FieldOfViewPanel.set_combo_str)
        fieldLayout.addRow('Type', self.field_combo.widget)

        self.field_edit = TextFieldWidget(gui_app, 'value', QLineEdit,
                                          FieldOfViewPanel.rootEvalStr,
                                          FieldOfViewPanel.evalStr[0],
                                          '{:12.5f}')
        fieldLayout.addRow('value', self.field_edit.widget)

        self.setLayout(fieldLayout)

    def update(self, opt_model):
        """ push backend data to widgets """
        self.field_combo.refresh()
        self.field_edit.refresh()


class SystemSpecPanel(QWidget):
    rootEvalStr = '.system_spec'
    evalStr = ('.dimensions.value', '.title', '.initials', '.temperature',
               '.pressure')
    comboItems = ["mm", "cm", "m", "inches"]
    set_combo_str = '.dimensions=DimensionType({})'
    set_text_str = ".title='{:s}'", ".initials='{:s}'"

    def __init__(self, gui_app, parent=None):
        super().__init__(parent)

        systemLayout = QFormLayout()

        self.dimension_combo = ChoiceWidget(gui_app, 'dimensions', QComboBox,
                                            SystemSpecPanel.rootEvalStr,
                                            SystemSpecPanel.evalStr[0],
                                            SystemSpecPanel.comboItems,
                                            set_eval_str=
                                            SystemSpecPanel.set_combo_str)
        systemLayout.addRow('system units', self.dimension_combo.widget)

        self.title_edit = TextFieldWidget(gui_app, 'title', QLineEdit,
                                          SystemSpecPanel.rootEvalStr,
                                          SystemSpecPanel.evalStr[1],
                                          '{:s}',
                                          set_eval_str=
                                          SystemSpecPanel.set_text_str[0])
        systemLayout.addRow('title', self.title_edit.widget)

        self.initials_edit = TextFieldWidget(gui_app, 'initials', QLineEdit,
                                             SystemSpecPanel.rootEvalStr,
                                             SystemSpecPanel.evalStr[2],
                                             '{:s}',
                                             set_eval_str=
                                             SystemSpecPanel.set_text_str[1])
        systemLayout.addRow('initials', self.initials_edit.widget)

        self.temp_edit = TextFieldWidget(gui_app, 'temperature', QLineEdit,
                                         SystemSpecPanel.rootEvalStr,
                                         SystemSpecPanel.evalStr[3],
                                         '{:7.1f}')
        systemLayout.addRow('temperature', self.temp_edit.widget)

        self.pressure_edit = TextFieldWidget(gui_app, 'pressure', QLineEdit,
                                             SystemSpecPanel.rootEvalStr,
                                             SystemSpecPanel.evalStr[4],
                                             '{:9.3f}')
        systemLayout.addRow('pressure', self.pressure_edit.widget)

        self.setLayout(systemLayout)

    def update(self, opt_model):
        """ push backend data to widgets """
        self.dimension_combo.refresh()
        self.title_edit.refresh()
        self.initials_edit.refresh()
        self.temp_edit.refresh()
        self.pressure_edit.refresh()
