#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

.. Created on Tue Nov 27 21:26:08 2018

.. codeauthor: Michael J. Hayford
"""

from collections import namedtuple

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtWidgets import (QCheckBox, QComboBox, QAction, QLineEdit,
                             QDockWidget, QFormLayout, QWidget)


PanelInfo = namedtuple('PanelInfo', ['dock', 'panel_widget', 'menu_action'])

panels = {}


def create_dock_windows(gui_app):
    panels['wavelengths'] = create_dock_widget(gui_app,
                                               'wavelengths', 'Wavelengths',
                                               SpectrumWavelengthsPanel, True)
    panels['aperture'] = create_dock_widget(gui_app,
                                            'aperture', 'Aperture',
                                            AperturePanel, True)
    panels['field'] = create_dock_widget(gui_app,
                                         'field', 'Field of View',
                                         FieldOfViewPanel, True)
    panels['defocus'] = create_dock_widget(gui_app,
                                           'defocus', 'Focus Range',
                                           FocusRangePanel, True)
    panels['system'] = create_dock_widget(gui_app,
                                          'system', 'System Info',
                                          SystemSpecPanel, True)


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


class ModelBinding():
    """
    ModelBinding the the base class for binding part of the optical model
    to a UI element. UI elements should extend this class.

    When more getters/setters are needed, overwrite the get/set functions
    to directly get/set the model part
    """
    def __init__(self, gui_app, get_parent, field):
        """
        Create a new ModelBinding. gui_app is the instance that should be
        refreshed on change, get_parent is a function which returns the parent
        of the optical model element being bound, and field is the name of
        element (or list index if get_parent() returns a list)

        E.g: if the model part being changed was foo.bar.baz, get_parent should
        return foo.bar and field should be the string 'baz'.

        get_parent() is a function to ensure the current optical model is
        updated. If the model itself were passed in and then a new model was
        loaded the binding would become stale. Using a function prevents this.
        """
        self.gui_app = gui_app
        self.get_parent = get_parent
        self.field = field

    def set(self, value):
        """ Updates the model with the new value """
        toset = self.get_parent()
        if isinstance(toset, list):
            toset.__setitem__(self.field, value)
        else:
            setattr(toset, self.field, value)

    def get(self):
        """ Retreives the model's current value """
        toget = self.get_parent()
        if isinstance(toget, list):
            return toget.__getitem__(self.field)
        else:
            return getattr(toget, self.field)


class EnumChoiceWidget(ModelBinding):
    def __init__(self, gui_app, get_parent, field, combo_items):
        super().__init__(gui_app, get_parent, field)
        w = QComboBox()
        for item in combo_items:
            w.addItem(item)
        w.currentIndexChanged.connect(self.currentIndexChanged)
        self.widget = w

    def currentIndexChanged(self):
        self.set(self.widget.currentIndex())
        self.gui_app.refresh_gui()

    def refresh(self):
        self.widget.setCurrentIndex(self.get())


class ListChoiceWidget(ModelBinding):
    def __init__(self, gui_app, get_parent, field, combo_items):
        super().__init__(gui_app, get_parent, field)
        w = QComboBox()
        for item in combo_items:
            w.addItem(item)
        w.currentIndexChanged.connect(self.currentIndexChanged)
        self.widget = w
        self.combo_items = combo_items

    def currentIndexChanged(self):
        self.set(self.widget.currentIndex())
        self.gui_app.refresh_gui()

    def refresh(self):
        idx = self.combo_items.index(self.get())
        self.widget.setCurrentIndex(idx)


class TextFieldWidget(ModelBinding):
    def __init__(self, gui_app, get_parent, field, valueFormat='{:s}'):
        super().__init__(gui_app, get_parent, field)
        w = QLineEdit()
        w.setAlignment(Qt.AlignLeft)
        w.editingFinished.connect(self.editingFinished)
        self.widget = w
        # valueFormat is how the data from the model is rendered in the textbox
        self.valueFormat = valueFormat
        # convert is called on the text to convert to the type the model uses
        self.convert = lambda a: a

    def editingFinished(self):
        self.set(self.convert(self.widget.text()))
        self.gui_app.refresh_gui()

    def refresh(self):
        self.widget.setText(self.valueFormat.format(self.get()))


class FloatFieldWidget(TextFieldWidget):
    """ FloatFieldWidget is like a TextFieldWidget but only for floats """
    def __init__(self, gui_app, root_fn, field, valueformat='{:.7g}'):
        super().__init__(gui_app, root_fn, field, valueformat)
        self.convert = float
        self.widget.setValidator(QDoubleValidator())


class SpectrumWavelengthsPanel(QWidget):
    def __init__(self, gui_app, parent=None):
        super().__init__(parent)

        self.gui_app = gui_app
        self.ctrl_wvl_edit = FloatFieldWidget(gui_app, self.root, 
                                              'central_wvl')
        self.red_wvl_edit = FloatFieldWidget(gui_app, 
                                             lambda: self.root().wavelengths, 
                                             -1)
        self.blue_wvl_edit = FloatFieldWidget(gui_app, 
                                              lambda: self.root().wavelengths, 
                                              0)
        wavlnsLayout = QFormLayout()
        wavlnsLayout.addRow('central', self.ctrl_wvl_edit.widget)
        wavlnsLayout.addRow('red', self.red_wvl_edit.widget)
        wavlnsLayout.addRow('blue', self.blue_wvl_edit.widget)

        self.setLayout(wavlnsLayout)

    def root(self):
        return self.gui_app.app_manager.model['optical_spec']['wvls']

    def update(self, opt_model):
        """ push backend data to widgets """
        num_wvls = len(opt_model['optical_spec']['wvls'].wavelengths)
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
    obj_img_items = ["object", "image"]
    pupil_type_items = ["epd", "NA", "f/#"]

    def __init__(self, gui_app, parent=None):
        super().__init__(parent)

        self.gui_app = gui_app
        apertureLayout = QFormLayout()

        self.obj_img_combo = ListChoiceWidget(gui_app, lambda: self.root(),
                                              None, self.obj_img_items)
        self.obj_img_combo.get = lambda: self.root().key[0]
        self.obj_img_combo.set = lambda value: \
            self.root().mutate_pupil_type((self.obj_img_items[value], 
                                           self.root().key[1]))
        apertureLayout.addRow('space', self.obj_img_combo.widget)

        self.pupil_type_combo = ListChoiceWidget(gui_app, self.root, None, 
                                                 self.pupil_type_items)
        self.pupil_type_combo.get = lambda: self.root().key[1]
        self.pupil_type_combo.set = lambda value: \
            self.root().mutate_pupil_type((self.root().key[0], 
                                           self.pupil_type_items[value]))
        apertureLayout.addRow('type', self.pupil_type_combo.widget)

        self.aperture_edit = FloatFieldWidget(gui_app, self.root, 'value')
        apertureLayout.addRow('value', self.aperture_edit.widget)

        self.setLayout(apertureLayout)

    def root(self):
        return self.gui_app.app_manager.model['optical_spec']['pupil']

    def update(self, opt_model):
        """ push backend data to widgets """
        self.obj_img_combo.refresh()
        self.pupil_type_combo.refresh()
        self.aperture_edit.refresh()


class FieldOfViewPanel(QWidget):
    obj_img_items = ["object", "image"]
    field_type_items = ["height", "angle", "real height"]

    def __init__(self, gui_app, parent=None):
        super().__init__(parent)

        self.gui_app = gui_app

        fieldLayout = QFormLayout()

        self.obj_img_combo = ListChoiceWidget(gui_app, lambda: self.root(),
                                              None, self.obj_img_items)
        self.obj_img_combo.get = lambda: self.root().key[0]
        self.obj_img_combo.set = lambda value: \
            self.root().mutate_field_type((self.obj_img_items[value], 
                                           self.root().key[1]))
        fieldLayout.addRow('space', self.obj_img_combo.widget)

        self.field_combo = ListChoiceWidget(gui_app, lambda: self.root(),
                                              None, self.field_type_items)
        self.field_combo.get = lambda: self.root().key[1]
        self.field_combo.set = lambda value: \
            self.root().mutate_field_type((self.root().key[0], 
                                           self.field_type_items[value]))
        fieldLayout.addRow('type', self.field_combo.widget)

        self.field_edit = FloatFieldWidget(gui_app, self.root, 'value')
        fieldLayout.addRow('value', self.field_edit.widget)

        self.setLayout(fieldLayout)

    def root(self):
        return self.gui_app.app_manager.model['optical_spec']['fov']

    def update(self, opt_model):
        """ push backend data to widgets """
        self.obj_img_combo.refresh()
        self.field_combo.refresh()
        self.field_edit.refresh()


class FocusRangePanel(QWidget):
    def __init__(self, gui_app, parent=None):
        super().__init__(parent)

        self.gui_app = gui_app

        self.focus_shift_edit = FloatFieldWidget(gui_app, self.root, 
                                                 'focus_shift')
        self.defocus_range_edit = FloatFieldWidget(gui_app, self.root,
                                                   'defocus_range')

        focusRangeLayout = QFormLayout()
        focusRangeLayout.addRow('focus shift', self.focus_shift_edit.widget)
        focusRangeLayout.addRow('defocus range', self.defocus_range_edit.widget)

        self.setLayout(focusRangeLayout)

    def root(self):
        return self.gui_app.app_manager.model['optical_spec']['focus']

    def update(self, opt_model):
        """ push backend data to widgets """
        self.focus_shift_edit.refresh()
        self.defocus_range_edit.refresh()


class SystemSpecPanel(QWidget):
    def __init__(self, gui_app, parent=None):
        super().__init__(parent)

        self.dim_combo_items = ["mm", "cm", "meters", "inches", "feet"]

        self.gui_app = gui_app

        systemLayout = QFormLayout()

        self.dimension_combo = ListChoiceWidget(gui_app, lambda: self.root(),
                                                'dimensions', self.dim_combo_items)

        def combo_set(i):
            setattr(self.root(), 'dimensions', self.dim_combo_items[i])
        self.dimension_combo.set = combo_set
        # self.dimension_combo.set = lambda value: setattr(self.root(), 'dimensions', value)

        self.title_edit = TextFieldWidget(gui_app, self.root, 'title')
        systemLayout.addRow('title', self.title_edit.widget)

        systemLayout.addRow('units', self.dimension_combo.widget)

        self.setLayout(systemLayout)

    def root(self):
        return self.gui_app.app_manager.model['system_spec']

    def update(self, opt_model):
        """ push backend data to widgets """
        self.dimension_combo.refresh()
        self.title_edit.refresh()
