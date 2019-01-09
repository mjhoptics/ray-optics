#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Support for fully featured QT windows for plotting/matplotlib

.. Created on Wed Nov  7 15:04:19 2018

.. codeauthor: Michael J. Hayford
"""
from PyQt5.QtCore import Qt as qt
from PyQt5.QtWidgets import (QHBoxLayout, QVBoxLayout, QWidget, QLineEdit,
                             QRadioButton, QGroupBox, QSizePolicy,
                             QListWidget, QListWidgetItem)

from matplotlib.backends.backend_qt5agg \
     import (NavigationToolbar2QT as NavigationToolbar)

from rayoptics.gui.appmanager import ModelInfo
#from rayoptics.gui.appcmds import update_figure_view
from rayoptics.qtgui.plotcanvas import PlotCanvas
from rayoptics.mpl.axisarrayfigure import Fit


def update_figure_view(plotFigure):
    plotFigure.refresh()


class CommandItem(QListWidgetItem):
    def __init__(self, parent, txt, cntxt):
        super().__init__(parent)
        self.setData(qt.DisplayRole, txt)
        self.setData(qt.EditRole, cntxt)

    def data(self, role):
        if role == qt.DisplayRole:
            return self.txt
        elif role == qt.EditRole:
            return self.cntxt
        else:
            return None

    def setData(self, role, value):
        if role == qt.DisplayRole:
            self.txt = value
            return True
        elif role == qt.EditRole:
            self.cntxt = value
            return True
        else:
            return False


def create_cmdplot_view(app, fig, title, view_width, view_ht,
                        commands=None, add_scale_panel=False,
                        add_nav_toolbar=False):
    # construct the top level widget
    widget = QWidget()

    top_layout = QHBoxLayout(widget)
    # set the layout on the widget
    widget.setLayout(top_layout)

    command_panel = create_command_panel(commands)
    top_layout.addWidget(command_panel)

    pc = PlotCanvas(app, fig)
    top_layout.addWidget(pc)

    mi = ModelInfo(app.app_manager.model, update_figure_view, (fig,))
    sub_window = app.add_subwindow(widget, mi)
    sub_window.setWindowTitle(title)
    orig_x, orig_y = app.initial_window_offset()
    sub_window.setGeometry(orig_x, orig_y, view_width, view_ht)

    sub_window.show()


def create_command_panel(commands):
    command_panel = QListWidget()

    for c in commands:
        cmd_txt, cntxt = c
        btn = CommandItem(command_panel, cmd_txt, cntxt)

    command_panel.setCurrentRow(0)
    command_panel.itemClicked.connect(on_command)
    width = command_panel.size()
    hint = command_panel.sizeHint()
    frame_width = command_panel.frameWidth() + 2
    column_width = command_panel.sizeHintForColumn(0) + 2*frame_width
    command_panel.setMinimumWidth(column_width)
    command_panel.setMaximumWidth(column_width)
    command_panel.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Minimum)

    return command_panel


def on_command(item):
    cntxt = item.data(qt.EditRole)
    fct, args, kwargs = cntxt
    fct(*args, **kwargs)

#        figure.plot()


def create_plot_view(app, fig, title, view_width, view_ht,
                     add_scale_panel=True, add_nav_toolbar=False):
    pc = PlotCanvas(app, fig)

    # construct the top level widget
    widget = QWidget()
    # construct the top level layout
    layout = QVBoxLayout(widget)

    # set the layout on the widget
    widget.setLayout(layout)

    if add_scale_panel:
        psp = create_plot_scale_panel(app, pc)
        layout.addWidget(psp)

    mi = ModelInfo(app.app_manager.model, update_figure_view, (fig,))
    sub_window = app.add_subwindow(widget, mi)
    sub_window.setWindowTitle(title)
    orig_x, orig_y = app.initial_window_offset()
    sub_window.setGeometry(orig_x, orig_y, view_width, view_ht)

    layout.addWidget(pc)

    if add_nav_toolbar:
        layout.addWidget(NavigationToolbar(pc, sub_window))

    sub_window.show()


def create_plot_scale_panel(app, pc):
    groupBox = QGroupBox("Plot Scale", app)

    user_scale_wdgt = QLineEdit()
    user_scale_wdgt.setReadOnly(True)
    pf = pc.figure
    cntxt = pf, user_scale_wdgt
    user_scale_wdgt.editingFinished.connect(lambda:
                                            on_plot_scale_changed(cntxt))
    fit_all_btn = QRadioButton("Fit All")
    fit_all_btn.setChecked(pf.scale_type == Fit.All)
    fit_all_btn.toggled.connect(lambda:
                                on_plot_scale_toggled(cntxt, Fit.All))
    fit_all_same_btn = QRadioButton("Fit All Same")
    fit_all_same_btn.setChecked(pf.scale_type == Fit.All_Same)
    fit_all_same_btn.toggled.connect(lambda: on_plot_scale_toggled(
                                             cntxt, Fit.All_Same))
    user_scale_btn = QRadioButton("User Scale")
    user_scale_btn.setChecked(pf.scale_type == Fit.User_Scale)
    user_scale_btn.toggled.connect(lambda: on_plot_scale_toggled(
                                           cntxt, Fit.User_Scale))
    box = QHBoxLayout()
    box.addWidget(fit_all_btn)
    box.addWidget(fit_all_same_btn)
    box.addWidget(user_scale_btn)
    box.addWidget(user_scale_wdgt)

    groupBox.setLayout(box)

    return groupBox


def on_plot_scale_toggled(cntxt, scale_type):
    plotFigure, scale_wdgt = cntxt
    plotFigure.scale_type = scale_type
    if scale_type == Fit.User_Scale:
        scale_wdgt.setReadOnly(False)
        scale_wdgt.setText('{:7.4f}'.format(plotFigure.user_scale_value))
    else:
        scale_wdgt.setReadOnly(True)

    plotFigure.plot()


def on_plot_scale_changed(cntxt):
    plotFigure, scale_wdgt = cntxt
    eval_str = scale_wdgt.text()
    try:
        val = eval(eval_str)
        plotFigure.user_scale_value = val
        scale_wdgt.setText('{:7.4f}'.format(val))
    except IndexError:
        return ''

    plotFigure.plot()
