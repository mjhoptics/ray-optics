#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Support for fully featured QT windows for plotting/matplotlib

.. Created on Wed Nov  7 15:04:19 2018

.. codeauthor: Michael J. Hayford
"""
from pathlib import Path

from PyQt5.QtCore import Qt as qt
from PyQt5.QtCore import QSize
from PyQt5 import QtGui
from PyQt5.QtWidgets import (QHBoxLayout, QVBoxLayout, QWidget, QLineEdit,
                             QRadioButton, QGroupBox, QSizePolicy, QCheckBox,
                             QListWidget, QListWidgetItem, QToolBar, QMenu)

from matplotlib.backends.backend_qt5agg \
     import (NavigationToolbar2QT as NavigationToolbar)
from matplotlib.backends.backend_qt5agg \
     import FigureCanvasQTAgg as FigureCanvas

from rayoptics.gui.appmanager import ModelInfo
from rayoptics.mpl.axisarrayfigure import Fit
from rayoptics.mpl.styledfigure import StyledFigure

from opticalglass import glassmap as gm
from opticalglass import glassmapviewer as gmv


class PlotCanvas(FigureCanvas):

    def __init__(self, parent, figure, accept_drops=True, drop_action=None):
        FigureCanvas.__init__(self, figure)
        self.setParent(parent)
        self.setAcceptDrops(True)
        self.drop_action = drop_action if drop_action else NullDropAction()

        # Next 2 lines are needed so that key press events are correctly
        #  passed with mouse events
        # https://github.com/matplotlib/matplotlib/issues/707/
        self.setFocusPolicy(qt.ClickFocus)
        self.setFocus()

        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

        self.figure.plot()

    def dragEnterEvent(self, event):
        if event.mimeData().hasFormat("text/plain"):
            self.drop_action.dragEnterEvent(self, event)
            event.acceptProposedAction()

    def dragMoveEvent(self, event):
        if event.mimeData().hasFormat("text/plain"):
            self.drop_action.dragMoveEvent(self, event)
            event.acceptProposedAction()

    def dragLeaveEvent(self, event):
        self.drop_action.dragLeaveEvent(self, event)

    def dropEvent(self, event):
        if event.mimeData().hasText():
            if self.drop_action.dropEvent(self, event):
                event.acceptProposedAction()
        else:
            event.ignore()


class NullDropAction():
    def dragEnterEvent(self, view, event):
        pass

    def dragMoveEvent(self, view, event):
        pass

    def dragLeaveEvent(self, view, event):
        pass

    def dropEvent(self, view, event):
        return False


def update_figure_view(plotFigure, **kwargs):
    plotFigure.refresh(**kwargs)


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


def create_command_panel(fig, commands):
    command_panel = QListWidget()

    for c in commands:
        cmd_txt, cntxt = c
        cntxt[2]['figure'] = fig
        btn = CommandItem(command_panel, cmd_txt, cntxt)

    command_panel.itemClicked.connect(on_command_clicked)
    width = command_panel.size()
    hint = command_panel.sizeHint()
    frame_width = command_panel.frameWidth() + 2
    column_width = command_panel.sizeHintForColumn(0) + 2*frame_width
    command_panel.setMinimumWidth(column_width)
    command_panel.setMaximumWidth(column_width)
    command_panel.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Minimum)

    return command_panel


def on_command_clicked(item):
    cntxt = item.data(qt.EditRole)
    fct, args, kwargs = cntxt
    fct(*args, **kwargs)


def create_plot_view(app, figure, title, view_width, view_ht, commands=None,
                     add_panel_fcts=None, add_nav_toolbar=False,
                     drop_action=None, context_menu=None):
    """ create a window hosting a (mpl) figure """

    def create_light_or_dark_callback(figure):
        def l_or_d(is_dark):
            figure.sync_light_or_dark(is_dark)
        return l_or_d

    # construct the top level widget
    widget = QWidget()

    top_layout = QHBoxLayout(widget)
    # set the layout on the widget
    widget.setLayout(top_layout)
    widget.figure = figure

    if commands:
        command_panel = create_command_panel(figure, commands)
        top_layout.addWidget(command_panel)

    # construct the top level layout
    plot_layout = QVBoxLayout()
    top_layout.addLayout(plot_layout)

    pc = PlotCanvas(app, figure, drop_action=drop_action)
    if add_panel_fcts is not None:
        panel_layout = QHBoxLayout()
        plot_layout.addLayout(panel_layout)
        for p in add_panel_fcts:
            panel = p(app, pc)
            panel_layout.addWidget(panel)
        panel_layout.addStretch(50)

    if context_menu is not None:
        handle_context_menu = create_handle_context_menu(app, pc, context_menu)
        pc.setContextMenuPolicy(qt.CustomContextMenu)
        pc.customContextMenuRequested.connect(handle_context_menu)
    mi = ModelInfo(app.app_manager.model, update_figure_view, (figure,))
    sub_window = app.add_subwindow(widget, mi)
    sub_window.sync_light_or_dark = create_light_or_dark_callback(figure)
    lens_title = app.app_manager.model.name()
    sub_window.setWindowTitle(title + ': ' + lens_title)
    orig_x, orig_y = app.initial_window_offset()
    sub_window.setGeometry(orig_x, orig_y, view_width, view_ht)

    plot_layout.addWidget(pc)

    if add_nav_toolbar:
        plot_layout.addWidget(NavigationToolbar(pc, sub_window))

    sub_window.show()


def create_glass_map_view(app, glass_db):

    def create_light_or_dark_callback(fig):
        def l_or_d(is_dark):
            fig.sync_light_or_dark(is_dark)
        return l_or_d

    title = 'Glass Map Viewer'

    width = 1100
    height = 650

    # glass_db = gm.GlassMapDB(['Schott', 'Hoya', 'Ohara'])
    pdt = "Refractive Index"
    # hotwire GlassMapFigure to inherit from StyledFigure
    gm.GlassMapFigure.__bases__ = (StyledFigure,)
    figure = gm.GlassMapFigure(glass_db, plot_display_type=pdt,
                            # width=5, height=4,
                            )

    widget, pick_model = gmv.init_UI(app, figure)
    widget.figure = figure

    def refresh_gui(**kwargs):
        pick_model.fill_table(figure.pick_list)
    figure.refresh_gui = refresh_gui
    mi = ModelInfo(app.app_manager.model, update_figure_view, (figure,))
    sub_window = app.add_subwindow(widget, mi)
    sub_window.sync_light_or_dark = create_light_or_dark_callback(figure)
    sub_window.setWindowTitle(title)
    orig_x, orig_y = app.initial_window_offset()
    sub_window.setGeometry(orig_x, orig_y, width, height)

    sub_window.show()


def on_plot_scale_toggled(cntxt, scale_type):
    plotFigure, scale_wdgt = cntxt
    plotFigure.scale_type = scale_type
    if scale_type == Fit.User_Scale:
        scale_wdgt.setReadOnly(False)
        scale_wdgt.setText('{:.7g}'.format(plotFigure.user_scale_value))
    else:
        scale_wdgt.setReadOnly(True)
    plotFigure.plot()


def on_plot_scale_changed(cntxt):
    plotFigure, scale_wdgt = cntxt
    try:
        val = float(scale_wdgt.text())
        plotFigure.user_scale_value = val
        scale_wdgt.setText('{:.7g}'.format(val))
    except ValueError:
        return ''
    plotFigure.plot()


def create_plot_scale_panel(app, pc):
    groupBox = QGroupBox("Plot Scale", app)

    user_scale_wdgt = QLineEdit()
    user_scale_wdgt.setReadOnly(True)
    pf = pc.figure
    cntxt = pf, user_scale_wdgt
    user_scale_wdgt.editingFinished.connect(lambda:
                                            on_plot_scale_changed(cntxt))
    fit_btn = QRadioButton("Fit")
    fit_btn.setChecked(pf.scale_type == Fit.All)
    fit_btn.toggled.connect(lambda:
                                on_plot_scale_toggled(cntxt, Fit.All))
    user_scale_btn = QRadioButton("User Scale")
    user_scale_btn.setChecked(pf.scale_type == Fit.User_Scale)
    user_scale_btn.toggled.connect(lambda: on_plot_scale_toggled(
                                           cntxt, Fit.User_Scale))
    box = QHBoxLayout()
    box.addWidget(fit_btn)
    box.addWidget(user_scale_btn)
    box.addWidget(user_scale_wdgt)

    groupBox.setLayout(box)

    return groupBox


def create_multi_plot_scale_panel(app, pc):
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


def get_icon(fig, icon_filepath, icon_size=48):
    pm = QtGui.QPixmap(str(icon_filepath)).scaled(icon_size, icon_size)
    if hasattr(pm, 'setDevicePixelRatio'):
        try:  # older mpl < 3.5.0
            pm.setDevicePixelRatio(fig.canvas._dpi_ratio)
        except AttributeError:
            pm.setDevicePixelRatio(fig.canvas.device_pixel_ratio)
    return QtGui.QIcon(pm)


def create_2d_figure_toolbar(app, pc):
    """ zoom, fit and pan commands for figures
    Pan
    Zoom Box
    Fit
    Zoom In, Out
    1:1
    """
    icon_size = 24
    tb = QToolBar()
    tb.setIconSize(QSize(icon_size, icon_size))
    tb.setStyleSheet("QToolBar{spacing:0px;}")

    toolitems = (
        ('Pan', 'Pan axes with mouse', 'pan', 'register_pan'),
        ('Zoom', 'Zoom to rectangle', 'zoom', 'register_zoom_box'),
        ('Fit', 'Fit to view', 'fit', 'fit'),
        ('Zoom In', 'Zoom in', 'zoom_in', 'zoom_in'),
        ('Zoom Out', 'Zoom out', 'zoom_out', 'zoom_out'),
      )

    pth = Path(__file__).resolve().parent
    image_dir = Path(pth / 'images')

    fig = pc.figure
    actions = {}
    for text, tooltip_text, image_file, callback in toolitems:
        if text is None:
            tb.addSeparator()
        else:
            def create_action_and_on_finished(fig, callback):
                def make_cb_fct():
                    def on_finished():
                        # unpress button when action is finished
                        actions[callback].setChecked(False)
                    getattr(fig, callback)(on_finished)
                return make_cb_fct

            image_path = image_dir / (image_file + '.png')
            icon = get_icon(fig, image_path, icon_size=icon_size)

            if callback in ['register_zoom_box', 'register_pan']:
                # action requires mouse input handling, unpress when finished
                a = tb.addAction(icon, text,
                                 create_action_and_on_finished(fig, callback))
                actions[callback] = a
                a.setCheckable(True)
            else:
                # immediate action button
                a = tb.addAction(icon, text, getattr(fig, callback))
                actions[callback] = a

            if tooltip_text is not None:
                a.setToolTip(tooltip_text)

    def attr_check(fig, attr, state):
        checked = state == qt.Checked
        setattr(fig, attr, checked)
        fig.refresh()

    # add checkbox for unit aspect ration display
    aratio_checkBox = QCheckBox("1:1")
    aratio_checkBox.setChecked(fig.is_unit_aspect_ratio)
    aratio_checkBox.stateChanged.connect(
        lambda checked: attr_check(fig, 'is_unit_aspect_ratio', checked))
    tb.addWidget(aratio_checkBox)

    return tb


def create_draw_rays_groupbox(app, pc):
    tb = QToolBar()
    fig = pc.figure

    def attr_check(fig, attr, state):
        checked = state == qt.Checked
#        cur_value = getattr(fig, attr, None)
        setattr(fig, attr, checked)
        fig.refresh()

    parax_checkBox = QCheckBox("&paraxial rays")
    parax_checkBox.setChecked(fig.do_paraxial_layout)
    parax_checkBox.stateChanged.connect(
        lambda checked: attr_check(fig, 'do_paraxial_layout', checked))
    beam_checkBox = QCheckBox("&beams")
    beam_checkBox.setChecked(fig.do_draw_beams)
    beam_checkBox.stateChanged.connect(
        lambda checked: attr_check(fig, 'do_draw_beams', checked))
    edge_checkBox = QCheckBox("&edge rays")
    edge_checkBox.setChecked(fig.do_draw_edge_rays)
    edge_checkBox.stateChanged.connect(
        lambda checked: attr_check(fig, 'do_draw_edge_rays', checked))
    fan_checkBox = QCheckBox("&ray fans")
    fan_checkBox.setChecked(fig.do_draw_ray_fans)
    fan_checkBox.stateChanged.connect(
        lambda checked: attr_check(fig, 'do_draw_ray_fans', checked))
    clip_rays_checkBox = QCheckBox("&clip rays")
    clip_rays_checkBox.setChecked(fig.clip_rays)
    clip_rays_checkBox.stateChanged.connect(
        lambda checked: attr_check(fig, 'clip_rays', checked))
    
    tb.addWidget(parax_checkBox)
    tb.addWidget(beam_checkBox)
    tb.addWidget(edge_checkBox)
    tb.addWidget(fan_checkBox)
    tb.addWidget(clip_rays_checkBox)

    return tb


def create_diagram_controls_groupbox(app, pc):
    def on_barrel_constraint_toggled(cntxt, state):
        fig, barrel_wdgt = cntxt
        diagram = fig.diagram
        checked = state == qt.Checked
        if checked:
            diagram.do_barrel_constraint = True
            barrel_wdgt.setReadOnly(False)
            barrel_wdgt.setText('{:.7g}'.format(diagram.barrel_constraint_radius))
        else:
            diagram.do_barrel_constraint = False
            barrel_wdgt.setReadOnly(True)
    
        fig.refresh()
    
    def on_barrel_constraint_changed(cntxt):
        fig, barrel_wdgt = cntxt
        try:
            val = float(barrel_wdgt.text())
            fig.diagram.barrel_constraint_radius = val
            barrel_wdgt.setText('{:.7g}'.format(val))
        except ValueError:
            return ''
    
        fig.refresh()

    groupBox = QGroupBox("Controls", app)

    def attr_check(fig, attr, state):
        checked = state == qt.Checked
        # cur_value = getattr(fig, attr, None)
        setattr(fig, attr, checked)
        # print('attr_check: {}={}'.format(attr, checked))
        fig.refresh()

    barrel_value_wdgt = QLineEdit()
    barrel_value_wdgt.setReadOnly(True)
    fig = pc.figure
    cntxt = fig, barrel_value_wdgt
    barrel_value_wdgt.editingFinished.connect(lambda:
                                              on_barrel_constraint_changed(
                                                      cntxt))

    slide_checkBox = QCheckBox("&slide")
    slide_checkBox.setChecked(fig.enable_slide)
    slide_checkBox.stateChanged.connect(
        lambda checked: attr_check(fig, 'enable_slide', checked))
    barrel_checkBox = QCheckBox("&barrel constraint")
    barrel_checkBox.setChecked(fig.diagram.do_barrel_constraint)
    barrel_checkBox.stateChanged.connect(
        lambda checked: on_barrel_constraint_toggled(cntxt, checked))

    vbox = QVBoxLayout()
    vbox.addWidget(slide_checkBox)
    vbox.addWidget(barrel_checkBox)
    vbox.addWidget(barrel_value_wdgt)

    groupBox.setLayout(vbox)

    return groupBox


def create_diagram_edge_actions_groupbox(app, pc):
    def on_bend_or_gap_toggled(diagram, radio_btn_id):
        diagram.bend_or_gap = radio_btn_id

    fig = pc.figure
    diagram = fig.diagram

    groupBox = QGroupBox("Edge Actions", app)
    groupBox.setMaximumWidth(190)

    bending_btn = QRadioButton("Bend")
    bending_btn.setChecked(diagram.bend_or_gap == 'bend')
    bending_btn.toggled.connect(lambda:
                                on_bend_or_gap_toggled(diagram, 'bend'))
    gap_btn = QRadioButton("Gap")
    gap_btn.setChecked(diagram.bend_or_gap == 'gap')
    gap_btn.toggled.connect(lambda: on_bend_or_gap_toggled(diagram, 'gap'))

    vbox = QVBoxLayout()
    vbox.addWidget(bending_btn)
    vbox.addWidget(gap_btn)

    groupBox.setLayout(vbox)

    return groupBox


def create_diagram_layers_groupbox(app, pc):
    def on_active_diagram_toggled(fig, layer_key):
        fig.diagram.set_active_layer(layer_key)
        fig.refresh(build='rebuild')

    fig = pc.figure
    diagram = fig.diagram

    groupBox = QGroupBox("Layers", app)
    groupBox.setMaximumWidth(190)

    ifcs_btn = QRadioButton("interfaces")
    ifcs_btn.setChecked(diagram.active_layer == 'ifcs')
    ifcs_btn.toggled.connect(lambda:
                                on_active_diagram_toggled(fig, 'ifcs'))
    ele_btn = QRadioButton("elements")
    ele_btn.setChecked(diagram.active_layer == 'eles')
    ele_btn.toggled.connect(lambda: on_active_diagram_toggled(fig, 'eles'))
    asm_btn = QRadioButton("assembly")
    asm_btn.setChecked(diagram.active_layer == 'asm')
    asm_btn.toggled.connect(lambda: on_active_diagram_toggled(fig, 'asm'))
    sys_btn = QRadioButton("system")
    sys_btn.setChecked(diagram.active_layer == 'sys')
    sys_btn.toggled.connect(lambda: on_active_diagram_toggled(fig, 'sys'))

    vbox = QVBoxLayout()
    vbox.addWidget(ifcs_btn)
    vbox.addWidget(ele_btn)
    vbox.addWidget(asm_btn)
    vbox.addWidget(sys_btn)

    groupBox.setLayout(vbox)

    return groupBox


def create_handle_context_menu(app, pc, menu_actions):
    """ Create a (Qt) context menu for a (mpl) figure. 
    
    Args:

    - app: the Qt main app
    - pc: the Qt-specific backend for mpl, CanvasFigure
    - menu_actions: list of potential context menu actions

    each menu action is a tuple consisting of:

        - action_desc:str the menu item text
        - action_handle:str the key in the actions dict for the action
        - action: the action class to be instantiated if picked

    The fig.selected_shape is what the actions are built against. The action 
    should raise a TypeError if an incompatible shape is supplied.
    """
    fig = pc.figure

    def do_action(fig, shape, handle, action, info):
        shape.actions[handle] = action
        fig.selected_shape = (shape, handle), info

    def handle_context_menu(point):
        # show menu about the row
        menu = QMenu(app)
        selected_shape = fig.selected_shape
        if selected_shape is None:
            selected_artist = (fig.artist_infos[0] 
                               if len(fig.artist_infos) > 0 else None)
            if selected_artist is not None:
                selected_shape, info = (selected_artist.artist.shape, 
                                        selected_artist.info)
        else:
            selected_shape, info = fig.selected_shape

        if selected_shape is not None:
            for ma in menu_actions:
                action_desc, action_handle, action = ma
                try:
                    action_obj = action(selected_shape[0])
                except TypeError as e:
                    pass
                else:
                    menu.addAction(action_desc,
                                lambda: do_action(fig, selected_shape[0], 
                                                    action_handle, action_obj,
                                                    info))

            menu.popup(pc.mapToGlobal(point))

    return handle_context_menu
