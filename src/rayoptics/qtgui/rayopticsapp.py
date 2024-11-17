#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Ray Optics GUI Application

Relies on PyQt5

.. Created on Mon Feb 12 09:24:01 2018

.. codeauthor: Michael J. Hayford
"""

import sys
import logging
from pathlib import Path

from PyQt5.QtCore import Qt
from PyQt5.QtCore import QEvent
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import (QApplication, QAction, QMainWindow, QMdiArea,
                             QFileDialog, QWidget, QMenu,
                             QVBoxLayout)
from PyQt5.QtCore import pyqtSlot
import qdarkstyle

from traitlets.config.configurable import MultipleInstanceError

import rayoptics
from rayoptics.raytr.trace import RaySeg
import rayoptics.gui.appcmds as cmds
from rayoptics.gui.appmanager import ModelInfo, AppManager
import rayoptics.qtgui.dockpanels as dock
from rayoptics.qtgui.ipyconsole import create_ipython_console
from rayoptics.qtgui.pytableview import TableView

from rayoptics.parax import firstorder
from rayoptics.raytr import trace
from rayoptics.raytr import traceerror as terr

logger = logging.getLogger(__name__)


class MainWindow(QMainWindow):
    count = 0

    def __init__(self, parent=None, qtapp=None):
        super().__init__(parent)
        self.qtapp = qtapp
        self.mdi = QMdiArea()
        self.setCentralWidget(self.mdi)

        self.app_manager = AppManager(None, gui_parent=self)
        self.mdi.subWindowActivated.connect(self.app_manager.
                                            on_view_activated)

        self.is_dark = self.light_or_dark(False)

        self.left = 100
        self.top = 50
        self.width = 2100
        self.height = 1200
        self.setGeometry(self.left, self.top, self.width, self.height)

        bar = self.menuBar()

        file_menu = bar.addMenu("File")
        file_menu.addAction("New")
        file_menu.addAction("New Diagram")
        # file_menu.addAction("New Spec Sheet")
        file_menu.addAction("Open...")
        file_menu.addSeparator()
        file_menu.addAction("Save")
        file_menu.addAction("Save As...")
        file_menu.addAction("Close")
        file_menu.triggered[QAction].connect(self.do_file_action)

        view_menu = bar.addMenu("Data View")
        view_menu.addAction("Spec Sheet")
        view_menu.addAction("Optical Layout")
        view_menu.addAction("Lens Table")
        view_menu.addAction("Element Table")
        view_menu.addAction("Glass Map")
        # view_menu.addAction("Lens View")
        view_menu.triggered[QAction].connect(self.do_view_action)

        parax_menu = bar.addMenu("Paraxial Model")
        parax_menu.addAction("Paraxial Model")
        parax_menu.addAction("y-ybar View")
        parax_menu.addAction("nu-nubar View")
        parax_menu.addAction("yui Ray Table")
        parax_menu.addAction("3rd Order Aberrations")
        parax_menu.triggered[QAction].connect(self.do_view_action)

        analysis_menu = bar.addMenu("Analysis")
        analysis_menu.addAction("Ray Table")
        analysis_menu.addAction("Ray Fans")
        analysis_menu.addAction("OPD Fans")
        analysis_menu.addAction("Spot Diagram")
        analysis_menu.addAction("Wavefront Map")
        analysis_menu.addAction("Astigmatism Curves")
        analysis_menu.triggered[QAction].connect(self.do_view_action)

        tools_menu = bar.addMenu("Tools")
        tools_menu.addAction("Refocus")
        tools_menu.addAction("Set Vignetting")
        tools_menu.addAction("Set Apertures")
        tools_menu.addAction("Set Pupil")
        tools_menu.triggered[QAction].connect(self.do_view_action)

        wnd_menu = bar.addMenu("Window")
        wnd_menu.addAction("Cascade")
        wnd_menu.addAction("Tiled")
        wnd_menu.addSeparator()
        wnd_menu.addAction("Light UI")
        wnd_menu.addAction("Dark UI")
        wnd_menu.addSeparator()

        dock.create_dock_windows(self)
        for pi in dock.panels.values():
            wnd_menu.addAction(pi.menu_action)

        wnd_menu.triggered[QAction].connect(self.do_window_action)

        self.setWindowTitle("Ray Optics")
        self.show()

        path = Path(rayoptics.__file__).parent
        self.cur_dir = path / "models"

        if True:
            # create new model
            # self.new_model()
            # self.new_model_via_diagram()
            # self.new_model_via_specsheet()
            self.new_empty_model()

        else:
            # restore a default model

            # self.cur_dir = path / "codev/tests"
            # self.open_file(path / "codev/tests/asp46.seq")
            # self.open_file(path / "codev/tests/achroMangin.seq")
            # self.open_file(path / "codev/tests/dar_test.seq")
            # self.open_file(path / "codev/tests/dec_test.seq")
            # self.open_file(path / "codev/tests/paraboloid.seq")
            # self.open_file(path / "codev/tests/paraboloid_f8.seq")
            # self.open_file(path / "codev/tests/schmidt.seq")
            # self.open_file(path / "codev/tests/questar35.seq")
            # self.open_file(path / "codev/tests/rc_f16.seq")
            # self.open_file(path / "codev/tests/ag_dblgauss.seq")
            # self.open_file(path / "codev/tests/threemir.seq")
            # self.open_file(path / "codev/tests/folded_lenses.seq")
            # self.open_file(path / "codev/tests/lens_reflection_test.seq")
            # self.open_file(path / "codev/tests/dec_tilt_test.seq")
            # self.open_file(path / "codev/tests/tilt_test.seq")
            # self.open_file(path / "codev/tests/landscape_lens.seq")
            # self.open_file(path / "codev/tests/mangin.seq")
            # self.open_file(path / "codev/tests/CODV_32327.seq")
            # self.open_file(path / "codev/tests/CODV_65988.seq")
            # self.open_file(path / "codev/tests/questar35.seq")
            # self.open_file(path / "codev/tests/dar_test.seq")

            # self.cur_dir = path / "optical/tests"
            # self.open_file(path / "optical/tests/achroMangin.roa")
            # self.open_file(path / "optical/tests/cell_phone_camera.roa")
            # self.open_file(path / "optical/tests/singlet_f3.roa")
            # self.open_file(path / "optical/tests/Nikon Nikkor Z 14-30mm f-4 S.roa")
            # self.open_file(path / "optical/tests/US007277232_Example04P.roa")
            # self.open_file(path / "optical/tests/US05831776-1.zmx")

            self.cur_dir = path / "models"
            # self.open_file(path / "models/Cassegrain.roa")
            # self.open_file(path / "models/collimator.roa")
            # self.open_file(path / "models/Dall-Kirkham.roa")
            # self.open_file(path / "models/HybridAchromat.roa")
            # self.open_file(path / "models/Newtonian with diagonal.roa")
            # self.open_file(path / "models/petzval.roa")
            # self.open_file(path / "models/Ritchey_Chretien.roa")
            self.open_file(path / "models/Sasian Triplet.roa")
            # self.open_file(path / "models/singlet_f5.roa")
            # self.open_file(path / "models/thinlens.roa")
            # self.open_file(path / "models/telephoto.roa")
            # self.open_file(path / "models/thin_triplet.roa")
            # self.open_file(path / "models/galilean.roa")
            # self.open_file(path / "models/TwoMirror.roa")
            # self.open_file(path / "models/TwoSphericalMirror.roa")

            # self.cur_dir = path / "zemax/tests"
            # self.open_file(path / "zemax/tests/US08427765-1.ZMX")
            # self.open_file(path / "zemax/tests/US00583336-2-scaled.zmx")
            # self.open_file(path / "zemax/tests/HoO-V2C18Ex03.zmx")
            # self.open_file(path / "zemax/tests/HoO-V2C18Ex27.zmx")
            # self.open_file(path / "zemax/tests/HoO-V2C18Ex46.zmx")
            # self.open_file(path / "zemax/tests/HoO-V2C18Ex66.zmx")
            # self.open_file(path / "zemax/tests/US05831776-1.zmx")  # litho lens
            # self.open_file(path / "zemax/tests/zmax_37992.zmx")
            # self.open_file(path / "zemax/tests/354710-C-Zemax(ZMX).zmx")

            # self.cur_dir = path / "zemax/models/telescopes"
            # self.open_file(path / "zemax/models/telescopes/Figure4.zmx")
            # self.open_file(path / "zemax/models/telescopes/HoO-V2C18Ex11.zmx")

            # self.cur_dir = path / "zemax/models/PhotoPrime"
            # self.open_file(path / "zemax/models/PhotoPrime/US05321554-4.ZMX")
            # self.open_file(path / "zemax/models/PhotoPrime/US06476982-1.ZMX")
            # self.open_file(path / "zemax/models/PhotoPrime/US07190532-1.ZMX")
            # self.open_file(path / "zemax/models/PhotoPrime/US04331391-1.zmx")
            # self.open_file(path / "zemax/models/PhotoPrime/US05331467-1.zmx")
            # self.open_file(path / "zemax/models/PhotoPrime/US05331467-1_asm.roa")

            # root_pth = Path("/Users/Mike/Developer/PyProjects") 
            # ro_test_files = root_pth / 'ro_test_files'
            # self.cur_dir = ro_test_files / "optical"
            # self.open_file(ro_test_files / "optical/Nikon Nikkor Z 14-30mm f-4 S real glasses.roa")

            # lens_designs_path = root_pth / 'lens-designs-dotcom'
            # self.cur_dir = lens_path = lens_designs_path / "PhotoPrime"
            # self.open_file(lens_path / "US007161746_Example09P.zmx")
            # self.open_file(lens_path / "US20190086648_Example02P.zmx")
            # self.open_file(lens_path / "US20130235467_Example01P.zmx")

    def add_subwindow(self, widget, model_info):
        sub_wind = self.mdi.addSubWindow(widget)
        sub_wind.installEventFilter(self)
        self.app_manager.add_view(sub_wind, widget, model_info)
        MainWindow.count += 1
        return sub_wind

    def delete_subwindow(self, sub_wind):
        self.app_manager.delete_view(sub_wind)
        self.mdi.removeSubWindow(sub_wind)
        MainWindow.count -= 1

    def add_ipython_subwindow(self, opt_model):
        try:
            title_bar = 'iPython console: '
            if opt_model is not None:
                title_bar = title_bar + opt_model.name()
            view_width, view_ht = 900, 600
            orig_x = dock.panels['wavelengths'].dock.pos().x()-view_width
            orig_y = 0
            create_ipython_console(self, opt_model, title_bar, 
                                   orig_x, orig_y, view_width, view_ht, ) 

        except MultipleInstanceError:
            logger.debug("Unable to open iPython console. "
                         "MultipleInstanceError")
        except Exception as inst:
            print(type(inst))    # the exception instance
            print(inst.args)     # arguments stored in .args
            print(inst)          # __str__ allows args to be printed directly,
            pass                 # but may be overridden in exception subclasses

    def initial_window_offset(self):
        offset_x = 50
        offset_y = 25
        orig_x = (MainWindow.count - 1)*offset_x
        orig_y = (MainWindow.count - 1)*offset_y
        return orig_x, orig_y

    def do_file_action(self, q):
        self.file_action(q.text())

    def file_action(self, action):
        if action == "New":
            self.new_empty_model()

        if action == "New Diagram":
            self.new_model_via_diagram()

        # if action == "New Spec Sheet":
        #     self.new_model_via_specsheet()

        if action == "Open...":
            options = QFileDialog.Options()
            # options |= QFileDialog.DontUseNativeDialog
            fileName, _ = QFileDialog.getOpenFileName(
                          self,
                          "Open File ...",
                          str(self.cur_dir),
                          "All files (*.seq *.zmx *.roa);;"
                          "CODE V files (*.seq);;"
                          "Ray-Optics files (*.roa);;"
                          "Zemax files (*.zmx)",
                          options=options)
            if fileName:
                logger.debug("open file: %s", fileName)
                filename = Path(fileName)
                self.cur_dir = filename.parent
                self.open_file(filename)

        if action == "Save":
            cur_model = self.app_manager.model
            fileName = self.app_manager.model_filenames[cur_model]
            if fileName:
                logger.debug("save file: %s", fileName)
                self.save_file(fileName)
            else:  # query for filename, i.e. do a SaveAs
                action = "Save As..."

        if action == "Save As...":
            options = QFileDialog.Options()
            # options |= QFileDialog.DontUseNativeDialog
            fileName, _ = QFileDialog.getSaveFileName(
                          self,
                          "Save File ...",
                          str(self.cur_dir),
                          "Ray-Optics Files (*.roa);;All Files (*)",
                          options=options)
            if fileName:
                logger.debug("save file: %s", fileName)
                filename = Path(fileName)
                self.cur_dir = filename.parent
                self.save_file(fileName)

        if action == "Close":
            self.close_model()

    def new_model(self, **kwargs):
        opt_model = cmds.create_new_optical_system(**kwargs)
        self.app_manager.set_model(opt_model)
        self.refresh_gui()

        self.create_lens_table()
        cmds.create_live_layout_view(opt_model, gui_parent=self)
        cmds.create_paraxial_design_view_v2(opt_model, 'ht',
                                            gui_parent=self)
        self.refresh_gui()

        self.new_ipython_console(opt_model)

    def new_empty_model(self, **kwargs):
        opt_model = cmds.create_empty_model(**kwargs)
        self.app_manager.set_model(opt_model)
        self.new_ipython_console(opt_model)

    def new_model_via_specsheet(self):
        cmds.create_new_ideal_imager_dialog(gui_parent=self,
                                            conjugate_type='infinite')
        self.new_ipython_console(None)

    def new_model_via_diagram(self, **kwargs):
        """ Define a new model using a |ybar| diagram. """
        opt_model = cmds.create_empty_model(**kwargs)
        self.app_manager.set_model(opt_model)

        cmds.create_paraxial_design_view_v2(opt_model, 'ht', gui_parent=self)
        self.new_ipython_console(opt_model)

    def new_ipython_console(self, opt_model):
        self.add_ipython_subwindow(opt_model)
        self.refresh_app_ui()

    def open_file(self, file_name, **kwargs):
        cur_model = cmds.open_model(file_name, **kwargs)
        self.app_manager.set_model(cur_model, filename=file_name)
        self.create_lens_table()
        cmds.create_live_layout_view(cur_model, gui_parent=self)
        self.add_ipython_subwindow(cur_model)
        self.refresh_app_ui()

    def save_file(self, file_name):
        cur_model = self.app_manager.model
        cur_model.save_model(file_name, version="0.9.0a1")
        self.app_manager.model_filenames[cur_model] = file_name

    def close_model(self):
        """ NOTE: this does not check to save a modified model """
        self.app_manager.close_model(self.delete_subwindow)

    def do_view_action(self, q):
        self.view_action(q.text())

    def view_action(self, action):
        opt_model = self.app_manager.model

        if action == "Spec Sheet":
            cmds.create_new_ideal_imager_dialog(opt_model=opt_model,
                                                gui_parent=self)

        if action == "Optical Layout":
            cmds.create_live_layout_view(opt_model, gui_parent=self)

        if action == "Lens Table":
            self.create_lens_table()

        if action == "Element Table":
            model = cmds.create_element_table_model(opt_model)
            self.create_table_view(model, "Element Table")

        if action == "Glass Map":
            cmds.create_glass_map_view(opt_model, gui_parent=self)

        if action == "Ray Fans":
            cmds.create_ray_fan_view(opt_model, "Ray", gui_parent=self)

        if action == "OPD Fans":
            cmds.create_ray_fan_view(opt_model, "OPD", gui_parent=self)

        if action == "Spot Diagram":
            cmds.create_ray_grid_view(opt_model, gui_parent=self)

        if action == "Wavefront Map":
            cmds.create_wavefront_view(opt_model, gui_parent=self)

        if action == "Astigmatism Curves":
            cmds.create_field_curves(opt_model, gui_parent=self)

        if action == "3rd Order Aberrations":
            cmds.create_3rd_order_bar_chart(opt_model, gui_parent=self)

        if action == "y-ybar View":
            cmds.create_paraxial_design_view_v2(opt_model, 'ht',
                                                gui_parent=self)

        if action == "nu-nubar View":
            cmds.create_paraxial_design_view_v2(opt_model, 'slp',
                                                gui_parent=self)

        if action == "yui Ray Table":
            model = cmds.create_parax_table_model(opt_model)
            self.create_table_view(model, "Paraxial Ray Table")

        if action == "Paraxial Model":
            model = cmds.create_parax_model_table(opt_model)
            self.create_table_view(model, "Paraxial Model")

        if action == "Ray Table":
            self.create_ray_table(opt_model)

        if action == "Set Vignetting":
            cmds.set_vignetting(opt_model, gui_parent=self)

        if action == "Set Apertures":
            cmds.set_apertures(opt_model, gui_parent=self)

        if action == "Set Pupil":
            cmds.set_pupil(opt_model, gui_parent=self)

        if action == "Refocus":
            cmds.refocus(opt_model, gui_parent=self)

    def do_window_action(self, q):
        self.window_action(q.text())

    def window_action(self, action):
        if action == "Cascade":
            self.mdi.cascadeSubWindows()

        if action == "Tiled":
            self.mdi.tileSubWindows()

        if action == "Light UI":
            self.is_dark = self.light_or_dark(False)
            self.app_manager.sync_light_or_dark(self.is_dark)

        if action == "Dark UI":
            self.is_dark = self.light_or_dark(True)
            self.app_manager.sync_light_or_dark(self.is_dark)

    def light_or_dark(self, is_dark):
        """ set the UI to a light or dark scheme.

        Qt doesn't seem to support controlling the MdiArea's background from a
        style sheet. Set the widget directly and save the original color
        to reset defaults.
        """
        if not hasattr(self, 'mdi_background'):
            self.mdi_background = self.mdi.background()

        if is_dark:
            self.qtapp.setStyleSheet(qdarkstyle.load_stylesheet(qt_api='pyqt5'))
            colors = qdarkstyle.dark.palette.DarkPalette.to_dict()
            self.mdi.setBackground(QColor(colors['COLOR_BACKGROUND_2']))
        else:
            self.qtapp.setStyleSheet('')
            self.mdi.setBackground(self.mdi_background)
        return is_dark

    def create_lens_table(self):
        seq_model = self.app_manager.model.seq_model

        def set_stop_surface(stop_surface):
            seq_model.stop_surface = stop_surface
            self.refresh_gui()

        def handle_context_menu(point):
            try:
                vheader = view.verticalHeader()
                row = vheader.logicalIndexAt(point.y())
            except NameError:
                pass
            else:
                # show menu about the row
                menu = QMenu(self)
                if row != seq_model.stop_surface:
                    menu.addAction('Set Stop Surface',
                                   lambda: set_stop_surface(row))
                if seq_model.stop_surface is not None:
                    menu.addAction('Float Stop Surface',
                                   lambda: set_stop_surface(None))
                menu.popup(vheader.mapToGlobal(point))

        model = cmds.create_lens_table_model(seq_model)
        view = self.create_table_view(model, "Surface Data Table")
        vheader = view.verticalHeader()
        vheader.setContextMenuPolicy(Qt.CustomContextMenu)
        vheader.customContextMenuRequested.connect(handle_context_menu)

    def create_ray_table(self, opt_model):
        osp = opt_model.optical_spec
        pupil = [0., 1.]
        fi = 0
        wl = osp.spectral_region.reference_wvl
        fld, wvl, foc = osp.lookup_fld_wvl_focus(fi, wl)
        ray_pkg, ray_error = trace.trace_ray(opt_model, pupil, fld, wvl)
 
#        ray, ray_op, wvl, opd = trace.trace_with_opd(opt_model, pupil,
#                                                     fld, wvl, foc)

#        cr = trace.RayPkg(ray, ray_op, wvl)
#        s, t = trace.trace_coddington_fan(opt_model, cr, foc)
        model = cmds.create_ray_table_model(opt_model, ray_pkg.ray)
        self.create_table_view(model, "Ray Table")

    def create_table_view(self, table_model, table_title, 
                          close_callback=None,
                          update_callback=None):
        # construct the top level widget
        widget = QWidget()
        # construct the top level layout
        layout = QVBoxLayout(widget)

        table_view = TableView(table_model)
        table_view.setAlternatingRowColors(True)

        # Add table to box layout
        layout.addWidget(table_view)

        # set the layout on the widget
        widget.setLayout(layout)

        sub = self.add_subwindow(widget, ModelInfo(self.app_manager.model,
                                                   cmds.update_table_view,
                                                   (table_view,), {}))
        lens_title = self.app_manager.model.name()
        sub.setWindowTitle(table_title + ': ' + lens_title)

        table_view.setMinimumWidth(table_view.horizontalHeader().length() +
                                   table_view.horizontalHeader().height())
#                                  The following line should work but returns 0
#                                  table_view.verticalHeader().width())

        view_width = table_view.width()
        view_ht = table_view.height()
        orig_x, orig_y = self.initial_window_offset()
        sub.setGeometry(orig_x, orig_y, view_width, view_ht)

        # table data updated successfully
        if update_callback is None:
            update_callback = self.on_data_changed
        table_model.update.connect(update_callback)

        sub.show()

        return table_view

    def eventFilter(self, obj, event):
        """Used by subwindows in response to installEventFilter."""
        if (event.type() == QEvent.Close):
            logger.debug(f"close event received: {obj}")
            self.delete_subwindow(obj)
        return False

    def refresh_gui(self, **kwargs):
        self.app_manager.refresh_gui(**kwargs)

    def refresh_app_ui(self):
        dock.update_dock_windows(self)

    def handle_ideal_imager_command(self, iid, command, specsheet):
        ''' link Ideal Imager Dialog buttons to model actions
        iid: ideal imager dialog
        command: text field with the action - same as button label
        specsheet: the input specsheet used to drive the actions
        '''
        if command == 'Apply':
            opt_model = self.app_manager.model
            opt_model.set_from_specsheet(specsheet)
            self.refresh_gui()
        elif command == 'Close':
            for view, info in self.app_manager.view_dict.items():
                if iid == info[0]:
                    self.delete_subwindow(view)
                    view.close()
                    break
        elif command == 'Update':
            opt_model = self.app_manager.model
            specsheet = opt_model.specsheet
            firstorder.specsheet_from_parax_data(opt_model, specsheet)
            iid.specsheet_dict[specsheet.conjugate_type] = specsheet
            iid.update_values()
        elif command == 'New':
            opt_model = cmds.create_new_optical_model_from_specsheet(specsheet)
            self.app_manager.set_model(opt_model)
            for view, info in self.app_manager.view_dict.items():
                if iid == info[0]:
                    w = iid
                    mi = info[1]
                    args = (iid, opt_model)
                    new_mi = ModelInfo(model=opt_model, fct=mi.fct,
                                       args=args, kwargs=mi.kwargs)
                    self.app_manager.view_dict[view] = w, new_mi
            self.refresh_gui()
            self.create_lens_table()
            cmds.create_live_layout_view(opt_model, gui_parent=self)
            cmds.create_paraxial_design_view_v2(opt_model, 'ht',
                                                gui_parent=self)
            self.refresh_gui()

    @pyqtSlot(object, int)
    def on_data_changed(self, rootObj, index):
        self.refresh_gui(src_model=rootObj)


def main():
    logging_level = logging.INFO
    try:
        logging.basicConfig(filename='rayoptics.log',
                            filemode='w',
                            level=logging_level)
    except:
        logging.basicConfig(filename=Path.home().joinpath('rayoptics.log'),
                            filemode='w',
                            level=logging_level)
        
    qtapp = QApplication(sys.argv)
    qtwnd = MainWindow(qtapp=qtapp)

    qtwnd.show()
    return qtapp.exec_()


if __name__ == '__main__':
    sys.exit(main())
