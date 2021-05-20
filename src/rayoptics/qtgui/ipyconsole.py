#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Support creation of an iPython console, with rayoptics environment

.. Created on Wed Nov 21 21:48:02 2018

.. codeauthor: Michael J. Hayford
"""

from PyQt5.QtGui import QColor
from qtconsole.rich_jupyter_widget import RichJupyterWidget
from qtconsole.inprocess import QtInProcessKernelManager
from qtconsole import styles
from IPython.lib import guisupport

import qdarkstyle

from rayoptics.gui.appmanager import ModelInfo
from rayoptics.util import colors

default_template = '''\
    QPlainTextEdit, QTextEdit {
        background-color: %(bgcolor)s;
        background-clip: padding;
        color: %(fgcolor)s;
        selection-background-color: %(select)s;
    }
    .inverted {
        background-color: %(fgcolor)s;
        color: %(bgcolor)s;
    }
    .error { color: red; }
    .in-prompt { color: %(i_color)s; }
    .in-prompt-number { font-weight: bold; }
    .out-prompt { color: %(o_color)s; }
    .out-prompt-number { font-weight: bold; }
'''


def create_ipython_console(gui_parent, opt_model, title, view_width, view_ht):
    """ create a iPython console with a rayoptics environment """

    def create_light_or_dark_callback(ipy_console):
        # if not hasattr(ipy_console, 'background'):
        #     ipy_console.background = ipy_console.background()

        def l_or_d(is_dark):
            accent = colors.accent_colors(is_dark)
            prompt_style = {
                'i_color': accent['cyan'],
                'o_color': accent['orange'],
            }
            if is_dark:
                ipy_console.setStyleSheet(
                    qdarkstyle.load_stylesheet(qt_api='pyqt5'))
                # color_defs = {**styles.get_colors('solarized-dark'),
                #               **prompt_style }
            else:
                ipy_console.setStyleSheet('')
                # color_defs = {**styles.get_colors('solarized-light'),
                #               **prompt_style }
            # ipy_console.style_sheet = default_template%color_defs
            # ipy_console._style_sheet_changed()
        return l_or_d

    if opt_model:
        ro_env = {
                'gui_parent': gui_parent,
                'app': gui_parent.app_manager,
                'opm': opt_model,
                'sm': opt_model.seq_model,
                'osp': opt_model.optical_spec,
                'pm': opt_model.parax_model,
                'em': opt_model.ele_model,
                'pt': opt_model.part_tree,
                }
    else:
        ro_env = {
            'gui_parent': gui_parent,
            'app': gui_parent.app_manager,
            }

    ro_setup = 'from rayoptics.environment import *'

    # construct the top level widget
    ipy_console = ConsoleWidget()

    # load the environment
    ipy_console.execute_command(ro_setup)
    ipy_console.push_vars(ro_env)

    mi = ModelInfo(opt_model)
    sub_window = gui_parent.add_subwindow(ipy_console, mi)
    sub_window.setWindowTitle(title)
    sub_window.sync_light_or_dark = create_light_or_dark_callback(ipy_console)
    orig_x, orig_y = gui_parent.initial_window_offset()
    sub_window.setGeometry(orig_x, orig_y, view_width, view_ht)

    sub_window.show()


class ConsoleWidget(RichJupyterWidget):

    def __init__(self, customBanner=None, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if customBanner is not None:
            self.banner = customBanner

        self.font_size = 6
        self.kernel_manager = kernel_manager = QtInProcessKernelManager()
        kernel_manager.start_kernel(show_banner=False)
        kernel_manager.kernel.gui = 'qt'
        self.kernel_client = kernel_client = self.kernel_manager.client()
        kernel_client.start_channels()

        def stop():
            kernel_client.stop_channels()
            kernel_manager.shutdown_kernel()
            guisupport.get_app_qt().exit()

        self.exit_requested.connect(stop)

    def push_vars(self, variableDict):
        """
        Given a dictionary containing name / value pairs, push those variables
        to the Jupyter console widget
        """
        self.kernel_manager.kernel.shell.push(variableDict)

    def clear(self):
        """
        Clears the terminal
        """
        self._control.clear()

        # self.kernel_manager

    def print_text(self, text):
        """
        Prints some plain text to the console
        """
        self._append_plain_text(text)

    def execute_command(self, command):
        """
        Execute a command in the frame of the console widget
        """
        self._execute(command, False)
