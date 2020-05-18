#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Support creation of an iPython console, with rayoptics environment

.. Created on Wed Nov 21 21:48:02 2018

.. codeauthor: Michael J. Hayford
"""

from qtconsole.rich_jupyter_widget import RichJupyterWidget
from qtconsole.inprocess import QtInProcessKernelManager
from IPython.lib import guisupport

from rayoptics.gui.appmanager import ModelInfo


def create_ipython_console(app, title, view_width, view_ht):
    """ create a iPython console with a rayoptics environment """
    opt_model = app.app_manager.model
    if opt_model:
        ro_env = {
                'app': app,
                'opm': opt_model,
                'sm': opt_model.seq_model,
                'osp': opt_model.optical_spec,
                'pm': opt_model.parax_model
                }
    else:
        ro_env = {
            'app': app,
            'opm': opt_model
            }

    ro_setup = 'from rayoptics.environment import *'

    # construct the top level widget
    ipy_console = ConsoleWidget()

    # load the environment
    ipy_console.execute_command(ro_setup)
    ipy_console.push_vars(ro_env)

    mi = ModelInfo(opt_model)
    sub_window = app.add_subwindow(ipy_console, mi)
    sub_window.setWindowTitle(title)
    orig_x, orig_y = app.initial_window_offset()
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
