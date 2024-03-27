#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Table model supporting data content via python eval() fct

.. Created on Wed Mar 14 21:59:43 2018

.. codeauthor: Michael J. Hayford
"""
import logging

from PyQt5.QtCore import Qt, QAbstractTableModel
from PyQt5.QtCore import pyqtSignal

from rayoptics.util.misc_math import isanumber

logger = logging.getLogger(__name__)


class PyTableModel(QAbstractTableModel):
    """Table model supporting data content via python eval() fct.

    Model interface for table view of list structures.

    Attributes:
        root: object or list at the root of the eval() string
        rootEvalStr: string that is concatentated to the root name and
                    passed to the eval() function. This will accomodate
                    dynamic name changes.
        colEvalStr: string that is concatentated to the root name and
                    passed to the eval() function. There should be a
                    replacement field, i.e. {} where the row value will
                    be substituted using the str.format() function.
        rowHeaders: list of strings, length defines number of rows in the
                    table
        colHeaders: list of strings, length defines number of columns in
                    the table
        colFormats: format strings to be used to format data in each column
        is_editable: if true, items are editable
        get_num_rows: if not None, a function that returns the number of
                      rows in the table
        get_row_headers: if not None, a function that returns the row
                         headers for the table
    """

    update = pyqtSignal(object, int)

    def __init__(self, root, rootEvalStr, colEvalStr, rowHeaders,
                 colHeaders, colFormats, is_editable=False, get_num_rows=None,
                 get_row_headers=None, drop_actions=None):
        """
    """
        super().__init__()
        self.root = root
        self.rootEvalStr = rootEvalStr
        self.colEvalStr = colEvalStr
        self.rowHeaders = rowHeaders
        self.colHeaders = colHeaders
        self.colFormats = colFormats
        self.is_editable = is_editable
        self.get_num_rows = get_num_rows
        self.get_row_headers = get_row_headers
        if drop_actions:
            self.drop_actions = drop_actions
        else:
            self.drop_actions = [None]*len(self.colHeaders)

    def rowCount(self, index):
        if self.get_num_rows is not None:
            return self.get_num_rows()
        else:
            return len(self.rowHeaders)

    def columnCount(self, index):
        return len(self.colHeaders)

    def headerData(self, section, orientation, role):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return self.colHeaders[section]
            elif orientation == Qt.Vertical:
                if self.get_row_headers is not None:
                    self.rowHeaders = self.get_row_headers()
                if len(self.rowHeaders) == 0:
                    return None
                return self.rowHeaders[section]
        else:
            return None

    def flags(self, index):
        base_flag = Qt.ItemIsEnabled | Qt.ItemIsSelectable
        if self.is_editable:
            return base_flag | Qt.ItemIsEditable
        else:
            return base_flag

    def get_root_object(self):
        if len(self.rootEvalStr) == 0:
            return self.root
        else:
            root_eval_str = ('self.root' + self.rootEvalStr)
            try:
                root = eval(root_eval_str)
                return root
            except IndexError:
                return self.root

    def data(self, index, role):
        root = self.get_root_object()
        if role == Qt.DisplayRole or role == Qt.EditRole:
            r = index.row()
            c = index.column()
            eval_str = ('root' + self.colEvalStr[c]).format(r)
            try:
                val = eval(eval_str)
                valStr = self.colFormats[c].format(val)
                return valStr
            except IndexError:
                return ''
            except TypeError:
                print('Data type error: ', eval_str, val)
                return ''
        else:
            return None

    def setData(self, index, value, role):
        root = self.get_root_object()
        if role == Qt.EditRole:
            r = index.row()
            c = index.column()
            exec_str = ('root' + self.colEvalStr[c]).format(r)
            if not isanumber(value):
                value = "'" + value + "'"
            exec_str = exec_str + '=' + value
            try:
                exec(exec_str)
                self.update.emit(root, r)
                return True
            except IndexError:
                return False
            except SyntaxError:
                logger.info('Syntax error: "%s"', value)
                return False
        else:
            return False
