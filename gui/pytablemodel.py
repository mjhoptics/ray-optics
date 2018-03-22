#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" Table model supporting data content via python eval() fct

Created on Wed Mar 14 21:59:43 2018

@author: Michael J. Hayford
"""
import logging

from PyQt5.QtCore import Qt, QAbstractTableModel
from PyQt5.QtCore import pyqtSignal


class PyTableModel(QAbstractTableModel):

    update = pyqtSignal(object, int)

    """ Model interface for table view of list structures """
    def __init__(self, rootObj, colEvalStr, rowHeaders, colHeaders,
                 colFormats, is_editable=False):
        """ Table model supporting data content via python eval() fct

        Initialization arguments:
            rootObj: object or list at the root of the eval() string
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
        """
        super(PyTableModel, self).__init__()
        self.root = rootObj
        self.colEvalStr = colEvalStr
        self.rowHeaders = rowHeaders
        self.colHeaders = colHeaders
        self.colFormats = colFormats
        self.is_editable = is_editable

    def rowCount(self, index):
        return len(self.rowHeaders)

    def columnCount(self, index):
        return len(self.colHeaders)

    def headerData(self, section, orientation, role):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return self.colHeaders[section]
            elif orientation == Qt.Vertical:
                return self.rowHeaders[section]
        else:
            return None

    def flags(self, index):
        base_flag = Qt.ItemIsEnabled | Qt.ItemIsSelectable
        if self.is_editable:
            return base_flag | Qt.ItemIsEditable
        else:
            return base_flag

    def data(self, index, role):
        if role == Qt.DisplayRole or role == Qt.EditRole:
            r = index.row()
            c = index.column()
            eval_str = ('self.root' + self.colEvalStr[c]).format(r)
            try:
                val = eval(eval_str)
                valStr = self.colFormats[c].format(val)
                return valStr
            except IndexError:
                return ''
        else:
            return None

    def setData(self, index, value, role):
        if role == Qt.EditRole:
            r = index.row()
            c = index.column()
            exec_str = ('self.root' + self.colEvalStr[c]).format(r)
            exec_str = exec_str + '=' + value
            try:
                exec(exec_str)
                self.update.emit(self.root, r)
                return True
            except IndexError:
                return False
            except SyntaxError:
                logging.info('Syntax error: "%s"', value)
                return False
        else:
            return False
