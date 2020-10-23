#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2020 Michael J. Hayford
"""

.. Created on Thu Oct  8 15:40:47 2020

.. codeauthor: Michael J. Hayford
"""


from PyQt5.QtCore import Qt as qt

from PyQt5.QtWidgets import (QTableView)


class TableView(QTableView):

    def __init__(self, table_model, accept_drops=True):
        super().__init__()
        self.setModel(table_model)
        self.setAcceptDrops(True)

        # Next 2 lines are needed so that key press events are correctly
        #  passed with mouse events
        # https://github.com/matplotlib/matplotlib/issues/707/
        self.setFocusPolicy(qt.ClickFocus)
        self.setFocus()

        self.drop_action = GlassDropAction()

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


class GlassDropAction():
    def __init__(self):
        self.drop_action = None

    def dragEnterEvent(self, view, event):
        self.drop_action = None

    def dragMoveEvent(self, view, event):
        indx = view.indexAt(event.pos())
        self.drop_action = view.model().drop_actions[indx.column()]

    def dragLeaveEvent(self, view, event):
        self.drop_action = None

    def dropEvent(self, view, event):
        if self.drop_action is not None:
            index = view.indexAt(event.pos())
            idx = index.row()
            self.drop_action(event, idx)
            model = view.model()
            model.update.emit(model.get_root_object(), idx)
            return True
        return False
