#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2023 Michael J. Hayford
""" Module for parsing a sequential model into elements

The grammar is as follows::

    seq_model   = object space optics image
    optics      = (part space)+

    part        = mangin / cemented / lens / mirror / surface / dummy / thin_lens
    space        = ~r"a|t"

    surface      = "i"
    lens         = "iti"
    mirror       = "r"
    air          = "a"
    cemented     = "ititi"("ti")*
    mangin       = ~r"it(?:r|(?R))*ti"
    thin_lens    = "l"
    dummy        = "d"
    object       = ~r"^d"
    image        = ~r"d$"

.. Created on Sun Sep 27 22:10:01 2023

.. codeauthor: Michael J. Hayford
"""

from parsimonious.grammar import Grammar
from parsimonious.nodes import NodeVisitor

sgz2ele_spec =  \
    r"""
    seq_model   = object space optics image
    optics      = (part space)*

    part        = mangin / cemented / lens / mirror / surface / dummy / thin_lens
    space        = ~r"a|t"

    surface      = "i"
    lens         = "iti"
    mirror       = "r"
    air          = "a"
    cemented     = "ititi"("ti")*
    mangin       = ~r"it(?:r|(?R))*ti"
    thin_lens    = "l"
    dummy        = "d"
    object       = ~r"^d"
    image        = ~r"d$"
    """


sgz2ele_grammar = Grammar(sgz2ele_spec)


def flatten_visit(visit_output):
    elements = []
    _flatten_visit_(elements, visit_output)
    return elements


def _flatten_visit_(elements, visit_output):
    for i,o in enumerate(visit_output):
        if isinstance(o, tuple):
            elements.append(o)
        if isinstance(o, list):
            _flatten_visit_(elements, o)


class SMVisitor(NodeVisitor):

    def __init__(self, do_print_visit=False):
        self.do_print_visit = do_print_visit

    def print_visit(self, node, part_name, idx_list, gap_list):
        if self.do_print_visit:
            # print(f"{part_name[0]}: {part_name[2]} {idx_list} {gap_list}")
            print(f"{part_name} {idx_list} {gap_list} {node.start} {node.end}")

    def visit_seq_model(self, node, visited_children):
        """ Returns the overall output. """
        sg2ele_visit = []
        for child in visited_children:
            sg2ele_visit.append(child)
        return sg2ele_visit

    def _visit_surface_(self, node, part_name):
        """ handle single interface tokens. """
        idx_list = (node.start >> 1,)
        gap_list = ()
        self.print_visit(node, part_name, idx_list, gap_list)
        part_def = part_name, idx_list, gap_list
        return part_def

    def _visit_space_(self, node, part_name):
        """ handle space tokens. """
        idx_list = ()
        gap_list = (node.start >> 1,)
        self.print_visit(node, part_name, idx_list, gap_list)
        part_def = part_name, idx_list, gap_list
        return part_def

    def _visit_element_(self, node, part_name):
        """ handle multi-surface, element, tokens. """
        idx1 = node.start >> 1
        idxk = node.end >> 1
        idx_list = tuple(idx for idx in range(idx1, idxk+1))
        gap_list = tuple(idx for idx in range(idx1, idxk))
        self.print_visit(node, part_name, idx_list, gap_list)
        part_def = part_name, idx_list, gap_list
        return part_def

    def visit_space(self, node, visited_children):
        """ Gets each key/value pair, returns a tuple. """
        space_token = node.full_text[node.start]
        if space_token == 'a':
            part_name = 'air', 'rayoptics.elem.elements', 'AirGap'
        else:
            part_name = 'space', 'rayoptics.elem.elements', 'Space'
        return self._visit_space_(node, part_name)

    def visit_surface(self, node, visited_children):
        """ Gets each key/value pair, returns a tuple. """
        part_name = 'surface', 'rayoptics.elem.elements', 'SurfaceInterface'
        return self._visit_surface_(node, part_name)
    
    def visit_lens(self, node, visited_children):
        """ Gets each key/value pair, returns a tuple. """
        part_name = 'lens', 'rayoptics.elem.elements', 'Element'
        return self._visit_element_(node, part_name)

    def visit_mirror(self, node, visited_children):
        """ Gets each key/value pair, returns a tuple. """
        part_name = 'mirror', 'rayoptics.elem.elements', 'Mirror'
        return self._visit_surface_(node, part_name)

    def visit_air(self, node, visited_children):
        """ Gets each key/value pair, returns a tuple. """
        part_name = 'air', 'rayoptics.elem.elements', 'AirGap'
        return self._visit_space_(node, part_name)

    def visit_cemented(self, node, visited_children):
        """ Gets each key/value pair, returns a tuple. """
        part_name = 'cemented', 'rayoptics.elem.elements', 'CementedElement'
        return self._visit_element_(node, part_name)

    def visit_mangin(self, node, visited_children):
        """ Gets each key/value pair, returns a tuple. """
        part_name = 'mangin', 'rayoptics.elem.elements', 'CementedElement'
        return self._visit_element_(node, part_name)

    def visit_thin_lens(self, node, visited_children):
        """ Gets each key/value pair, returns a tuple. """
        part_name = 'thin_lens', 'rayoptics.elem.elements', 'ThinElement'
        return self._visit_surface_(node, part_name)

    def visit_dummy(self, node, visited_children):
        """ Gets each key/value pair, returns a tuple. """
        part_name = 'dummy', 'rayoptics.elem.elements', 'DummyInterface'
        return self._visit_surface_(node, part_name)

    def visit_object(self, node, visited_children):
        """ Gets each key/value pair, returns a tuple. """
        part_name = 'object', 'rayoptics.elem.elements', 'DummyInterface'
        return self._visit_surface_(node, part_name)

    def visit_image(self, node, visited_children):
        """ Gets each key/value pair, returns a tuple. """
        part_name = 'image', 'rayoptics.elem.elements', 'DummyInterface'
        # idx_list = (-1,)
        # gap_list = ()
        # self.print_visit(node, part_name, idx_list, gap_list)
        # part_def = part_name, idx_list, gap_list
        # return part_def
        return self._visit_surface_(node, part_name)

    def generic_visit(self, node, visited_children):
        """ The generic visit method. """
        return visited_children or node
