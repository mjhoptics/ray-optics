#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2023 Michael J. Hayford
""" Module for parsing a sequential model into elements

The grammar is as follows:

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
    """


sgz2ele_grammar = Grammar(sgz2ele_spec)


def flatten_visit(elements, output):
    for i,o in enumerate(output):
        if isinstance(o, tuple):
            elements.append(o)
        if isinstance(o, list):
            flatten_visit(elements, o)

class SMVisitor(NodeVisitor):
    def print_visit(self, node, part_name, idx_list, gap_list):
        print(f"{part_name[0]}: {part_name[1]} {idx_list} {gap_list}")
#        print(f"{part_name} {idx_list} {gap_list} {node.start} {node.end}")

    def visit_seq_model(self, node, visited_children):
        """ Returns the overall output. """
        output = []
        for child in visited_children:
            output.append(child)
        return output

    def visit_space(self, node, visited_children):
        """ Gets each key/value pair, returns a tuple. """
        idx_list = ()
        gap_list = (node.start >> 1,)
        space_token = node.full_text[node.start]
        part_name = 'space', ('AirGap' if space_token == 'a' else 'Space')
        self.print_visit(node, part_name, idx_list, gap_list)
        part_def = part_name[1], idx_list, gap_list
        return part_def

    def visit_surface(self, node, visited_children):
        """ Gets each key/value pair, returns a tuple. """
        part_name = 'surface', 'SurfaceInterface'
        idx_list = (node.start >> 1,)
        gap_list = ()
        self.print_visit(node, part_name, idx_list, gap_list)
        part_def = part_name[1], idx_list, gap_list
        return part_def
    
    def visit_lens(self, node, visited_children):
        """ Gets each key/value pair, returns a tuple. """
        part_name = 'lens', 'Element'
        idx_list = (idx1 := node.start >> 1, node.end >> 1)
        gap_list = (idx1,)
        self.print_visit(node, part_name, idx_list, gap_list)
        part_def = part_name[1], idx_list, gap_list
        return part_def

    def visit_mirror(self, node, visited_children):
        """ Gets each key/value pair, returns a tuple. """
        part_name = 'mirror', 'Mirror'
        idx_list = (idx := node.start >> 1,)
        gap_list = (idx,)
        self.print_visit(node, part_name, idx_list, gap_list)
        part_def = part_name[1], idx_list, gap_list
        return part_def

    def visit_air(self, node, visited_children):
        """ Gets each key/value pair, returns a tuple. """
        part_name = 'air', 'AirGap'
        idx_list = ()
        gap_list = (node.start >> 1,)
        self.print_visit(node, part_name, idx_list, gap_list)
        part_def = part_name[1], idx_list, gap_list
        return part_def

    def visit_cemented(self, node, visited_children):
        """ Gets each key/value pair, returns a tuple. """
        part_name = 'cemented', 'CementedElement'
        idx1 = node.start >> 1
        idxk = node.end >> 1
        idx_list = (idx for idx in range(idx1, idxk+1))
        gap_list = (idx for idx in range(idx1, idxk))
        self.print_visit(node, part_name, idx_list, gap_list)
        part_def = part_name[1], idx_list, gap_list
        return part_def

    def visit_mangin(self, node, visited_children):
        """ Gets each key/value pair, returns a tuple. """
        part_name = 'mangin', 'CementedElement'
        idx1 = node.start >> 1
        idxk = node.end >> 1
        idx_list = (idx for idx in range(idx1, idxk+1))
        gap_list = (idx for idx in range(idx1, idxk))
        self.print_visit(node, part_name, idx_list, gap_list)
        part_def = part_name[1], idx_list, gap_list
        return part_def

    def visit_thin_lens(self, node, visited_children):
        """ Gets each key/value pair, returns a tuple. """
        part_name = 'thin_lens', 'ThinElement'
        idx_list = (node.start >> 1,)
        gap_list = ()
        self.print_visit(node, part_name, idx_list, gap_list)
        part_def = part_name[1], idx_list, gap_list
        return part_def

    def visit_dummy(self, node, visited_children):
        """ Gets each key/value pair, returns a tuple. """
        part_name = 'dummy', 'DummyInterface'
        idx_list = (node.start >> 1,)
        gap_list = ()
        self.print_visit(node, part_name, idx_list, gap_list)
        part_def = part_name[1], idx_list, gap_list
        return part_def

    def visit_object(self, node, visited_children):
        """ Gets each key/value pair, returns a tuple. """
        part_name = 'object', 'DummyInterface'
        idx_list = (node.start >> 1,)
        gap_list = ()
        self.print_visit(node, part_name, idx_list, gap_list)
        part_def = part_name[1], idx_list, gap_list
        return part_def

    def visit_image(self, node, visited_children):
        """ Gets each key/value pair, returns a tuple. """
        part_name = 'image', 'DummyInterface'
        idx_list = (node.start >> 1,)
        gap_list = ()
        self.print_visit(node, part_name, idx_list, gap_list)
        part_def = part_name[1], idx_list, gap_list
        return part_def

    def generic_visit(self, node, visited_children):
        """ The generic visit method. """
        return visited_children or node
