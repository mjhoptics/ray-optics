#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2017 Michael J. Hayford
""" Support for reading CODE V .seq files

.. codeauthor: Michael J. Hayford
"""

import re


def tokenize_command(cmd):
    tkns = re.findall(r"[^'\"\s]\S*|\".+?\"|'.+?'", cmd)
    tokens = []
    for t in tkns:
        if t[:1] == '"':
            tokens.append(t.strip('"'))
        elif t[:1] == "'":
            tokens.append(t.strip("'"))
        else:
            tokens.append(t)
    return tokens


def next_line(it):
    ln = next(it)
    contIndex = ln.rfind('&')
    if contIndex >= 0:
        return ln[:contIndex] + next_line(it)
    else:
        return ln


def strip_comments(textLine):
    commentIndex = textLine.find('!')
    if commentIndex >= 0:
        return textLine[:commentIndex]
    else:
        return textLine


def read_seq_buffer(inputLines, tokenize=True):
    inputs = []
    it = iter(inputLines)
    while True:
        try:
            ln = next_line(it)
            ln = strip_comments(ln)
            lnList = ln.split(';')
            for line in lnList:
                line = line.strip()
                if len(line) > 0:
                    if tokenize:
                        cmd = tokenize_command(line)
                        inputs.append(cmd)
                    else:
                        inputs.append(line)

        except StopIteration:
            break

    return inputs


def read_seq_file(filename, tokenize=True):
    file = open(filename, 'r')
    inpt = file.read()
    file.close()
    inputLines = inpt.splitlines()
    return read_seq_buffer(inputLines, tokenize)
