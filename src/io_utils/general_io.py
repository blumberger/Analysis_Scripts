#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module containing methods relevant to general input and output operations.

@author: oem
"""

import os

# Reads a file and closes it
def open_read(filename, throw_error=True):
    """
    A slightly unnecessary function to open a file and read it safely, with the choice to raise an error.

    Inputs:
        * filename <str> => The path to the file that needs opening
        * throw_error <bool> OPTIONAL => Choose whether to throw an error if the file can't be found. Default: True

    Outputs:
        <str> The file text.
    """
    if os.path.isfile(filename):
        with open(filename, 'r') as f:
            txt = f.read()
        return txt
    else:
        if throw_error:
            raise SystemExit("The %s file doesn't exist!" % filename)
        return False
