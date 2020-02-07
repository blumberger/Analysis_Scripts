#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module that provides useful functions for checking the type of certain variables.

@author: mellis
"""

def eval_type(String):
    """
    Will convert a string to a number if possible. If it isn't return a string.

    Inputs:
        * String  =>  Any string

    Outpus:
        * If the string can be converted to a number it will be with ints being
          preferable to floats. Else will return the same string
    """
    if is_float(String):
        return float(String)
    elif String.isdigit():
        return int(String)
    else:
        return String


def is_float(Str):
    """
    Check whether a string can be represented as a non-integer float

    Inputs:
      * Str <str> => A string to check

    Outputs:
      <bool> Whether the string can be represented as a non-integer float or not.
    """
    if type(Str) == str:
       if is_num(Str):
          if not float(num).is_integer():
             return True
    return False

def is_num(Str):
    """
    Check whether a string can be represented as a number

    Inputs:
      * Str <str> => A string to check

    Outputs:
      <bool> Whether the string can be represented as a number or not.
    """
    try:
        float(Str)
        return True
    except ValueError:
        return False
