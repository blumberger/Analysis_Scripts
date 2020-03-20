#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module that provides useful functions for checking the type of certain variables.
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
    elif is_num(String):
        return int(String)
    else:
        return String


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
          if not float(Str).is_integer():
             return True
    return False


def is_int(num, msg=False):
    """
    Will check whether a parameter is a float or not.

    If a msg is supplied a TypeError will be thrown if the float isn't an int.

    Inputs:
        * num <float> => Number to check
        * msg <str> => The msg to report if it isn't an int
    Outputs:
        <bool>
    """
    if not num.is_integer():
        if msg is not False:
            raise TypeError(msg)
        else:
            return False
    else:
        return True
