#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module created to hold some general parsing functions such as some used to
remove comments from lines in files or get the index of a close of bracket or
to get a string between delimeter (quotation marks).
"""
import re

from src.system import type_checking as type_check


def get_str_between_delims(string, delim='"'):
    """
    Will get the string between 2 delimeters.

    E.g. if a string = 'bob "alice"' this function would return
    ('bob ', 'alice')

    Inputs:
        * string <str> => The txt to search through
        * delim <str> => The delimeter
    Outputs:
        (<str>, <str>) The line without the text within the delimeter and the text within
    """
    start_ind = string.find(delim)
    if start_ind == -1:
        return "", string
    start_ind += 1

    end_ind = get_bracket_close(string[start_ind:], start_delim=delim, end_delim=delim)
    end_ind += start_ind

    txt_within_delim = string[start_ind: end_ind]
    txt_without_delim = string[:start_ind-1] + string[end_ind+1:]

    return txt_within_delim, txt_without_delim



def rm_comment_from_line(line, comment_str='#'):
    """
    Will remove any comments in a line for parsing.

    Inputs:
        * line   =>  line from input file

    Ouputs:
        The line with comments removed
    """
    # Split the line by the comment_str and join the bit after the comment delim
    words = line.split(comment_str)
    if len(words) >= 1:
        line = words[0]
        comment = comment_str.join(words[1:])
    else:
        line = ''.join(words[:-1])
        comment = ""

    return line, comment


def get_bracket_close(txt, start_delim='(', end_delim=')'):
    """
    Get the close of the bracket of a delimeter in a string.

    This will work for nested and non-nested delimeters e.g. "(1 - (n+1))" or
    "(1 - n)" would return the end index of the string.

    Inputs:
        * txt <str> => A string with a bracket to be closed including the opening
                       bracket.
        * delim <str> OPTIONAL => A delimeter (by default it is an open bracket)
    Outputs:
        <int> The index of the corresponding end_delim
    """
    start_ind = txt.find(start_delim)

    brack_num = 0
    for ichar in range(len(txt[start_ind:])):
        char = txt[ichar]
        if char == start_delim: brack_num += 1
        elif char == end_delim: brack_num -= 1

        if brack_num == 0:
            return ichar + start_ind

    else:
        return -1
