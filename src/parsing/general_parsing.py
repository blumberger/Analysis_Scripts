#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module created to hold some general parsing functions such as some used to
remove comments from lines in files or get the index of a close of bracket or
to get a string between delimeter (quotation marks).
"""
import re

from src.system import type_checking as type_check


def get_str_between_delims(string, start_delim='"', end_delim=False):
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
    if end_delim is False:  end_delim = start_delim

    start_ind = string.find(start_delim)
    if start_ind == -1:
        return "", string

    end_ind = get_bracket_close(string[start_ind:], start_delim=start_delim,
                                end_delim=end_delim)
    end_ind += start_ind

    txt_within_delim = string[start_ind+1: end_ind]
    txt_without_delim = string[:start_ind+1] + string[end_ind:]
    return txt_within_delim, txt_without_delim

def split_str_by_multiple_splitters(string, splitters):
    """
    Will split a string by multiple values.

    For example if a string was 'a,b.c-d' and the splitters were ',.-' then this
    would return ['a', 'b', 'c']

    Inputs:
        * splitters <iterable> => The values to split the string by.
    Outputs:
        <list<str>> The string split by many splitters
    """
    split_parts = []
    build_up = ""
    for i in string:
        if i in splitters:
            if build_up:
                split_parts.append(build_up)
            build_up = ""
        else:
            build_up += i
    if build_up: split_parts.append(build_up)
    return split_parts

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
    brack_num = 1
    for ichar in range(start_ind+1, len(txt)):
        char = txt[ichar]
        if char == end_delim: brack_num -= 1
        elif char == start_delim: brack_num += 1

        if brack_num == 0:
            return ichar + start_ind

    else:
        if txt[-1] == end_delim:
            return len(txt) - 1
        return -1



def get_nums_in_str(string, blank_is_0=False):
    """
    Will get only the numbers from a string

    Inputs:
        * string <str> => The string to get numbers from
    Outputs:
        <list<float>> The numbers from a string
    """
    all_nums = re.findall("[0-9]+", string)
    if blank_is_0 and len(all_nums) == 0:
        return [0.0]
    else:
        return [float(i) for i in all_nums]


def remove_num_from_str(string):
    """
    Will remove any numbers from a string object.

    Inputs:
        * string <str> => The string to remove numbers from
    Outputs:
        <str> A string without numbers
    """
    all_nums = re.findall("[0-9]", string)
    for i in all_non_nums:
        string = string.replace(i, "")

    return string
