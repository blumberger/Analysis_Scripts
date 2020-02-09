#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module created to hold some general parsing functions such as some used to
remove comments from lines in files or to ...
"""
import re

from src.utils import type_checking as type_check


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


def get_bracket_close(txt):
    """
    Get the close of the bracket ')' in a string.

    Inputs:
        * txt <str> => A string with a bracket to be closed including the opening
                       bracket.
    """
    start_ind = txt.find('(')

    brack_num = 0
    for ichar in range(len(txt[start_ind:])):
        char = txt[ichar]
        if char == '(': brack_num += 1
        elif char == ')': brack_num -= 1

        if brack_num == 0:
            return ichar + start_ind


def parse_math_expressions(txt):
    """
    Will parse mathematical expressions into a list with each element being a
    mathematical object.

    Inputs:
        * txt <str> => The string to parse
    Outputs:
        <list<str>> Each element is a separate mathematical object.
    """
    all_exp = []  # stores all math expressions
    build_up_str = ''  # Accumulates to get variable names etc..
    i = 0  # Character counter
    while i < len(txt):
        char = txt[i]

        # If the character is blank skip it
        if char == " ":
            if build_up_str:  all_exp.append(build_up_str)
            build_up_str = ""
            i += 1
            continue

        elif char == '*' and i < len(txt)-1 and txt[i+1] == '*':
            if build_up_str:  all_exp.append(build_up_str)
            all_exp.append('^')
            build_up_str = ""
            i += 2
            continue

        elif char in '^/*+-':
            if build_up_str:  all_exp.append(build_up_str)
            all_exp.append(char)
            build_up_str = ""
            i += 1
            continue

        build_up_str += char
        i += 1

    if build_up_str: all_exp.append(build_up_str)

    # Check the expression, we need alternate operators and variables.
    if all_exp[0] in '^/*+-':
        raise SystemExit("Invalid mathematical expression, '%s'" % txt)
    for i in range(len(all_exp)-1):
        prev = all_exp[i] in '^/*+-'
        next = all_exp[i + 1] in '^/*+-'
        if next == prev:
            raise SystemExit("Invalid mathematical expression, '%s'" % txt)

    return [type_check.eval_type(i) for i in all_exp]


def eval_maths(txt, var_dict={}, val=False):
    """
    Will evaluate a mathematical expression from a string.

    This function will evaluate a mathematical function (given certain variables
    from the var_dict dict). It will use the BIDMAS order of operations.

    For example if '1 + 3' was entered, 4 would be returned.

    Inputs:
        * txt <str> => The maths to be parsed.
        * var_dict <dict> => The variable values that appear in the maths

    Outputs
        <numeric> The value of the maths after parsing.
    """
    # Have a check for non allowed characters (anything that isn't [a-zA-Z0-9.^/*+-])
    all_variables = re.findall("[a-zA-Z]+", txt)
    if not all_variables:
        return eval(txt)
    # Basic error checking
    if txt.count('(') != txt.count(')'): raise SystemExit("Error: bracket mis-match")

    # Deal with brackets
    valI = 0
    while '(' in txt:
        start_ind = txt.find('(')
        close_ind = get_bracket_close(txt[start_ind:])
        new_text = txt[start_ind+1: close_ind+start_ind]
        var_dict['$VAL%i$' % valI] = eval_maths(new_text, var_dict, val)
        txt = txt.replace("(%s)" % new_text, '$VAL%i$' % valI)
        valI += 1

    # Parse mathematical objects into a list
    all_exp = parse_math_expressions(txt)

    # Use BIDMAS to carry out operations
    evaluated_exp=False
    for operator in '^/*+-':
        for op_count in range(all_exp.count(operator)):
            # Get the objects to work with
            op_ind = all_exp.index(operator)
            var1, var2 = all_exp[op_ind - 1], all_exp[op_ind + 1]

            # Set the values of the variables
            if type(var1) == str: var1 = var_dict[var1]
            if type(var2) == str: var2 = var_dict[var2]

            # Decide how to manipilate the objects
            if operator == '+':    all_exp[op_ind - 1] = var1 + var2
            elif operator == '*':  all_exp[op_ind - 1] = var1 * var2
            elif operator == '/':  all_exp[op_ind - 1] = var1 / var2
            elif operator == '^':  all_exp[op_ind - 1] = var1 ** var2
            elif operator == '-':  all_exp[op_ind - 1] = var1 - var2

            all_exp = all_exp[:op_ind] + all_exp[op_ind+2:]

    if len(all_exp) > 1: raise SystemExit("Not all arguments parsed in fnc 'eval_maths'")

    return all_exp[0]
