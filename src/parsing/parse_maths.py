#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module created to hold some general parsing functions such as some used to
remove comments from lines in files or get the index of a close of bracket or
to get a string between delimeter (quotation marks).
"""
import re

from src.system import type_checking as type_check
from src.parsing import general_parsing as gen_parse

math_chars = '-+=()/^*'
operator_chars = '^/*+-'  # need to be ordered to satisfy BIDMAS

def parse_math_expressions(txt):
    """
    Will parse mathematical expressions into a list with each element being a
    mathematical object.

    Inputs:
        * txt <str> => The string to parse
    Outputs:
        <list<str>> Each element is a separate mathematical object.
    """
    txt = txt.strip()
    all_exp = []  # stores all math expressions
    build_up_str = ''  # Accumulates to get variable names etc..
    i = 0  # Character counter

    # Replace the scientific notation strings (e.g. 1e-6 etc...)
    floats = re.findall('[0-9.]*[0-9]+ *e *-[0-9]+', txt)
    for i, f in enumerate(floats):
        txt = txt.replace(f, f"FLOAT_{i}")

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

        elif char in math_chars:
            if build_up_str:  all_exp.append(build_up_str)
            all_exp.append(char)
            build_up_str = ""
            i += 1
            continue

        build_up_str += char
        i += 1


    if build_up_str: all_exp.append(build_up_str)

    # Error checking
    check_math_exp(all_exp)

    # Add those floats back in
    for i in range(len(floats)):
        ind = all_exp.index(f"FLOAT_{i}")
        all_exp[ind] = floats[i]

    return [type_check.eval_type(i) for i in all_exp]

def check_math_exp(all_exp):
    """
    Will perform some error checking on a parsed mathematical expression.

    These are:
        1) Check the first character is not a mathematical operator (*, +, -, ...)
        2) Check that mathematical operators and variables alternate

    Inputs:
        * all_exp <list<str>> => A parsed list of mathetmatical objects.
    """
    test_exp = all_exp[:]
    for i in range(test_exp.count("(")):
        test_exp.remove('(')
        test_exp.remove(')')

    if test_exp[0] in math_chars and test_exp[0] not in '()':
        raise SystemExit("Invalid mathematical expression, '%s'" % txt)

    for i in range(len(test_exp)-1):
        prev = test_exp[i] in math_chars
        next = test_exp[i + 1] in math_chars
        if next == prev:
            raise SystemExit("Invalid mathematical expression, '%s'" % txt)

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
        close_ind = gen_parse.get_bracket_close(txt[start_ind:])
        new_text = txt[start_ind+1: close_ind+start_ind]
        var_dict['$VAL%i$' % valI] = eval_maths(new_text, var_dict, val)
        txt = txt.replace("(%s)" % new_text, '$VAL%i$' % valI)
        valI += 1

    # Parse mathematical objects into a list
    all_exp = parse_math_expressions(txt)
    # Use BIDMAS to carry out operations
    evaluated_exp=False
    for operator in operator_chars:
        for op_count in range(all_exp.count(operator)):
            # Get the objects to work with
            op_ind = all_exp.index(operator)
            var1, var2 = all_exp[op_ind - 1], all_exp[op_ind + 1]

            # Set the values of the variables
            if type(var1) == str: var1 = var_dict[var1]
            if type(var2) == str: var2 = var_dict[var2]

            # Decide how to manipilate the objects
            if operator == '+':    tmp = var1 + var2
            elif operator == '*':  tmp = var1 * var2
            elif operator == '/':  tmp = var1 / var2
            elif operator == '^':  tmp = var1 ** var2
            elif operator == '-':  tmp = var1 - var2

            all_exp = all_exp[:op_ind-1] + all_exp[op_ind+1:]
            all_exp[op_ind-1] = tmp

    if len(all_exp) > 1: raise SystemExit("Not all arguments parsed in fnc 'eval_maths'")
    if type(all_exp[0]) == str: all_exp[0] = var_dict[all_exp[0]]

    return all_exp[0]
