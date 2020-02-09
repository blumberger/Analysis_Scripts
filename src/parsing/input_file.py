#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains all the relevant code that parses the input file that tells
this code what to do with it.

The class 'INP_File' will handle the parsing of the input file, line by line. To
use it call a new instance of it passing the filepath of the input file as an
argument.
"""

import re
import os

from src.io_utils import general_io as gen_io
from src.parsing import general_parsing as gen_parse
from src.utils import type_checking as type_check

# File loading functions
from src.io_utils import CP2K_inp_files as CP2K_inp
from src.io_utils import xyz_files as xyz

class INP_File(object):
    """
    The class that parses each line in the input file and decide what each line

    The parser works by first cleaning up the file text by removing any comments
    and whitespace etc...

    The file will then be checked for errors, before acting on any commands so
    the user doesn't need to wait until the very end of the analysis for the code
    to crash and have to start again.

    The method 'parse_lines' is then called. This is the method that does the
    majority of the work. It loops over each line and calls the relevant function
    (depending on what syntax it finds on each line) to parse it.

    Inputs:
        * inp_file <str> => The filepath that points to the file in question.

    Attributes:
        * inp_filepath <str>     =>  The filepath of the input file.
        * file_txt <str>         =>  The file text.
        * file_ltxt <list<str>>  =>  The file text split by "\n" (and cleaned/editted).
        * file_ltxt_orig <list<str>>  =>  The file text split by "\n".
        * line_declarations <dict> => A dictionary of functions to check if a line
                                      is a certain type of line.
        * load_fncs <dict>       =>  Functions/classes to load certain types of
                                     files. Key = type, value = class/function.
        * variables <dict>       =>  A dictionary containing the name and value
                                     of each variable in the input file.
        * line_nums <list<int>>  =>  A list of the line numbers for each element
                                     in the file_ltxt.
        * load_data <dict>       =>  Stores the data loaded by a load command.
    """
    variables = {}
    line_num = 0
    file_ltxt_orig = {}
    load_data = {}
    load_fncs = {'inp': CP2K_inp.parse_inp_file, 'xyz': xyz.read_xyz_file}
    write_fncs = {'inp': CP2K_inp.write_inp}#, 'xyz': xyz.write_xyz_file}
    line_declarations = {'variable': lambda x: '=' in x and not any([j in x for j in ('>', '<')]),
                         'load': lambda x: len(re.findall("^load", x)) > 0,
                         'write': lambda x: len(re.findall("^write", x)) > 0,
                        }
    def __init__(self, inp_filepath):
        self.inp_filepath = inp_filepath
        # Read the file
        self.file_txt = gen_io.open_read(inp_filepath)
        self.file_ltxt = self.file_txt.split("\n")
        self.line_nums = list(range(1, len(self.file_ltxt)+1))
        self.inp_filename = self.inp_filepath[self.inp_filepath.rfind('/')+1:]

        # Get line nums for error messages
        for line_num in self.line_nums:
            self.file_ltxt_orig[line_num] = self.file_ltxt[line_num - 1]

        # Now do the parsing
        self.__clean_inp()
        self.__check_all_lines()
        self.__parse_lines()

    def __clean_inp(self):
        """
        Will loop over the lines and remove comments and whitespace etc...
        """
        del_lines = []
        # Check for unnecessary lines in the file.
        for line_num, line in enumerate(self.file_ltxt):
            edit_line, comment = gen_parse.rm_comment_from_line(line)
            edit_line = edit_line.strip()
            if not edit_line:
                del_lines.append(line_num)
            self.file_ltxt[line_num] = edit_line

        # Remove any unnecessary lines
        for line_num in reversed(del_lines):
            del self.file_ltxt[line_num]
            del self.line_nums[line_num]  # For printing of error messages


    #############      Error Checking Methods       #############

    def __check_all_lines(self):
        """
        Will check the lines in the input file for any errors.
        """
        for line in self.file_ltxt:
            # Error check any variables
            if self.line_declarations['variable'](line):
                self.__check_variable_line(line)
                self.__parse_variable_line(line)

            # Error check any file loading commands
            elif self.line_declarations['load'](line):
                self.__check_load_command(line)

            # Error check any file loading commands
            elif self.line_declarations['write'](line):
                self.__check_write_command(line)

            self.line_num += 1

        self.line_num = 0

    def __check_load_command(self, line):
        """
        Will check a load file command line for errors.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        err_msg = "The load command takes the syntax 'load <file> <file_type> as <name>'"

        words = line.split()
        words = self.__fix_words(words)

        if len(words) != 5:
            self.__print_error(err_msg, self.line_num)

        if words[3] != 'as':
            self.__print_error(err_msg + "\n\nThe 4th word should be 'as'")

        if not os.path.isfile(words[1]):
            self.__print_error("Can't find file '%s'" % words[1])

        if words[2] not in self.load_fncs:
            self.__print_error("I don't know how to load files of type '%s'." % words[2])

    def __check_write_command(self, line):
        """
        Will check a write file command line for errors.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        err_msg = "The write command takes the syntax 'write <data> <filepath>'"

        words = line.split()
        words = self.__fix_words(words)
        if len(words) != 3:
            self.__print_error(err_msg)


    def __check_variable_line(self, line):
        """
        Will check a variable setting line for any errors

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        words = [i for i in line.split('=') if i]
        if len(words) < 2:
            self.__print_error("The syntax for declaring variables is <name> = <value>",
                               line_num)


    #############      Parsing Methods       #############

    def __parse_lines(self):
        """
        Will loop over all lines and produce
        """
        for line in self.file_ltxt:
            # Parse any variables
            if self.line_declarations['variable'](line):
                name, variable = self.__parse_variable_line(line)

            # Parse any file loading commands
            elif self.line_declarations['load'](line):
                self.__parse_load_cmd(line)

            # Parse any file loading commands
            elif self.line_declarations['write'](line):
                self.__parse_write_cmd(line)

            self.line_num += 1

    def __parse_variable_line(self, line):
        """
        Will parse a line that sets a variable and save teh name variable pair in self.variables.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        Outputs:
            (<str>, <str|int|float>) The name, variable.
        """
        words = [i for i in line.split('=') if i]

        name = words[0]
        var = type_check.eval_type('='.join(words[1:]))
        if type(var) == str: var = type_check.remove_quotation_marks(var)
        self.variables[name] = var

        return name, var

    def __parse_load_cmd(self, line):
        """
        Will parse a load file command.

        This function will call the relevant file loading function to parse the
        requested file and store it in the dictionary self.load_data.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        Outputs:
            None
        """
        # Split the line by whitespace and get the values of any variables
        words = line.split()
        words = self.__fix_words(words)

        # Read the data
        _, fpath, dtype, _, data_name = words
        fpath = type_check.remove_quotation_marks(fpath)
        self.load_data[data_name] = (self.load_fncs[dtype](fpath), dtype)

    def __parse_write_cmd(self, line):
        """
        Will parse a load file command.

        This function will call the relevant file loading function to parse the
        requested file and store it in the dictionary self.load_data.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        Outputs:
            None
        """
        # Split the line by whitespace and get the values of any variables
        words = line.split()
        words = self.__fix_words(words)

        # Write the data
        _, dname, fpath = words
        fpath = type_check.remove_quotation_marks(fpath)

        data, dtype = self.load_data[dname]
        self.write_fncs[dtype](data, fpath)


    def __fix_words(self, words):
        """
        Will parse individual words from a line.

        Will replace any words in the line with the variable, remove quotation
        marks from strings and convert the words to the correct type.

        Inputs:
            * words <list<str>> => A list containing the line split by whitespace
        Outputs:
            <list<str|int|float>> The words list but with variables replaced.
        """
        new_words = []
        for word in words:
            # Get word type and remove quotation marks from strings
            word = type_check.eval_type(word)
            if type(word) == str: word = type_check.remove_quotation_marks(word)

            # Replace any variables
            if "$" == word[0]:
                for var in self.variables:
                    if word[1:] in self.variables:
                        new_words.append(self.variables[var])
                        break
                    else:
                        self.__print_error("Can't find variable $'%s'" % word[1:])
            else:
                new_words.append(word)

        return new_words

    def __print_error(self, msg, line_num=False, errorFunc=SystemExit):
        """
        Will print an error message in a consistent format.

        Inputs:
            * msg <str> => The error message to be displayed.
            * line_num <int> OPTIONAL => If set: the index of the line. Else it is self.line_num
        """
        if line_num is False: line_num = self.line_num
        bad_line_ind = self.line_nums[line_num]

        err_msg = "\n\n\n############  ERROR  #############\n"
        err_msg += "Error in input_file '%s'\n\n---\n" % self.inp_filename
        err_msg += msg.strip("\n")
        err_msg += "\n---\n\nline number: %i\n" % self.line_nums[line_num]
        err_msg += "line: '%s'" % self.file_ltxt_orig[bad_line_ind]
        err_msg += "\n#################################\n\n"
        raise errorFunc(err_msg)
