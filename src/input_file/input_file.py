#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains all the relevant code that parses the input file that tells
this code what to do with it.

The class 'INP_File' will handle the parsing of the input file, line by line. To
use it call a new instance of it passing the filepath of the input file as an
argument.

IMPORTANT VARIABLES:
    * INP_File    ->    The class that does all the parsing of the input file.

    * CMD_LIST    ->    This lets the input file know what commands are available
                        in the input file. It is a tuple near the top of this file.

    * write_fncs  ->    This is a dictionary that contains a link to all the functions
                        used to write data. The keys give all the file types in the
                        input file that can be written.

    * read_fncs   ->    This is a dictionary that contains a link to all the functions
                        used to read/load data. The keys give all the file types in the
                        input file that can be read.

    * calc_fncs   ->    This is a dictionary that contains a link to all the functions
                        used to calculate data. The keys give all the values that can
                        be specified in the input file to be calculated.
"""

import re
import os

# Fundamental system functions
from src.system import type_checking as type_check

# File loading functions
from src.io_utils import general_io as gen_io
from src.io_utils import CP2K_inp_files as CP2K_inp
from src.io_utils import xyz_files as xyz
from src.io_utils import lammps
from src.io_utils import json_files as json

# Parsing
from src.parsing import general_parsing as gen_parse
from src.parsing import parse_maths

# Calculator functions
from src.calc import pvecs as pvec_lib
from src.calc import NN
from src.calc import density as dens

# Input file functions
from src.input_file import input_file_types as inp_types


CMD_LIST = ('echo', 'write', 'read', 'load', 'calc')


def is_var_line(line):
    """
    Will check if the given line is a variable line or not

    Inputs:
      * line <str> => A string containing the cleaned line from a input file.
    Output
      <bool> True if line is a variable line, else False
    """
    if '=' in line:
        # Check it isn't some other command
        for cmd in CMD_LIST:
            if re.findall(f"^{cmd} ", line):
                return False

        str_txt, non_str = gen_parse.get_str_between_delims(line, '"')
        if any(j in non_str for j in '<>-+/*^'):
            return False
        return True
    else:
        return False


def is_math_line(line):
    """
    Will check if the given line is a line containing maths or not

    Inputs:
      * line <str> => A string containing the cleaned line from a input file.
    Output
      <bool> True if line is a math line, else False
    """
    if '=' in line:
        # Check it isn't some other command
        for cmd in CMD_LIST:
            if re.findall(f"^{cmd} ", line):
                return False

        str_txt, non_str = gen_parse.get_str_between_delims(line, '"')
        if any(j in non_str for j in '<>-+/*^'):
            return True
    return False

###########################################################################
# I'm sure this is terrible practice but I couldn't find a way around it.
###########################################################################
s = "LINE_DECLARATIONS = {"
for cmd in CMD_LIST:
    s += f"'{cmd}': lambda line: len(re.findall('^{cmd} ', line)) > 0,"
s += "}"
exec(s)
###########################################################################

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

    Important Attributes:
        * inp_filepath <str>     =>  The filepath of the input file.
        * file_txt <str>         =>  The file text.
        * file_ltxt <list<str>>  =>  The file text split by "\n" (and cleaned/editted).
        * file_ltxt_orig <list<str>>  =>  The file text split by "\n".

        * line_declarations <dict> => A dictionary of functions to check if a line
                                      is a certain type of line.

        * load_fncs <dict>       =>  Functions/classes to load certain types of
                                     files. Key = type, value = class/function.
                                     These should take 1 argument, the filepath.

        * write_fncs <dict>      =>  Functions/classes to write certain types of
                                     files. Key = type, value = class/function.
                                     All functions/class should take 2 arguments
                                     the data, then the filepath.

        * calc_fncs <dict>       =>  Functions/classes to calculate certain properties
                                     Key = type, value = class/function.

        * variables <list<str>>  =>  A list of all variable names from the file
        * line_nums <list<int>>  =>  A list of the line numbers for each element
                                     in the file_ltxt.
    """
    line_num = 0
    file_ltxt_orig = {}
    variables = []
    E_str = ""   # A variable to count the errors

    load_fncs = {
                 'cp2k_inp': CP2K_inp.parse_inp_file, 'xyz': xyz.XYZ_File,
                 'json': json.read_json, 'lammps_log': lammps.Lammps_Log_File,
                 'txt': gen_io.DataFileStorage, 'lammps_data': lammps.Lammps_Data_File,
                }
    write_fncs = {
                  'cp2k_inp': CP2K_inp.write_inp, 'xyz': xyz.Write_XYZ_File,
                  'json': json.write_json,
                 }
    calc_fncs = {'pvecs': pvec_lib.PVecs, 'NN': NN.NN, 'density': dens.Density}

    line_declarations = LINE_DECLARATIONS
    line_declarations['variable'] = is_var_line
    line_declarations['load'] = lambda x: len(re.findall("^load |^read ", x)) > 0
    line_declarations['math'] = is_math_line

    def __init__(self, inp_filepath):
        self.E_str = "__init__"
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


    #############      Error Checking Methods       #############

    def __check_all_lines(self):
        """
        Will check the lines in the input file for any errors.
        """
        self.E_str = "__check_all_lines"
        variables = []
        for line in self.file_ltxt:
            # Error check any variables
            if self.line_declarations['variable'](line):
                self.__check_variable_line(line)
                name, _ = self.__parse_variable_line__(line)
                variables.append(name)

            # Error check any file loading commands
            elif self.line_declarations['load'](line):
                var = self.__check_load_command(line)
                variables.append(var)

            # Error check any file loading commands
            elif self.line_declarations['write'](line):
                self.__check_write_command(line)

            # Error check any file loading commands
            elif self.line_declarations['math'](line):
                var = self.__check_math_line__(line)
                variables.append(var)

            # Error check any echo commands
            elif self.line_declarations['echo'](line):
                self.__check_echo_command__(line)

            # Error check any calc commands
            elif self.line_declarations['calc'](line):
                var = self.__check_calc_command(line)
                if var != "": variables.append(var)

            self.line_num += 1

        # Reset the inp file variables and line number
        self.line_num = 0
        for var in set(variables):
            delattr(self, var)
            self.variables.remove(var)

    def __check_load_command(self, line):
        """
        Will check a load file command line for errors.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        self.E_str = "__check_load_command"
        err_msg = "The load command takes the syntax 'load <file> <file_type> as <name>'"

        line, any_vars = self.__find_vars_in_line(line)

        words = line.split()
        words = self.__fix_words__(words)
        if len(words) != 5:
            self.__print_error(err_msg, self.line_num, self.E_str)

        try:
            words[1] = type_check.remove_quotation_marks(words[1])
            filepath = gen_io.get_abs_path(words[1])
        except IOError as e:
            self.__print_error(str(e), self.E_str)

        if words[2] not in self.load_fncs:
            err_msg = "I don't know how to load files of type '{words[2]}'."
            err_msg += "\n\nFor a full list of file types that can be loaded see below:\n\t* "
            err_msg += "\n\t* ".join(list(self.load_fncs.keys()))
            self.__print_error(err_msg, self.E_str)

        # Save the variable name for error checking later
        var = inp_types.Variable(words[4], "")
        setattr(self, words[4], var)
        if words[4] not in self.variables:
            self.variables.append(words[4])
        return words[4]

    def __check_write_command(self, line):
        """
        Will check a write file command line for errors.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        self.E_str = "__check_write_command"
        err_msg = "The write command takes the syntax:\n\n\twrite <data_name> <filepath>"
        err_msg += "\n\nor you could specify the type of file to write via:\n\n\t"
        err_msg += "write <data_name> <filepath> as <file_type>"

        line, any_vars = self.__find_vars_in_line(line)
        words = line.split()
        words = self.__fix_words__(words)
        if len(words) != 3 and len(words) != 5:
            self.__print_error(err_msg, self.E_str)

        # Check the variable to be written actually exists
        if words[1] not in self.variables:
            self.__print_error(f"I can't find the data named: '{words[1]}'",
                               self.E_str)

        # Check we know how to write the requested filetype
        if len(words) == 5:
           if words[4] not in self.write_fncs:
               err_msg = "I don't know how to write that type of file.\n\n"
               err_msg += "Please use one of:\n\t*"
               err_msg += "\n\t*".join(list(self.write_fncs.keys()))
               self.__print_error(err_msg, self.E_str)

    def __check_variable_line(self, line):
        """
        Will check a variable setting line for any errors

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        self.E_str = "__check_variable_line"
        line, any_vars = self.__find_vars_in_line(line)
        words = [i for i in line.split('=') if i]
        words = self.__fix_words__(words)

        if len(words) < 2:
            self.__print_error("The syntax for declaring variables is: "
                               + "'<name> = <value>'", self.E_str)

    def __check_math_line__(self, line):
        """
        Will check a math line for any errors

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        self.E_str = "__check_math_line__"
        err_msg = "The syntax for a math command is: math <var> = <arthimetic operation>."
        err_msg += "\nFor example: 'x = x / (1 - x)'\n\n"

        # Check we don't too have many equals signs
        if line.count('=') > 1:
            self.__print_error("Too many '=' found!\n\n" + err_msg, self.E_str)
        elif line.count('=') == 0:
            self.__print_error(f"I can't find a '=' on the math line!{err_msg}",
                                self.E_str)

        # Set the variable for error checking later
        new_var_name, metadata_name = self.__get_variable_name__(line)

        # Split the line by =
        words = line.split('=')
        _, maths = words
        md_var_names, md_names, new_line = self.__parse_metadata_line__(maths,
                                                                        "get")

        # Check metadata and their variables have been declared
        metadata = {}
        for v_name, m_name in zip(md_var_names, md_names):
            if v_name not in self.variables:
                self.__print_error(f"Undeclared variable '{var}'", self.E_str)
            Var = getattr(self, v_name)
            if m_name not in Var.metadata:
                e_msg = f"Undeclared metadata '{m_name}' in variable '{v_name}''"
                self.__print_error(e_msg, self.E_str)
            metadata[m_name] = ""

        # Check if there are any undeclared variables
        line, any_vars = self.__find_vars_in_line(new_line)

        # Check if there are any unwanted characters
        bad_chars = "%£\"!&}{[]}:;@'^~#<,>?¬`|"
        for j in bad_chars:
            if j in new_line:
                self.__print_error(f"Illegal character '{j}' in math",
                                   self.E_str)

        # Check all brackets are closed
        if new_line.count("(") != line.count(")"):
            err_msg = "You've not closed one of the brackets you opened.\n"
            err_msg += "Num of '(' = %i\n" % line.count("(")
            err_msg += "Num of ')' = %i\n" % line.count(")")
            self.__print_error(err_msg, self.E_str)

        Var = inp_types.Variable("", "", {metadata_name: ""})
        setattr(self, new_var_name, Var)
        if new_var_name not in self.variables:
            self.variables.append(new_var_name)

        return new_var_name

    def __check_echo_command__(self, line):
        """
        Will check an echo command line for any errors.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        self.E_str = "__check_echo_command__"

        # Check for any undeclared variables
        self.__find_vars_in_line(line)

        # Check for any undeclared metadata
        line = re.sub("^echo ", "", line)
        md_var_names, md_names, _ = self.__parse_metadata_line__(line, 'get')
        for v_name, m_name in zip(md_var_names, md_names):
            if v_name not in self.variables:
                self.__print_error(f"Undeclared variable '{v_name}'", self.E_str)
            Var = getattr(self, v_name)
            if m_name not in Var.metadata:
                self.__print_error(f"Undeclared metadata '{m_name}'", self.E_str)

    def __check_calc_command(self, line):
        """
        Will check a calc command line for errors.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        self.E_str = "__check_calc_command"
        err_msg =  "The calc command takes the syntax 'calc <property to calc>"
        err_msg += " from <variable name> as <new variable name>'"

        # Clean up the line
        line, any_vars = self.__find_vars_in_line(line)
        words = line.split()
        words = self.__fix_words__(words)

        if len(words) != 6:
            self.__print_error(err_msg, self.E_str)

        _, calc_type, _, var_name, _, new_var_name = words

        # Check the variable to be written actually exists
        if var_name not in self.variables:
            self.__print_error(f"I can't find the data named: '{var_name}'",
                                self.E_str)

        # Check we have a function to calculate the variable
        if calc_type not in self.calc_fncs:
            self.__print_error(f"Can't calculate '{calc_type}' yet choose from"
                                + '\n\t* '
                                + '\n\t* '.join(list(self.calc_fncs.keys())),
                                self.E_str)

        # Check the required metadata has been set
        required_metadata = self.calc_fncs[calc_type].required_metadata
        for attr in required_metadata:
            if attr not in getattr(self, var_name).metadata:
                err_msg = f"'{attr}' required for calculation of '{calc_type}'"
                err_msg += "\n\nPlease set it with the following syntax:\n\t"
                err_msg += f"{var_name}['{attr}'] = <value>"
                self.__print_error(err_msg, self.E_str)

        # Check the required_calc data can be calculated
        required_calc = self.calc_fncs[calc_type].required_calc
        for calc in required_calc:
            if calc not in self.calc_fncs:
                err_msg = f"Calculation of '{calc}' required for calculation of "
                err_msg += f"'{calc_type}'. This currently can't be calculated."
                self.__print_error(err_msg, self.E_str)

        # Set the new var attribute to help error checking later.
        new_var = inp_types.Variable(new_var_name, "")
        setattr(self, new_var_name, new_var)
        if new_var_name not in self.variables:
            self.variables.append(new_var_name)

        return new_var_name


    #############      Parsing Methods       #############

    def __parse_lines(self):
        """
        Will loop over all lines and produce
        """
        self.E_str = "__parse_lines"
        for line in self.file_ltxt:
            if line == "echo": print("")
            # Parse any variables
            if self.line_declarations['variable'](line):
                self.__parse_variable_line__(line)

            # Parse any file loading commands
            elif self.line_declarations['load'](line):
                self.__parse_load_cmd__(line)

            # Parse any file loading commands
            elif self.line_declarations['write'](line):
                self.__parse_write_cmd__(line)

            # Parse any math commands
            elif self.line_declarations['math'](line):
                self.__parse_math_cmd__(line)

            # Parse any echo commands
            elif self.line_declarations['echo'](line):
                self.__parse_echo_cmd__(line)

            # Parse any echo commands
            elif self.line_declarations['calc'](line):
                self.__parse_calc_cmd__(line)

            self.line_num += 1

    def __parse_variable_line__(self, line):
        """
        Will parse a line that sets a variable and save the name in self.variables.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        Outputs:
            (<str>, <str|int|float>) The name, variable.
        """
        self.E_str = "__parse_variable_line__"
        corr_syn = "<name> = <value>      OR      <name>['<metadata'] = <value>"

        # Get the variable name
        var_name, metadata_name = self.__get_variable_name__(line)

        # Get the value
        words = line.split('=')
        value = self.__parse_variable_value__('='.join(words[1:]))

        # If we are setting some metadata
        if metadata_name is not False:
            Var = getattr(self, var_name)
            Var.metadata[metadata_name] = value

        # Just setting a normal variable
        else:
            Var = inp_types.Variable(var_name, value)
            setattr(self, var_name, Var)
            if var_name not in self.variables: self.variables.append(var_name)

        return var_name, Var

    def __get_variable_name__(self, line):
        """
        Will get a variable name from a line in the input file.

        This will also return a metadata name too.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        Outputs:
            <*> The value of a variable`
        """
        self.E_str = "__get_variable_name__"
        corr_syn = "<variable> = <value>   OR   <variable[<metadata>] = <value>"
        words = line.split('=')
        md_var_names, md_names, _ = self.__parse_metadata_line__(words[0],
                                                                get_set="set")
        # Some error checking
        if len(md_var_names) > 1:
            err_msg = "Syntax Error: Can only declare 1 variable per line\n\n"
            self.__print_error(err_msg+corr_syn, self.E_str)
            if md_var_names[0] not in self.variables:
                v_name = md_var_names[0]
                self.__print_error(f"Undeclared variable {v_name}", self.E_str)

        # If we are setting some metadata
        metadata_name = False
        name = words[0].strip()
        if len(md_var_names) == 1:
            name = md_var_names[0].strip()
            metadata_name = md_names[0]

        return name, metadata_name

    def __parse_variable_value__(self, line):
        """
        Will parse the part of the line that is the variable value.

        That is the <value> bit after the equals sign i.e. <var_name> = <value>.

        Inputs:
            * line <str> => A string containing the value from the line (anything after the =)
        Outputs:
            <*> The value of a variable
        """
        self.E_str = "__parse_variable_value__"
        tmp = self.__parse_metadata_line__(line, 'get')
        md_var_names, md_names, md_line = tmp

        # Error Checking
        for v_name, m_name in zip(md_var_names, md_names):
            if v_name not in self.variables:
                self.__print_error(f"Undeclared variable '{v_name}'", self.E_str)
            Var = getattr(self, v_name)
            if m_name not in Var.metadata:
                self.__print_error(f"Undeclared metadata '{m_name}'",
                                   self.E_str)

        # If there is no metadata just set the value
        if len(md_var_names) == 0:
            value = type_check.eval_type(line)
            if type(value) == str:
                value = type_check.remove_quotation_marks(value)
                value, _ = self.__find_vars_in_line(value)

        else:
            # Replace all metadata instances with their values
            for var_i, (v_name, m_name) in enumerate(zip(md_var_names,
                                                         md_names)):
                Var = getattr(self, v_name)
                metadata = Var.metadata[m_name]
                md_line = md_line.replace(f"METADATA_{var_i}", str(metadata))
            value = md_line

        return type_check.eval_type(value)

    def __parse_metadata_line__(self, line, get_set=False):
        """
        Will parse a line that contains some metadata.

        This function will return all the variable names that include some
        metadata syntax, the metadata names and the line with the metadata
        replaced by indices.

        Inputs:
            * line <str> => A string containing the cleaned line from a input
                            file.
        Outputs:
            (<list<str>>, <list<str>>, <str>) var_names, metadata_names, new_line
        """
        self.E_str = "__parse_metadata_line__"
        metadata_regex = "\$[A-Za-z_-]+\['[A-Za-z_-]+'\]|"
        metadata_regex += "\$[A-Za-z_-]+\[\"[A-Za-z_-]+\"\]"
        corr_syn = "<var_name>['<metadata_name>']"

        # First get the all the parts of the line that contain metadata references
        value = False
        str_part, non_str = gen_parse.get_str_between_delims(line, '"')
        metadata_parts = re.findall(metadata_regex, line)
        # Also parse without the dollar sign
        for i in re.findall(metadata_regex.replace("\$", ""), line):
            metadata_parts.append(i)

        # Now get the metadata name and variable name for each reference
        new_line = line
        all_var_names, all_metadata_names = [], []
        for var_i, var in enumerate(metadata_parts):
            new_line = new_line.replace(var, f"METADATA_{var_i}")

            # First get the variable name
            var = var.strip(' $')
            words = var.split('[')
            if len(words) != 2:
                self.__print_error(f"Syntax Error: Correct metadata syntax self.E_strs {corr_syn}", 19)
            var_name = words[0]

            # Is the metadata setting data or getting it?
            if get_set is not False:
                get_set = "get"
                if var_name in non_str.split('=')[0]:
                    get_set = "set"

            # Next get the metadata name
            metadata_name, _ = gen_parse.get_str_between_delims(words[1], "'")
            if not metadata_name:
                metadata_name, _ = gen_parse.get_str_between_delims(words[1], '"')
            if not metadata_name:
                self.__print_error(f"Correct syntax for metadata is: {corr_syn}",
                                   self.E_str)

            # Now get the metadata value
            if var_name not in self.variables:
                self.__print_error(f"Undeclared variable '{var_name}'",
                                   self.E_str)
            Var = getattr(self, var_name)

            err = metadata_name not in Var.metadata and get_set is not False
            err *= get_set == "get"
            if  err:
                err_msg = f"Can't find metadata '{metadata_name}' in variable '{var_name}'"
                err_msg += "\n\n"
                err_msg += f"{var_name} metadata: {Var.metadata}"
                self.__print_error(err_msg, self.E_str)

            all_var_names.append(var_name)
            all_metadata_names.append(metadata_name)

        return all_var_names, all_metadata_names, new_line

    def __parse_math_cmd__(self, line):
        """
        Will parse maths from a variable line.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        Outputs:
            <numeric> The output from the math operations
        """
        self.E_str = "__parse_math_cmd__"
        line = line.strip()

        new_var_name, metadata_name = self.__get_variable_name__(line)

        # Split the line by =
        words = line.split('=')
        _, maths = words
        maths = maths.replace("$", "")
        md_var_names, md_names, new_line = self.__parse_metadata_line__(maths, "get")

        # Create a dictionary of variables to give to the eval_maths func
        variables = {}
        for var_name in self.variables:
            Var = getattr(self, var_name)
            variables[var_name] = Var

        # Add required metadata
        for i, (var_name, md_name) in enumerate(zip(md_var_names, md_names)):
            Var = getattr(self, var_name)
            variables[f"METADATA_{i}"] = Var.metadata[md_name]

        # Actually do the maths
        new_var = parse_maths.eval_maths(new_line, variables)

        # Store the result of the maths
        if metadata_name:
            Var = getattr(self, new_var_name)
            Var.metadata[metadata_name] = new_var
        else:
            setattr(self, new_var_name, new_var)
            if new_var_name not in self.variables:
               self.variables.append(new_var_name)

    def __parse_load_cmd__(self, line):
        """
        Will parse a load file command.

        This function will call the relevant file loading function to parse the
        requested file and store it in the dictionary self.load_data.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        Outputs:
            None
        """
        self.E_str = "__parse_load_cmd__"
        # Split the line by whitespace and get the values of any variables
        line, any_vars = self.__find_vars_in_line(line)
        words = line.split()
        words = self.__fix_words__(words)

        # Read the data
        _, fpath, dtype, _, data_name = words
        fpath = type_check.remove_quotation_marks(fpath)
        fpath = gen_io.get_abs_path(fpath)

        # Create the variable object and save it
        var = inp_types.Variable(data_name, self.load_fncs[dtype](fpath), {'file_type': dtype})
        setattr(self, data_name, var)
        if data_name not in self.variables:
           self.variables.append(data_name)

    def __parse_write_cmd__(self, line):
        """
        Will parse a load file command.

        This function will call the relevant file loading function to parse the
        requested file and store it in the dictionary self.load_data.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        self.E_str = "__parse_write_cmd__"
        # Split the line by whitespace and get the values of any variables
        line, any_vars = self.__find_vars_in_line(line)
        words = line.split()
        words = self.__fix_words__(words)

        # Write the data
        if len(words) == 3:
           _, dname, fpath = words
        elif len(words) == 5:
           _, dname, fpath, _, ftype = words
        fpath = type_check.remove_quotation_marks(fpath)
        fpath = os.path.abspath(os.path.expanduser(fpath))

        # Get the data to be written
        var = getattr(self, dname)
        if len(words) == 3:
           ftype = var.metadata['file_type']

        # Write the data
        self.write_fncs[ftype](var.data, fpath)

    def __parse_echo_cmd__(self, line):
        """
        Will parse an echo command.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        self.E_str = "__parse_echo_cmd__"
        echo_cmd = re.sub("^echo ", "", line).strip()

        # Get any metdata mentioned
        md_var_names, md_names, md_line = self.__parse_metadata_line__(echo_cmd,
                                                                       'get')
        for var_i, (v_name, m_name) in enumerate(zip(md_var_names, md_names)):
            Var = getattr(self, v_name)
            metadata = Var.metadata[m_name]
            md_line = md_line.replace(f"METADATA_{var_i}", str(metadata))
        echo_cmd = md_line

        # Get any variables mentioned
        any_vars = re.findall("\$[A-Za-z_-]+", echo_cmd)
        for check_var in any_vars:
            var = getattr(self, check_var.lstrip('$'))
            echo_cmd = echo_cmd.replace(check_var, str(var))

        # Print the statement
        print(type_check.remove_quotation_marks(echo_cmd))

    def __parse_calc_cmd__(self, line):
        """
        Will parse and run the a calc command

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        self.E_str = "__parse_calc_cmd__"
        # Clean up the line
        line, any_vars = self.__find_vars_in_line(line)
        words = line.split()
        _, calc_type, _, var_name, _, new_var_name = words

        # Get the variable to calculate the property with
        Var = getattr(self, var_name)

        # Initialise the Calc_Obj
        Calc_Obj = self.calc_fncs[calc_type](Var)
        # Calculate and prerequisites and save them in the new object
        for calc in Calc_Obj.required_calc:
            setattr(Calc_Obj, calc, self.calc_fncs[calc](Var).calc())
        Calc_Obj.calc()

        # Create a new variable type
        New_Var = inp_types.Variable(Calc_Obj.name, Calc_Obj, Calc_Obj.metadata)
        setattr(self, new_var_name, New_Var)
        if new_var_name not in self.variables:
            self.variables.append(new_var_name)



    ##### Inp Cleaning methods #################################

    def __clean_inp(self):
        """
        Will loop over the lines and remove comments and whitespace etc...
        """
        self.E_str = "__clean_inp"
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

    def __find_vars_in_line(self, line):
        """
        Will find the variables that have been used in a line of text.

        If any variables found aren't declared an error will be thrown

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        Outputs:
            * (<str>, <list<*>>) The editted line with variables replaced and the variables
        """
        self.E_str = "__find_vars_in_line"
        any_vars = [i[1:] for i in re.findall("\$[A-Za-z_-]+", line)]
        for check_var in any_vars:
            # Check the variable exists
            if check_var not in self.variables:
                self.__print_error(f"Can't find variable '{check_var}'",
                                   self.E_str)

            var = getattr(self, check_var)
            line = line.replace(f"${check_var}", str(var))

        return line, any_vars

    def __fix_words__(self, words):
        """
        Will parse individual words from a line.

        Will remove quotation marks from strings and convert the words to the correct
        type.

        Inputs:
            * words <list<str>> => A list containing the line split by whitespace
        Outputs:
            <list<str|int|float>> The words list but with variables replaced.
        """
        self.E_str = "__fix_words__"
        new_words = []
        for word in words:
            # Get word type and remove quotation marks from strings
            word = word.strip()
            word = type_check.eval_type(word)
            new_word = word

            if type(word) == str:
                word = type_check.remove_quotation_marks(word)

            # Make new words list with old word fixed up
            new_words.append(new_word)

        return new_words



    ############# Error Handling ####################################################

    def __print_error(self, msg, err_str=False, line_num=False,
                      errorFunc=SystemExit):
        """
        Will print an error message in a consistent format.

        Inputs:
            * msg <str> => The error message to be displayed.
            * line_num <int> OPTIONAL => If set: the index of the line. Else it
                                         is self.line_num
            * err_str <str> OPTIONAL => A string to identify which error
                                        message is complaining.
        """
        if line_num is False: line_num = self.line_num
        bad_line_ind = self.line_nums[line_num]

        err_msg = "\n\n\n############  ERROR  #############\n"
        err_msg += "Error in input_file '%s'\n\n---\n" % self.inp_filename
        err_msg += msg.strip("\n")
        err_msg += "\n---\n\nline number: %i\n" % self.line_nums[line_num]
        err_msg += f"line: '{self.file_ltxt_orig[bad_line_ind]}'"
        if err_str is not False:
            err_msg += "\n"
            err_msg += f"err id: {err_str}"
        err_msg += "\n#################################\n\n"
        raise errorFunc(err_msg)
