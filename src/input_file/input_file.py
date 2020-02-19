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
from IPython import embed

# Fundamental system functions
from src.system import type_checking as type_check

# File loading functions
from src.io_utils import general_io as gen_io
from src.io_utils import CP2K_inp_files as CP2K_inp
from src.io_utils import xyz_files as xyz
from src.io_utils import lammps
from src.io_utils import json_files as json
from src.io_utils import csv_files

# Parsing
from src.parsing import general_parsing as gen_parse
from src.parsing import parse_maths

# Calculator functions
from src.calc import pvecs as pvec_lib
from src.calc import NN
from src.calc import angular_distribution as ang_dist
from src.calc import density as dens

# Input file functions
from src.input_file import input_file_types as inp_types

CMD_LIST = ('echo', 'write', 'read', 'load', 'calc', 'set', 'shell')
SET_FOLDERPATH = "src/data/set"



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

        * load_fncs <dict>       =>  Classes to load certain types of files.
                                     Key = type, value = class.
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
                 'cp2k_inp': CP2K_inp.Read_INP, 'xyz': xyz.XYZ_File,
                 'json': json.read_json, 'lammps_log': lammps.Lammps_Log_File,
                 'txt': gen_io.DataFileStorage, 'lammps_data': lammps.Lammps_Data_File,
                 'lammps_dump': lammps.Lammps_Dump,
                }
    write_fncs = {
                  'cp2k_inp': CP2K_inp.Write_INP, 'xyz': xyz.Write_XYZ_File,
                  'json': json.write_json, 'csv': csv_files.Write_CSV,
                 }
    calc_fncs = {
                 'pvecs': pvec_lib.PVecs, 'NN': NN.NN, 'density': dens.Density,
                 'angular_dist': ang_dist.Angular_Dist,
                }

    line_declarations = LINE_DECLARATIONS
    line_declarations['variable'] = is_var_line
    line_declarations['load'] = lambda x: len(re.findall("^load |^read ", x)) > 0
    line_declarations['math'] = is_math_line
    line_declarations['shell'] = lambda x: x.strip().lower() == "shell"

    def __init__(self, inp_filepath):
        self.E_str = "init"
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
        self.clean_inp()
        self.check_all_lines()
        self.parse_lines()


    #############      Error Checking Methods       #############

    def check_all_lines(self):
        """
        Will check the lines in the input file for any errors.
        """
        self.E_str = "check_all_lines"
        variables = []
        for line in self.file_ltxt:
            # Error check any variables
            if self.line_declarations['variable'](line):
                self.check_variable_line(line)
                name, _ = self.parse_variable_line(line)
                variables.append(name)

            # Error check any file loading commands
            elif self.line_declarations['load'](line):
                var = self.check_load_command(line)
                variables.append(var)

            # Error check any file loading commands
            elif self.line_declarations['write'](line):
                self.check_write_command(line)

            # Error check any file loading commands
            elif self.line_declarations['math'](line):
                var = self.check_math_line(line)
                variables.append(var)

            # Error check any echo commands
            elif self.line_declarations['echo'](line):
                self.check_echo_command(line)

            # Error check any calc commands
            elif self.line_declarations['calc'](line):
                var = self.check_calc_command(line)
                if var != "": variables.append(var)

            # Error check any calc commands
            elif self.line_declarations['set'](line):
                var = self.check_set_command(line)
                # if var != "": variables.append(var)

            self.line_num += 1

        # Reset the inp file variables and line number
        self.line_num = 0
        for var in set(variables):
            delattr(self, var)
            self.variables.remove(var)

    def check_load_command(self, line):
        """
        Will check a load file command line for errors.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        self.E_str = "check_load_command"
        err_msg = "The load command takes the syntax 'load <file> <file_type> as <name>'"

        line, any_vars = self.find_vars_in_line(line)

        words = line.split()
        words = self.fix_words(words)
        self.E_str = "check_load_command"
        if len(words) != 5:
            self.print_error(err_msg)

        try:
            words[1] = type_check.remove_quotation_marks(words[1])
            filepaths = gen_io.fix_filepath(words[1])
        except IOError as e:
            self.E_str = "check_load_command"
            self.print_error(str(e))

        self.E_str = "check_load_command"
        if words[2] not in self.load_fncs:
            err_msg = "I don't know how to load files of type '{words[2]}'."
            err_msg += "\n\nFor a full list of file types that can be loaded see below:\n\t* "
            err_msg += "\n\t* ".join(list(self.load_fncs.keys()))
            self.print_error(err_msg)

        # Save the variable name for error checking later
        metadata = self.load_fncs[words[2]].metadata
        self.set_var(words[4], "", metadata)
        return words[4]

    def check_write_command(self, line):
        """
        Will check a write file command line for errors.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        self.E_str = "check_write_command"
        err_msg = "The write command takes the syntax:\n\n\twrite <data_name> <filepath>"
        err_msg += "\n\nor you could specify the type of file to write via:\n\n\t"
        err_msg += "write <data_name> <filepath> as <file_type>"

        line, any_vars = self.find_vars_in_line(line)
        words = line.split()
        words = self.fix_words(words)
        if len(words) != 3 and len(words) != 5:
            self.print_error(err_msg)

        # Check the variable to be written actually exists
        if words[1] not in self.variables:
            self.print_error(f"I can't find the data named: '{words[1]}'")

        # Check we know how to write the requested filetype
        if len(words) == 5:
           if words[4] not in self.write_fncs:
               err_msg = "I don't know how to write that type of file.\n\n"
               err_msg += "Please use one of:\n\t*"
               err_msg += "\n\t*".join(list(self.write_fncs.keys()))
               self.print_error(err_msg)

    def check_variable_line(self, line):
        """
        Will check a variable setting line for any errors

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        self.E_str = "check_variable_line"
        line, any_vars = self.find_vars_in_line(line)
        words = [i for i in line.split('=') if i]
        words = self.fix_words(words)

        if len(words) < 2:
            self.print_error("The syntax for declaring variables is: "
                               + "'<name> = <value>'")

    def check_math_line(self, line):
        """
        Will check a math line for any errors

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        self.E_str = "check_math_line"
        err_msg = "The syntax for a math command is: math <var> = <arthimetic operation>."
        err_msg += "\nFor example: 'x = x / (1 - x)'\n\n"

        # Check we don't too have many equals signs
        if line.count('=') > 1:
            self.print_error("Too many '=' found!\n\n" + err_msg)
        elif line.count('=') == 0:
            self.print_error(f"I can't find a '=' on the math line!{err_msg}")

        # Set the variable for error checking later
        new_var_name, metadata_name = self.get_variable_name(line)

        # Split the line by =
        words = line.split('=')
        _, maths = words
        md_var_names, md_names, new_line = self.parse_metadata_line(maths,
                                                                        "get")
        # Check metadata and their variables have been declared
        metadata = {}
        self.E_str = "check_math_line"
        for v_name, m_name in zip(md_var_names, md_names):
            if v_name not in self.variables:
                self.print_error(f"Undeclared variable '{var}'")
            Var = getattr(self, v_name)
            if m_name not in Var.metadata:
                e_msg = f"Undeclared metadata '{m_name}' in variable '{v_name}''"
                self.print_error(e_msg)
            metadata[m_name] = ""

        # Check if there are any undeclared variables
        line, any_vars = self.find_vars_in_line(new_line)
        self.E_str = "check_math_line"

        # Check if there are any unwanted characters
        bad_chars = "%£\"!&}{[]}:;@'^~#<,>?¬`|"
        for j in bad_chars:
            if j in new_line:
                self.print_error(f"Illegal character '{j}' in math")

        # Check all brackets are closed
        if new_line.count("(") != line.count(")"):
            err_msg = "You've not closed one of the brackets you opened.\n"
            err_msg += "Num of '(' = %i\n" % line.count("(")
            err_msg += "Num of ')' = %i\n" % line.count(")")
            self.print_error(err_msg)

        self.set_var(new_var_name, "", {metadata_name: ""})
        return new_var_name

    def check_echo_command(self, line):
        """
        Will check an echo command line for any errors.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        self.E_str = "check_echo_command"

        # Check for any undeclared variables
        self.find_vars_in_line(line)

        # Check for any undeclared metadata
        line = re.sub("^echo ", "", line)
        md_var_names, md_names, _ = self.parse_metadata_line(line, 'get')
        for v_name, m_name in zip(md_var_names, md_names):
            if v_name not in self.variables:
                self.print_error(f"Undeclared variable '{v_name}'")
            Var = getattr(self, v_name)
            if m_name not in Var.metadata:
                self.print_error(f"Undeclared metadata '{m_name}'")

    def check_calc_command(self, line):
        """
        Will check a calc command line for errors.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        self.E_str = "check_calc_command"
        err_msg =  "The calc command takes the syntax 'calc <property to calc>"
        err_msg += " from <variable name> as <new variable name>'"

        # Clean up the line
        line, any_vars = self.find_vars_in_line(line)
        words = line.split()
        words = self.fix_words(words)

        if len(words) != 6:
            self.print_error(err_msg)

        _, calc_type, _, var_name, _, new_var_name = words

        # Check the variable to be written actually exists
        if var_name not in self.variables:
            self.print_error(f"I can't find the data named: '{var_name}'")

        # Check we have a function to calculate the variable
        if calc_type not in self.calc_fncs:
            self.print_error(f"Can't calculate '{calc_type}' yet choose from"
                                + '\n\t* '
                                + '\n\t* '.join(list(self.calc_fncs.keys())))

        # Check the required metadata has been set
        required_metadata = self.calc_fncs[calc_type].required_metadata
        Var = getattr(self, var_name)
        for attr in required_metadata:
            if attr not in Var.metadata:
                err_msg = f"'{attr}' required for calculation of '{calc_type}'"
                err_msg += "\n\nPlease set it with the following syntax:\n\t"
                err_msg += f"{var_name}['{attr}'] = <value>"
                self.print_error(err_msg)

        # Check the required_calc data can be calculated
        required_calc = self.calc_fncs[calc_type].required_calc
        for calc in required_calc:
            if calc not in self.calc_fncs:
                err_msg = f"Calculation of '{calc}' required for calculation of "
                err_msg += f"'{calc_type}'. This currently can't be calculated."
                self.print_error(err_msg)

        # Set the new var attribute to help error checking later.
        self.set_var(new_var_name, "")
        return new_var_name


    def check_set_command(self, line):
        """
        Will check a set command line for errors.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        self.E_str = "check_set_command"
        err_msg =  "The set command takes the syntax 'set <set type> <data name>"
        err_msg += " to <set value>' e.g. 'set system data to pentacene'"

        # Check syntax
        words = line.split()
        if len(words) != 5:
            self.print_error("Too many words found!\n\n" + err_msg)

        # Check undeclared variables
        _, set_folder, var_name, _, set_name = words
        if var_name not in self.variables:
            self.print_error(f"Undeclared variable {var_name}")

        # Run set command for error checking variables later
        self.parse_set_cmd(line)

    #############      Parsing Methods       #############

    def parse_lines(self):
        """
        Will loop over all lines and produce
        """
        self.E_str = "parse_lines"
        for line in self.file_ltxt:
            if line == "echo": print("")

            # Parse any variables
            elif self.line_declarations['variable'](line):
                self.parse_variable_line(line)

            # Parse any file loading commands
            elif self.line_declarations['load'](line):
                self.parse_load_cmd(line)

            # Parse any file loading commands
            elif self.line_declarations['write'](line):
                self.parse_write_cmd(line)

            # Parse any math commands
            elif self.line_declarations['math'](line):
                self.parse_math_cmd(line)

            # Parse any echo commands
            elif self.line_declarations['echo'](line):
                self.parse_echo_cmd(line)

            # Parse any echo commands
            elif self.line_declarations['calc'](line):
                self.parse_calc_cmd(line)

            # Parse any echo commands
            elif self.line_declarations['set'](line):
                self.parse_set_cmd(line)

            # Parse any echo commands
            elif self.line_declarations['shell'](line):
                self.parse_shell_cmd()

            # Print a warning about unknown line
            else:
                self.print_warning("I don't understand a line!")

            self.line_num += 1

    def parse_variable_line(self, line):
        """
        Will parse a line that sets a variable and save the name in self.variables.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        Outputs:
            (<str>, <str|int|float>) The name, variable.
        """
        self.E_str = "parse_variable_line"
        corr_syn = "<name> = <value>      OR      <name>['<metadata'] = <value>"

        # Get the variable name
        var_name, metadata_name = self.get_variable_name(line)

        # Get the value
        words = line.split('=')
        value = self.parse_variable_value('='.join(words[1:]))

        # If we are setting some metadata
        if metadata_name is not False:
            Var = getattr(self, var_name)
            Var[metadata_name] = value

        # Just setting a normal variable
        else:
            Var = self.set_var(var_name, value)

        return var_name, Var

    def get_variable_name(self, line):
        """
        Will get a variable name from a line in the input file.

        This will also return a metadata name too.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        Outputs:
            <*> The value of a variable`
        """
        self.E_str = "get_variable_name"
        corr_syn = "<variable> = <value>   OR   <variable[<metadata>] = <value>"
        words = line.split('=')
        md_var_names, md_names, _ = self.parse_metadata_line(words[0],
                                                                get_set="set")

        # Some error checking
        if len(md_var_names) > 1:
            err_msg = "Syntax Error: Can only declare 1 variable per line\n\n"
            self.print_error(err_msg+corr_syn)
            if md_var_names[0] not in self.variables:
                v_name = md_var_names[0]
                self.print_error(f"Undeclared variable {v_name}")

        # If we are setting some metadata
        metadata_name = False
        name = words[0].strip()
        if len(md_var_names) == 1:
            name = md_var_names[0].strip()
            metadata_name = md_names[0]

        return name, metadata_name

    def parse_variable_value(self, line):
        """
        Will parse the part of the line that is the variable value.

        That is the <value> bit after the equals sign i.e. <var_name> = <value>.

        Inputs:
            * line <str> => A string containing the value from the line (anything after the =)
        Outputs:
            <*> The value of a variable
        """
        tmp = self.parse_metadata_line(line, 'get')
        md_var_names, md_names, md_line = tmp

        # Error Checking
        for v_name, m_name in zip(md_var_names, md_names):
            if v_name not in self.variables:
                self.E_str = "parse_variable_value"
                self.print_error(f"Undeclared variable '{v_name}'")

            Var = getattr(self, v_name)
            if m_name not in Var.metadata:
                self.E_str = "parse_variable_value"
                self.print_error(f"Undeclared metadata '{m_name}'")

        # Split any lists up
        str_part, non_str = gen_parse.get_str_between_delims(line)
        words = [md_line]
        if ',' in non_str:
            words = md_line.split(",")

        # Loop over all values in lists
        values = []
        for word in words:
            # If there is no metadata just set the value
            if len(md_var_names) == 0:
                value = type_check.eval_type(word)
                if type(value) == str:
                    value = type_check.remove_quotation_marks(value)
                    value, _ = self.find_vars_in_line(value)

            # Replace all metadata instances with their values
            else:
                # Loop over all metadata
                for var_i, (v_name, m_name) in enumerate(zip(md_var_names,
                                                             md_names)):
                    Var = getattr(self, v_name)
                    metadata = Var[m_name]
                    word = word.replace(f"METADATA_{var_i}", str(metadata))
                value = word

            value = type_check.eval_type(value)
            values.append(value)

        value = values
        if len(value) == 1:
            value = values[0]

        return value

    def parse_metadata_line(self, line, get_set=False):
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
        self.E_str = "parse_metadata_line"

        metadata_name_regex = "[A-Za-z0-9_-]+"
        var_regex = "\$*[A-Za-z]+[A-Za-z0-9_-]*"
        metadata_regex = f"{var_regex}\['{metadata_name_regex}'\]"
        metadata_regex += f"|{var_regex}\['{metadata_name_regex}'\]"

        corr_syn = "\n\tTo Set: <var_name>['<metadata_name>'] = <value>\n"
        corr_syn = "\n\tTo Get: <var_name>['<metadata_name>']"

        # First get the all the parts of the line that contain metadata references
        value = False
        str_part, non_str = gen_parse.get_str_between_delims(line, '"')
        metadata_parts = re.findall(metadata_regex, line)

        # Now get the metadata name and variable name for each reference
        new_line = line
        all_var_names, all_metadata_names = [], []
        for var_i, var in enumerate(metadata_parts):
            new_line = new_line.replace(var, f"METADATA_{var_i}")

            # First get the variable name
            var = var.strip(' $')
            words = var.split('[')
            if len(words) != 2:
                self.print_error(f"Syntax Error: Correct metadata syntax {corr_syn}")
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
                self.print_error(f"Correct syntax for metadata is: {corr_syn}")

            # Now get the metadata value
            if var_name not in self.variables:
                self.print_error(f"Undeclared variable '{var_name}'")
            Var = getattr(self, var_name)

            err = metadata_name not in Var.metadata and get_set is not False
            err *= get_set == "get"
            if  err:
                err_msg = f"Can't find metadata '{metadata_name}' in variable '{var_name}'"
                err_msg += "\n\n"
                err_msg += f"{var_name} metadata: {Var.metadata}"
                self.print_error(err_msg)

            all_var_names.append(var_name)
            all_metadata_names.append(metadata_name)

        return all_var_names, all_metadata_names, new_line

    def parse_math_cmd(self, line):
        """
        Will parse maths from a variable line.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        Outputs:
            <numeric> The output from the math operations
        """
        self.E_str = "parse_math_cmd"
        line = line.strip()

        new_var_name, metadata_name = self.get_variable_name(line)

        # Split the line by =
        words = line.split('=')
        _, maths = words
        maths = maths.replace("$", "")
        md_var_names, md_names, new_line = self.parse_metadata_line(maths, "get")

        # Create a dictionary of variables to give to the eval_maths func
        variables = {}
        for var_name in self.variables:
            Var = getattr(self, var_name)
            variables[var_name] = Var

        # Add required metadata
        for i, (var_name, md_name) in enumerate(zip(md_var_names, md_names)):
            Var = getattr(self, var_name)
            variables[f"METADATA_{i}"] = Var[md_name]

        # Actually do the maths
        New_Var = parse_maths.eval_maths(new_line, variables)

        # Store the result of the maths
        if metadata_name:
            Var = getattr(self, new_var_name)
            Var[metadata_name] = New_Var
        else:
            setattr(self, new_var_name, New_Var)
            if new_var_name not in self.variables:
               self.variables.append(new_var_name)

    def parse_load_cmd(self, line):
        """
        Will parse a load file command.

        This function will call the relevant file loading function to parse the
        requested file and store it in the dictionary self.load_data.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        Outputs:
            None
        """
        self.E_str = "parse_load_cmd"
        # Split the line by whitespace and get the values of any variables
        line, any_vars = self.find_vars_in_line(line)
        words = line.split()
        words = self.fix_words(words)

        # Read the data
        _, fpath, dtype, _, data_name = words
        fpath = type_check.remove_quotation_marks(fpath)
        fpaths = gen_io.fix_filepath(words[1])

        # Create the variable object and save it
        loaded_data = self.load_fncs[dtype](fpaths)
        metadata = {'file_type': dtype}
        for key in loaded_data.metadata:
            if key not in metadata: metadata[key] = loaded_data.metadata[key]

        self.set_var(data_name, loaded_data, metadata)

    def parse_write_cmd(self, line):
        """
        Will parse a load file command.

        This function will call the relevant file loading function to parse the
        requested file and store it in the dictionary self.load_data.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        self.E_str = "parse_write_cmd"
        # Split the line by whitespace and get the values of any variables
        line, any_vars = self.find_vars_in_line(line)
        words = line.split()
        words = self.fix_words(words)

        # Write the data
        if len(words) == 3:
           _, dname, fpath = words
        elif len(words) == 5:
           _, dname, fpath, _, ftype = words
        fpath = type_check.remove_quotation_marks(fpath)
        fpath = os.path.abspath(os.path.expanduser(fpath))

        # Get the data to be written
        Var = getattr(self, dname)
        if len(words) == 3:
           ftype = Var['file_type']

        # Write the data
        self.write_fncs[ftype](Var.data, fpath)

    def parse_echo_cmd(self, line):
        """
        Will parse an echo command.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        self.E_str = "parse_echo_cmd"
        echo_cmd = re.sub("^echo ", "", line).strip()

        # Get any metdata mentioned
        md_var_names, md_names, md_line = self.parse_metadata_line(echo_cmd,
                                                                       'get')
        for var_i, (v_name, m_name) in enumerate(zip(md_var_names, md_names)):
            Var = getattr(self, v_name)
            metadata = Var[m_name]
            md_line = md_line.replace(f"METADATA_{var_i}", str(metadata))
        echo_cmd = md_line

        # Get any variables mentioned
        any_vars = re.findall("\$[A-Za-z_-]+", echo_cmd)
        for check_var in any_vars:
            var = getattr(self, check_var.lstrip('$'))
            echo_cmd = echo_cmd.replace(check_var, str(var))

        # Print the statement
        print(type_check.remove_quotation_marks(echo_cmd))

    def parse_calc_cmd(self, line):
        """
        Will parse and run the a calc command

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        self.E_str = "parse_calc_cmd"
        # Clean up the line
        line, any_vars = self.find_vars_in_line(line)
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

    def parse_set_cmd(self, line):
        """
        Will set system data to the metadata in a Variable.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        self.E_str = "parse_set_cmd"

        # Get parameters and fix them slightly
        _, set_foldername, var_name, _, set_name = line.split()
        set_name = set_name.replace(" ", "_").lower()
        set_name = f"{set_name}.json"

        # Check folder exists and get its path
        set_folder = os.path.join(SET_FOLDERPATH, set_foldername)
        set_filepath = os.path.join(set_folder, set_name)
        if not os.path.isdir(set_folder):
            set_folders = os.listdir(SET_FOLDERPATH)
            self.print_error(f"Currently can't set '{set_foldername}'."
                             + " Please choose from:\n\t* "
                             + "\n\t* ".join(set_folders)
                             + "\n\n" + "Or create a file at {set_filepath}"
                             )

        # Check the filepath
        if not os.path.isfile(set_filepath):
            set_files = os.listdir(set_folder)
            self.print_error(f"I don't recognise set file '{set_name}'."
                                 + " Please choose from:\n\t* "
                                 + "\n\t* ".join(set_files)
                                 + "\nOr you could create your own in the folder"
                                 + f": {set_folder}.")

        # Read the data
        set_data = json.read_json(set_filepath)

        # Get the variable and add metadata
        Var = getattr(self, var_name)
        for key in set_data:
            Var[key] = set_data[key]

    def parse_shell_cmd(self):
        """
        Will run the IPython shell.
        """
        self.E_str = "parse_set_cmd"

        # Print the welcome message
        print("\n" + "-"*76 + "\n\n\n")
        print("\nWelcome to the IPython Shell.\n")
        print("\nFrom here you will be able to manipulate the variables you"
              + " have declared in your input file. In order to do this it may"
              + " help to know a bit about them first!")
        print("\nEach variable is saved as a 'Variable' type (defined in "
              + " src/input_file/input_file_types.py) and has 3 important "
              + "properties. The data, name and metadata.")
        print("\nThe data attribute:")
        print("The data attribute (accessed via <var>.data) stored the value"
              + " of the data. E.g. for an xyz file it stores the XYZ_File"
              + " object. For the number 2 this would be 2.")
        print("\nThe name attribute:")
        print("The name will give the variable's name (accessed via <var>.name)")
        print("\nThe metadata attribute:")
        print("The metadata attribute (accessed via <var>.metadata) stores the"
              + " metadata for a variable (i.e. <var>['metadata']).")
        print("\n\n\nVariables loaded:")
        for var in self.variables:
            Var = getattr(self, var)
            attrs = [i for i in dir(Var.data) if '_' != i[0] and not callable(getattr(Var.data, i))]
            n_attr = 3
            print("-"*(len(var) + 5))
            print(f"| {var}: |")
            print("-"*(len(var) + 5))
            if attrs:
                print("Attributes:\n")
                for i in range((len(attrs) // n_attr) + 1):
                    print('\t'.join([i.ljust(20) for i in attrs[i*n_attr: (i+1)*n_attr]]))
                print("-" * 76)

        # print("\n" + "-"*76 + "\n\n\n")


        # Declare all the variables in the global scope so the user can use them
        for var_name in self.variables:
            globals()[var_name] = getattr(self, var_name)

        # Open the IPython shell
        embed(colors="Linux")



    ##### Inp Cleaning methods #################################

    def clean_inp(self):
        """
        Will loop over the lines and remove comments and whitespace etc...
        """
        self.E_str = "clean_inp"
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

    def find_vars_in_line(self, line):
        """
        Will find the variables that have been used in a line of text.

        If any variables found aren't declared an error will be thrown

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        Outputs:
            * (<str>, <list<*>>) The editted line with variables replaced and the variables
        """
        self.E_str = "find_vars_in_line"
        any_vars = [i[1:] for i in re.findall("\$[A-Za-z_-]+", line)]
        for check_var in any_vars:
            # Check the variable exists
            if check_var not in self.variables:
                self.print_error(f"Can't find variable '{check_var}'")

            var = getattr(self, check_var)
            line = line.replace(f"${check_var}", str(var))

        return line, any_vars

    def fix_words(self, words):
        """
        Will parse individual words from a line.

        Will remove quotation marks from strings and convert the words to the correct
        type.

        Inputs:
            * words <list<str>> => A list containing the line split by whitespace
        Outputs:
            <list<str|int|float>> The words list but with variables replaced.
        """
        self.E_str = "fix_words"
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

    def fix_filepath(self, filepath):
        """
        Will fix a filepath by expanding a glob and getting the absolute path.
        """
        filepaths = gen_io.expand_glob_filepath(filepath)
        filepaths = [gen_io.get_abs_path(f) for f in filepaths]
        return filepaths


    ############# Utilities ####################################################
    def set_var(self, var_name, var_data, metadata={}):
        """
        Will set a variable in a consistent way.

        Inputs:
            * var_name <str> => Name of the variable to be set
            * var_data <*> => The data to set in the variable
            * metadata <dict> => Any metadata for the variable
        Outputs:
            * <Variable> An instance of the Variable class
        """
        var_name = var_name.strip()

        # Get any old metadata
        if var_name in self.variables:
            Old_Var = getattr(self, var_name)
            md_old = Old_Var.metadata
            for key in md_old:
                if key not in metadata:
                    metadata[key] = md_old[key]

        # Create a new Variable and set it to self
        Var = inp_types.Variable(var_name, var_data, metadata)
        setattr(self, var_name, Var)

        # Append the variable to the list of variables
        if var_name not in self.variables:
            self.variables.append(var_name)

        return Var

    ############# Error Handling ###############################################

    def print_error(self, msg, line_num=False, errorFunc=SystemExit):
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
        err_msg += "\n"
        err_msg += f"err id: {self.E_str}"
        err_msg += "\n#################################\n\n"
        raise errorFunc(err_msg)


    def print_warning(self, msg, line_num=False):
        """
        Will print an warning message in a consistent format.

        Inputs:
            * msg <str> => The error message to be displayed.
            * line_num <int> OPTIONAL => If set: the index of the line. Else it
                                         is self.line_num
        """
        msg = msg.strip('\n')
        if line_num is False: line_num = self.line_num
        bad_line_ind = self.line_nums[line_num]

        warn_msg = "\n_____\nWarning: "
        warn_msg += f"{msg}"
        warn_msg += "\nfile: %s" % self.inp_filename
        warn_msg += "\nline number: %i\n" % self.line_nums[line_num]
        warn_msg += f"line: '{self.file_ltxt_orig[bad_line_ind]}'"
        warn_msg += "\n-----\n\n"
        print(warn_msg)
