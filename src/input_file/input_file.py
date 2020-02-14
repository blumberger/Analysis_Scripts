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

        * write_fncs <dict>      =>  Functions/classes to write certain types of
                                     files. Key = type, value = class/function.

        * calc_fncs <dict>       =>  Functions/classes to calculate certain properties
                                     Key = type, value = class/function.

        * variables <list<str>>  =>  A list of all variable names from the file
        * line_nums <list<int>>  =>  A list of the line numbers for each element
                                     in the file_ltxt.
    """
    line_num = 0
    file_ltxt_orig = {}
    variables = []

    load_fncs = {
                 'cp2k_inp': CP2K_inp.parse_inp_file, 'xyz': xyz.XYZ_File,
                 'json': json.read_json, 'lammps_log': lammps.Lammps_Log_File,
                 'txt': gen_io.DataFileStorage, 'lammps_data': lammps.Lammps_Data_File,
                }
    write_fncs = {
                  'cp2k_inp': CP2K_inp.write_inp, 'xyz': xyz.write_xyz_file,
                  'json': json.write_json,
                 }
    calc_fncs = {'pvecs': pvec_lib.PVecs, 'NN': NN.NN}

    line_declarations = LINE_DECLARATIONS
    line_declarations['variable'] = is_var_line
    line_declarations['load'] = lambda x: len(re.findall("^load |^read ", x)) > 0
    line_declarations['math'] = is_math_line

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


    #############      Error Checking Methods       #############

    def __check_all_lines(self):
        """
        Will check the lines in the input file for any errors.
        """
        variables = []
        for line in self.file_ltxt:
            # Error check any variables
            if self.line_declarations['variable'](line):
                self.__check_variable_line(line)
                name, _ = self.__parse_variable_line(line)
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
                var = self.__check_math_line(line)
                if var != "": variables.append(var)

            # Error check any echo commands
            elif self.line_declarations['echo'](line):
                self.__check_echo_command(line)

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
        err_msg = "The load command takes the syntax 'load <file> <file_type> as <name>'"

        line, any_vars = self.__find_vars_in_line(line)
        words = line.split()
        words = self.__fix_words(words)

        if len(words) != 5:
            self.__print_error(err_msg, self.line_num)

        if words[3] != 'as':
            self.__print_error(err_msg + "\n\nThe 4th word should be 'as'")

        try:
           filepath = gen_io.get_abs_path(words[1])
        except IOError as e:
            self.__print_error(str(e))

        if words[2] not in self.load_fncs:
            err_msg = "I don't know how to load files of type '{words[2]}'."
            err_msg += "\n\nFor a full list of file types that can be loaded see below:\n\t* "
            err_msg += "\n\t* ".join(list(self.load_fncs.keys()))
            self.__print_error(err_msg)

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
        err_msg = "The write command takes the syntax:\n\n\twrite <data_name> <filepath>"
        err_msg += "\n\nor you could specify the type of file to write via:\n\n\t"
        err_msg += "write <data_name> <filepath> as <file_type>"

        line, any_vars = self.__find_vars_in_line(line)
        words = line.split()
        words = self.__fix_words(words)
        if len(words) != 3 and len(words) != 5:
            self.__print_error(err_msg)

        # Check the variable to be written actually exists
        if words[1] not in self.variables:
            self.__print_error("I can't find the data named: '%s'" % words[1])

        # Check we know how to write the requested filetype
        if len(words) == 5:
           if words[4] not in self.write_fncs:
               err_msg = "I don't know how to write that type of file.\n\n"
               err_msg += "Please use one of:\n\t*"
               err_msg += "\n\t*".join(list(self.write_fncs.keys()))
               self.__print_error(err_msg)

    def __check_variable_line(self, line):
        """
        Will check a variable setting line for any errors

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        line, any_vars = self.__find_vars_in_line(line)
        words = [i for i in line.split('=') if i]
        words = self.__fix_words(words)

        if len(words) < 2:
            self.__print_error("The syntax for declaring variables is <name> = <value>")

    def __check_math_line(self, line):
        """
        Will check a math line for any errors

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        err_msg = "The syntax for a math command is: math <var> = <arthimetic operation>."
        err_msg += "\nFor example: 'x = x / (1 - x)'\n\n"

        # Check we don't too have many equals signs
        if line.count('=') > 1:
            self.__print_error("Too many '=' found!\n\n" + err_msg)
        elif line.count('=') == 0:
            self.__print_error("I can't find a '=' on the math line!\n\n" + err_msg)

        # Set the variable for error checking later
        new_var, line = line.split("=")
        new_var = new_var.strip()
        if new_var not in self.variables:
            self.variables.append(new_var)
        setattr(self, new_var, "")

        # Check if all the variables are initialised
        any_vars = re.findall("[a-zA-Z_-]+", line)
        any_vars = [i for i in any_vars if i not in '-']
        for var in set(any_vars):
            if var not in self.variables:
                self.__print_error("Can't find variable '%s'" % var)

        # Check if there are any unwanted characters
        bad_chars = "%£\"!&}{[]}:;@'^~#<,>.?¬`|"
        for j in bad_chars:
            if j in line:
                self.__print_error(f"Illegal character '{j}' in mathematical expression.")

        # Check all brackets are closed
        if line.count("(") != line.count(")"):
            err_msg = "You've not closed one of the brackets you opened.\n"
            err_msg += "Num of '(' = %i\n" % line.count("(")
            err_msg += "Num of ')' = %i\n" % line.count(")")
            self.__print_error(err_msg)

        return new_var

    def __check_echo_command(self, line):
        """
        Will check an echo command line for any errors.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        self.__find_vars_in_line(line)

    def __check_calc_command(self, line):
        """
        Will check a calc command line for errors.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        err_msg = "The calc command takes the syntax 'calc <property to calc> from <variable name> as <new variable name>'"

        # Clean up the line
        line, any_vars = self.__find_vars_in_line(line)
        words = line.split()
        words = self.__fix_words(words)

        if len(words) != 6:
            self.__print_error(err_msg)

        _, calc_type, _, var_name, _, new_var_name = words

        # Check the variable to be written actually exists
        if var_name not in self.variables:
            self.__print_error(f"I can't find the data named: '{var_name}'")

        # Check we have a function to calculate the variable
        if calc_type not in self.calc_fncs:
            self.__print_error(f"I cannot yet calculate '{calc_type}'. Please choose from" +
                               '\n\t* ' + '\n\t* '.join(list(self.calc_fncs.keys())))

        # Check the required metadata has been set
        required_metadata = self.calc_fncs[calc_type].required_metadata
        for attr in required_metadata:
            if attr not in getattr(self, var_name).metadata:
               err_msg = f"The attribute '{attr}' is required for the calculation of '{calc_type}'"
               err_msg += "\n\nPlease set it with the following syntax in the input file:\n\t"
               err_msg += f"{var_name}['{attr}'] = <value>"
               self.__print_error(err_msg)

        # Check the required_calc data can be calculated
        required_calc = self.calc_fncs[calc_type].required_calc
        for calc in required_calc:
            if calc not in self.calc_fncs:
               err_msg = f"The attribute '{calc}' is required for the calculation of '{calc_type}'"
               err_msg += "\n\nThis attribute currently can't be calculated so the calculation of"
               err_msg += f"{calc_type} can't be performed."
               self.__print_error(err_msg)

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
        for line in self.file_ltxt:
            # Parse any variables
            if self.line_declarations['variable'](line):
                self.__parse_variable_line(line)

            # Parse any file loading commands
            elif self.line_declarations['load'](line):
                self.__parse_load_cmd(line)

            # Parse any file loading commands
            elif self.line_declarations['write'](line):
                self.__parse_write_cmd(line)

            # Parse any math commands
            elif self.line_declarations['math'](line):
                self.__parse_math_cmd(line)

            # Parse any echo commands
            elif self.line_declarations['echo'](line):
                self.__parse_echo_cmd(line)

            # Parse any echo commands
            elif self.line_declarations['calc'](line):
                self.__parse_calc_cmd(line)

            self.line_num += 1

    def __parse_variable_line(self, line):
        """
        Will parse a line that sets a variable and save teh name variable pair in self.variables.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        Outputs:
            (<str>, <str|int|float>) The name, variable.
        """
        name, metadata_name, var = self.__parse_metadata_indexer(line)

        # If there isn't any metadata then create a variable
        if metadata_name is False:
           # Save the variable
           var = inp_types.Variable(name, var)
           setattr(self, name, var)
           if name not in self.variables:
               self.variables.append(name)
           return name, var

        # Add the metadata to a variable
        else:
            name = name[:name.find('[')]
            if name not in self.variables:
                self.__print_error("Can't find variable '%s'" % var)

            change_var = getattr(self, name)
            change_var.metadata[metadata_name] = var
            return name, var


    def __parse_metadata_indexer(self, line):
        """
        Will parse any metadata names from a variable declaration.

        The metadata in the input file for variables is set by using the syntax:
          `<var_name>['<metadata_name>'] = <value>`

        This function will get the var_name, metadata_name and metadata if there is the
        syntax for metadata setting.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        Outputs:
            (<str>, <str>, <*>) var_name, metadata_name, metadata
        """
        # Clean the line up a bit
        line, any_vars = self.__find_vars_in_line(line)
        words = [i for i in line.split('=') if i]
        words = self.__fix_words(words)

        # Parse the variable
        name = words[0]
        value = '='.join([str(i) for i in words[1:]])
        value = type_check.eval_type(value)
        if type(value) == str: value = type_check.remove_quotation_marks(value)

        # Check we aren't indexing a current variable
        poss_metadata = re.findall("\['[A-Za-z_-]+'\]|\[\"[A-Za-z_-]+\"\]", name)

        # Return the name and value -there is no metadata.
        if len(poss_metadata) == 0:
           return name, False, value

        # Get the metadata name
        elif len(poss_metadata) == 1:
           metadata = type_check.remove_quotation_marks(poss_metadata[0].strip('[]'))
           return name, metadata, value

        # Trying to index metadata currently isn't supported
        else:
            err_msg = "The correct syntax for attributing data to variables is:\n\t`<var_name>['<metadata_name>'] = <value>`"
            self.__print_error(f"Trying to ascribe metadata to metadata currently isn't supported\n\n{err_msg}")

    def __parse_math_cmd(self, line):
        """
        Will parse maths from a variable line.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        Outputs:
            <numeric> The output from the math operations
        """
        line = line.replace("math", "").strip()

        # Split the line by =
        words = line.split('=')
        new_var_name, maths = words
        new_var_name = new_var_name.strip()
        maths = maths.replace("$", "")

        # Create a dictionary of variables to give to the eval_maths func
        variables = {}
        for var_name in self.variables:
            var = getattr(self, var_name)
            variables[var_name] = var

        # Actually do the maths
        new_var = parse_maths.eval_maths(maths, variables)

        # Store the result of the maths
        setattr(self, new_var_name, new_var)
        if new_var_name not in self.variables:
           self.variables.append(new_var_name)

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
        line, any_vars = self.__find_vars_in_line(line)
        words = line.split()
        words = self.__fix_words(words)

        # Read the data
        _, fpath, dtype, _, data_name = words
        fpath = type_check.remove_quotation_marks(fpath)
        fpath = gen_io.get_abs_path(fpath)

        # Create the variable object and save it
        var = inp_types.Variable(data_name, self.load_fncs[dtype](fpath), {'file_type': dtype})
        setattr(self, data_name, var)
        if data_name not in self.variables:
           self.variables.append(data_name)

    def __parse_write_cmd(self, line):
        """
        Will parse a load file command.

        This function will call the relevant file loading function to parse the
        requested file and store it in the dictionary self.load_data.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        # Split the line by whitespace and get the values of any variables
        line, any_vars = self.__find_vars_in_line(line)
        words = line.split()
        words = self.__fix_words(words)

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

    def __parse_echo_cmd(self, line):
        """
        Will parse an echo command.

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
        echo_cmd = line.replace("echo ", "").strip()

        # Get any variables mentioned
        any_vars = re.findall("\$[A-Za-z_-]+", echo_cmd)
        for check_var in any_vars:
            var = getattr(self, check_var.lstrip('$'))
            echo_cmd = echo_cmd.replace(check_var, str(var))

        # Print the statement
        print(type_check.remove_quotation_marks(echo_cmd))

    def __parse_calc_cmd(self, line):
        """
        Will parse and run the a calc command

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        """
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
        Will find the variables that have been used in a line of text

        Inputs:
            * line <str> => A string containing the cleaned line from a input file.
        Outputs:
            * (<str>, <list<*>>) The editted line with variables replaced and the variables
        """
        any_vars = [i[1:] for i in re.findall("\$[A-Za-z_-]+", line)]
        for check_var in any_vars:
            # Check the variable exists
            if check_var not in self.variables:
                self.__print_error("Can't find variable '%s'" % check_var)

            var = getattr(self, check_var)
            line = line.replace(f"${check_var}", str(var))

        return line, any_vars

    def __fix_words(self, words):
        """
        Will parse individual words from a line.

        Will remove quotation marks from strings and convert the words to the correct
        type.

        Inputs:
            * words <list<str>> => A list containing the line split by whitespace
        Outputs:
            <list<str|int|float>> The words list but with variables replaced.
        """
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
