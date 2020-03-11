#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 11:00:42 2018

A module to parse and edit and save inp files used with CP2K.

To load inp files use the function `parse_inp_file(filepath <str>)`.

To write them with correct indents etc use the function `write_inp(parsed_inp <dict>, filepath)`
"""
import re
import os
from collections import OrderedDict
import copy

from src.parsing import general_parsing as gen_parse
from src.io_utils import general_io as gen_io

class Read_INP(gen_io.DataFileStorage):
    """
    A class to store loaded INP file data.

    This class is just a wrapper for the parse_inp_file function to allow it to
    be used as a file loading class.

    Inputs:
        * filepath <str> => The path to the file to be loaded.
    """
    metadata = {'file_type': 'CP2K_inp'}
    def parse(self):
        self.all_data = parse_inp_file(self.filepath)

    def __str__(self):
        """
        Overload the str magic function.
        """
        return write_inp(self.all_data)


class Write_INP(gen_io.Write_File):
    """
    Inherit from gen_io.Write_File to write CP2K input file.

    This is basically the same as the gen_io.Write_File class see that for more
    details.

     Inputs:
        * Data_Class <class> => The class containing all the data to be written
        * filepath <str>     => The path to the file to be written.
    """
    def __init__(self, Data_Class, filepath):
        super().__init__(Data_Class, filepath, "inp")


class INP_Line(object):
   """
   A class to parse and store data describing a single line in an inp file.

   All methods are private.

   Inputs:
        * line <str> => The line that needs parsing
        * line_num <int> OPTIONAL => The line number within the inp file.

   Attributes:
        * comment => <str> the comment on a line
        * extra_paramters => <list> Any misc info about a section (i.e. &PRINT HIGH, the HIGH would be extra)
        * unit => <str> the unit the parameter is given in
        * parameter => <str> the parameter name
        * value => <str> the value of the parameter
        * whitespace => <bool> whether the line is just whitespace or not
        * is_include => <bool> if the line is an @INCLUDE line
        * is_parameter => <bool> if the line is an @parameter line
        * is_section => <bool> if the line is a section line
        * is_section_start => <bool> if the line is a section_start line
        * is_section_end => <bool> if the line is a section_end line
   """
   def  __init__(self, line, curr_section, line_num=False):
       self.line = line
       self.curr_section = curr_section
       self.line_num = line_num

       self.comment, self.extra_parameters = "", []
       self.unit, self.parameter, self.value = "", "", ""
       self.whitespace, self.is_include = False, False
       self.is_parameter, self.is_section = False, False
       self.is_section_end, self.is_section_start = False, False
       self.is_set, self.set_txt, self.is_coord = False, "", False


       self.clean_line()
       self.parse_line()

   def clean_line(self):
       """
       Will remove spaces and parse any comments saving the cleaned line into the
       attribute self.edit_line.
       """
       self.edit_line = self.line.strip()
       self.edit_line, self.comment = gen_parse.rm_comment_from_line(self.line)
       self.comment = self.comment.strip()
       self.edit_line = self.edit_line.strip()

   def parse_line(self):
       """
       Will choose which line parsing function to use in order to parse the inp line.
       """
       if self.edit_line == "" or self.edit_line.isspace():
          self.whitespace = True
          return
       if '&' == self.edit_line[0]:
           self.is_section = True
           self.parse_section_line()
       elif '@' == self.edit_line[0]:
           words = self.edit_line.split()
           if 'include' in words[0].lower():
               self.is_include = True
               self.include_file = self.edit_line.split()[1:]
           elif 'set' in words[0].lower():
               self.is_set = True
               self.set_txt = '  '.join(self.edit_line.split()[1:])
           else:
               print("WARNING: I don't understand the line '%s'" % self.line)
       else:
           # Check for coord lines
           regex_check = re.findall("[0-9]+\.[0-9]+", self.edit_line)
           if len(regex_check) == 3:
               non_coord_strs = ('abc ', 'alpha_beta_gamma ',)
               if any(j in self.edit_line.lower() for j in non_coord_strs):
                  self.is_parameter = True
                  self.parse_paramter_line()
               else:
                  self.is_coord = True
                  self.parse_coord_line()

           else:
              self.is_parameter = True
              self.parse_paramter_line()


   def parse_section_line(self):
       """
       Will parse any line starting with the "&" as a section line.
       """
       words = self.edit_line.lstrip("&").strip().split()
       if 'end' in words[0].lower():
          self.is_section_end = True
          self.section = self.curr_section
          if len(words) > 1:
              self.section = words[1]
          if len(words) > 2:
              self.extra_parameters = words[2:]
       else:
          self.is_section_start = True
          self.section = words[0]
          self.extra_parameters = words[1:]

   def parse_paramter_line(self):
       """
       Will parse a parameter line into the variables, self.parameter, self.value
       and self.units.
       """
       words = self.edit_line.split()

       # Get a simple paramter value pair
       if len(words) == 2:  self.parameter, self.value = words

       # Try to find any units
       elif re.findall("\[[a-zA-Z_]+\]", self.edit_line):
            unit_ind = [i for i, word in enumerate(words) if '[' in word and ']' in word]
            if len(unit_ind) == 1:
                self.unit = words[unit_ind[0]].strip('[]')
            else:
                self.unit = [words[i].strip('[]') for i in unit_ind]

            self.parameter = words[0]
            self.value = '   '.join([words[i] for i in range(1, len(words)) if i not in unit_ind])

       # Simply save other types of parameters (i.e. parameters with list values)
       else:
            self.parameter = words[0]
            self.value = "    ".join(words[1:])

   def parse_coord_line(self):
      """
      Will parse a line containing coordinates in the inp file.
      """
      words = self.edit_line.split()

      self.elm_name = words[0]
      self.coords = [float(i) for i in words[1:]]

   def __str__(self):
       """
       Override how the in-built str() acts on this class.
       """
       return self.edit_line


def find_section_end(parsed_inp_lines, section, curr_line):
    """
    Will loop over the lines in an inp_file and return the index of the line that
    has the section end in.

    Inputs:
        * inp_file <list<INP_Line>> => The text from the inp file with every line parsed into an INP_Line object.
        * section <str>        => The name of the section to find the end of.
        * curr_line <int>      => The index of the line to start looking from.

    Outputs:
        An integer with the line in the inp file that contains the end section.
    """
    # First get the start of the section
    for line_num, line in enumerate(parsed_inp_lines[curr_line:]):
        if line.is_section and line.is_section_start and line.section == section:
            curr_line = curr_line + line_num + 1
            break

    # if section == 'KIND':
    #     print(curr_line)
    #     [print(i) for i in parsed_inp_lines[curr_line:]]
    #     raise SystemExit("BREAK")
    # Now find the section end
    # Sect num is used to keep track of how many sections have been started and ended
    sect_num = -1
    for line_num, line in enumerate(parsed_inp_lines[curr_line:]):
        # Start a section +1 and End a section -1
        if line.is_section and line.is_section_end:
           sect_num += 1
        elif line.is_section and line.is_section_start:
           sect_num -= 1

        # When we have ended as many sections as we have started we have found the end
        if sect_num == 0:
            line.section = section
            return line_num + curr_line
    else:
        raise SystemError("Can't find the end of the section '%s'" % section)


def parse_inp_file(inp_file, inp_dict=False, all_lines=False, full_data_dict=False):
    """
    A replacement function to parse large inp files.

    This function relies on a for loop to iterate over all steps then
    everytime we hit a section the function calls itself to parse it.
    The previous function was purely recursive and didn't use a for loop,
    for larger inp files (over ~950 lines) the maximum recursion limit
    is reached to prevent stack overflows.

    Inputs:
        * inp_file <str>  => The path to the file to be parsed.

        other parameters shouldn't be changed, they are used to pass
        data from one recursive iteration to the next.

    Outputs:
        A dictionary containing the each parsed line and the nested
        dictionary containing all settings.
    """

    # Read if the inp_file variable isn't a list.
    if type(inp_file) == str:
        with open(inp_file, 'r') as f:
            inp_file = f.read().split("\n")

    # First parse all the lines and create storage
    if inp_dict is False: inp_dict = {}
    if full_data_dict is False: full_data_dict = {}
    if all_lines is False:
        all_lines, curr_section = [], ""
        for i in inp_file:
            parsed_line = INP_Line(i, curr_section)
            if parsed_line.is_section:
                curr_section = parsed_line.section
            all_lines.append(parsed_line)

    # Loop over all lines
    line_num = 0
    while line_num < len(all_lines):
        curr_line = all_lines[line_num]

        # Handle the parsing of sections
        if curr_line.is_section:
            if curr_line.is_section_start:
                # Create new section dicts
                section = curr_line.section.upper()
                if section not in inp_dict:   # If the section doesn't already exist
                    new_inp_dict = inp_dict.setdefault(section, {})
                    new_full_data_dict = full_data_dict.setdefault(section, {'': curr_line})
                else:                         # If the section is repeated pop it in a list
                    if type(inp_dict[section]) != list: inp_dict[section] = [inp_dict[section]]
                    if type(full_data_dict[section]) != list: full_data_dict[section] = [full_data_dict[section]]
                    inp_dict[section].append({})
                    full_data_dict[section].append({})
                    new_inp_dict = inp_dict[section][-1]
                    new_full_data_dict = full_data_dict[section][-1]

                # Find the end of the section
                end_ind = find_section_end(all_lines, curr_line.section, line_num)
                new_inp_file = inp_file[line_num+1: end_ind]

                # Parse the section
                parse_inp_file(new_inp_file, new_inp_dict,
                               all_lines[line_num+1:end_ind],
                               new_full_data_dict)

                # Skip past the section as it is being parsed by the line above
                line_num += (end_ind - line_num)   # This also ignores the unecessary END ... line

        # Handle the parsing of paramters
        elif curr_line.is_parameter:
            inp_dict[curr_line.parameter.upper()] = curr_line.value
            full_data_dict[curr_line.parameter.upper()] = curr_line

        elif curr_line.is_include:
             inp_dict.setdefault('INCLUDE', []).append(curr_line.include_file)
             full_data_dict.setdefault('INCLUDE', []).append(curr_line)

        else:
            if curr_line.edit_line:
                full_data_dict.setdefault("LINES", []).append(curr_line)

        line_num += 1

    return {'params': inp_dict, 'full_data': full_data_dict, 'lines': all_lines}


def get_max_parameter_len_in_section(lines, curr_line_ind):
    """
    Will get the length of the longest parameter in a section within an inp file.

    This is useful when writing the inp file in a pretty way and allow the file
    to have the parameters all written with the same indentation.

    Inputs:
        * lines <list<INP_Line>> => All lines in the inp file with each one parsed
                                  into an INP_Line object within a list.
        * curr_line_ind <int> => An integer giving the index where the section
                                 line appears in the lines list.
    Ouputs:
        <int> The maximum length of the parameter string in a section of an inp file.
    """
    # Check if we have the correct input type
    if not lines[curr_line_ind].is_section:
        raise SystemError("The input line to 'get_max_parameter_len_in_section' should be a section line!")

    # Get the section name
    curr_section = lines[curr_line_ind].section.upper()
    max_param_len = 0

    # Loop over lines until we hit the end section and record the max_param_len
    for line in lines[curr_line_ind:]:
        if line.is_section_end and line.section.upper() == curr_section:
            break
        if line.is_parameter:
            max_param_len = max([len(line.parameter),
                                 max_param_len])

    return max_param_len


def write_inp(inp_dict, filename=False):
    """
    Will create a string containing the inp file.

    Inputs:
        * inp_dict <dict> => The output of the 'parse_inp_file' function.
        * filename <str> OPTIONAL => A string specifying where to output the inp_file (if given).

    Outputs:
        * A string containing the inp file.
    """
    if type(inp_dict) != dict or 'lines' not in inp_dict:
        msg = "Argument 'inp_dict' must be a dictionary as outputted by the function"
        msg += " 'parse_inp_file'.\n\n"
        msg += "Error: Wrong input to fucntion 'write_inp'. Bad 'inp_dict' argument."
        raise SystemError(msg)

    inp_txt = ""
    indent_level = 0
    prev_line_type = ""

    # Loop over every line and add it to the
    param_len = 0
    for line_num, line in enumerate(inp_dict['lines']):
        indent = "   " * indent_level  # the indentation for each line
        line_txt = ""
        if not line.whitespace:
            # Handle writing the setions
            if line.is_section:
                if line.is_section_start:
                    param_len = get_max_parameter_len_in_section(inp_dict['lines'], line_num)
                    if prev_line_type == "parameter" or prev_line_type == "end_section":
                        line_txt += "\n"

                    line_txt += "%s&%s\t%s" % (indent,
                                              line.section.upper(),
                                              '  '.join(line.extra_parameters))
                    indent_level += 1
                    prev_line_type = "start_section"

                elif line.is_section_end:
                    indent_level -= 1
                    indent = "   " * indent_level # Need to change indent here to take effect on this line
                    line_txt += "%s&END %s\t%s" % (indent,
                                                   line.section.upper(),
                                                   '  '.join(line.extra_parameters))
                    prev_line_type = "end_section"


            # Write any parameters (i.e. MD_TIMESTEP   [fs]   0.1)
            elif line.is_parameter:
                unit = "[%s]" % line.unit if line.unit else ""
                line_txt += "%s%s  %s  %s" % (indent,
                                              line.parameter.ljust(param_len),
                                              unit,
                                              line.value)
                prev_line_type = "parameter"

            # If the line has an @include in write it
            elif line.is_include:
                if prev_line_type == "end_section":
                    line_txt += "\n"
                line_txt += "%s%s    %s" % (indent,
                                            "@INCLUDE".ljust(param_len),
                                            '  '.join(line.include_file))
                prev_line_type = "include"

            # If the line has an @set in write it
            elif line.is_set:
                if prev_line_type == "end_section":
                    line_txt += "\n"
                line_txt += "%s%s    %s" % (indent,
                                            "@SET".ljust(param_len),
                                            line.set_txt)
                prev_line_type = "set"

            elif line.is_coord:
               name = line.elm_name
               coords = [("%.6f" % i).ljust(10) for i in line.coords]
               line_txt += "%s%s     %s    %s    %s" % (indent, name.ljust(6),
                                                        coords[0], coords[1], coords[2])

            comment = "# %s" % line.comment if line.comment else ""
            inp_txt += "%s    %s\n" % (line_txt, comment)

        # If there are any comments then write them
        elif line.comment:
            if prev_line_type == "end_section":
                inp_txt += "\n%s# %s\n" % (indent, line.comment)
            else: inp_txt += "%s# %s\n" % (indent, line.comment)
            prev_line_type = "comment"

    # If there is a filename given then write the file
    if filename:
        with open(filename, "w") as f:
            f.write(inp_txt)

    return inp_txt

def check_section_path(section_path, inp_dict):
   """
   Will check a section path exists with the nested dict structure

   Inputs:
      * section_path => A string that points to the section to remove with '%' splitting each subsection e.g.
                        MOTION%MD would remove the MD subsection within the MOTION section.
      * inp_dict => The nested dict containin the input parameters
   """
   sections = section_path.split("%")

   # Check if the path is specified correctly
   if len(sections) == 1 and sections[0] not in inp_dict:
      print("Can't find the section '%s' from '%s'" % (sections[0], section_path))
      print("INP File keys: ",  inp_dict.keys())
      return False
   else:
      return sections

   new_dict = copy.deepcopy(inp_dict)
   for sect in sections:
      if new_dict.get(sect) is None:
         print("Can't find the section '%s' from '%s'" % (sect, section_path))
         print("INP File keys: ",  inp_dict.keys())
         return False
      new_dict = new_dict[sect]

   return sections


def remove_section(section_path, inp_dict):
   """
   Will remove a section from the nested dictionary.

   Inputs:
      * section_path => A string that points to the section to remove with '%' splitting each subsection e.g.
                        MOTION%MD would remove the MD subsection within the MOTION section.
      * inp_dict => The nested dict containin the input parameters

   Output:
      * The dictionary is changed inplace so there is no output
   """
   section_path = section_path.upper()
   sections = check_section_path(section_path, inp_dict)
   if not sections: return

   # Remove the section
   if len(sections) == 1:
      inp_dict.pop(sections[0])
   else:
      new_dict = inp_dict[sections[0]]
      for sect in sections[1:-1]:
         new_dict = inp_dict[sect]
      new_dict.pop(sections[-1])


def create_section(section_path, inp_dict, section):
   """
   Will create a section in the nested dictionary at a certain point.

   Inputs:
      * section_path => A string that points to the section to remove with '%' splitting each subsection e.g.
                        MOTION%MD would remove the MD subsection within the MOTION section.
      * inp_dict => The nested dict containin the input parameters

   Output:
      * The dictionary is changed inplace so there is no output
   """
   if section_path == "":
      for key in section:
         inp_dict[key] = section[key]
   else:
      section_path = section_path.upper()
      sections = check_section_path(section_path, inp_dict)
      if not sections: return

      # First check if we have a section or a parameter
      new_dict = copy.deepcopy(inp_dict)
      last_section = 0
      for sect in sections:
         if isinstance(new_dict[sect], (dict, type(OrderedDict()))):
            new_dict = new_dict[sect]
            last_section += 1
         else:
            break

      # Now get the section dict to change
      new_dict = inp_dict[sections[0]]
      for sect in sections[1:last_section]:
         new_dict = new_dict[sect]

      # If we want the new section after a parameter
      for key in section:
         new_dict[key] = section[key]
