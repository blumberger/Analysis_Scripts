#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 11:00:42 2018


A module to parse and edit and save inp files used with CP2K.

To load inp files use the function `parsed_inp_files(filepath <str>)`.

To write them with correct indents etc use the function `write_inp(parsed_inp <dict>, filepath)`


@author: mellis
"""
import re
import os
from collections import OrderedDict
import copy


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
   def  __init__(self, line, line_num=False):
      self.line_num = line_num
      self.line = line

      self.comment, self.extra_parameters = "", []
      self.unit, self.parameter, self.value = "", "", ""
      self.whitespace, self.is_include = False, False
      self.is_parameter, self.is_section = False, False
      self.is_section_end, self.is_section_start = False, False

      self.__clean_line()
      self.__parse_line()

   def __clean_line(self):
       """
       Will remove spaces and parse any comments saving the cleaned line into the
       attribute self.edit_line.
       """
       self.edit_line = self.line.strip()
       self.edit_line, self.comment = rm_comment_from_line(self.line)
       self.edit_line = self.edit_line.strip()

   def __parse_line(self):
       """
       Will choose which line parsing function to use in order to parse the inp line.
       """

       if self.edit_line == "" or self.edit_line.isspace():
          self.whitespace = True
          return
       if '&' == self.edit_line[0]:
           self.is_section = True
           self.__parse_section_line()
       elif '@' == self.edit_line[0]:
           if 'include' in self.edit_line.lower():
               self.is_include = True
               self.include_file = self.edit_line.split()[1:]
       else:
           self.is_parameter = True
           self.__parse_paramter_line()


   def __parse_section_line(self):
       """
       Will parse any line starting with the "&" as a section line.
       """
       words = self.edit_line.lstrip("&").strip().split()
       if 'end' in words[0].lower():
          self.is_section_end = True
          self.section = words[1]
          self.extra_parameters = words[2:]
       else:
          self.is_section_start = True
          self.section = words[0]
          self.extra_parameters = words[1:]

   def __parse_paramter_line(self):
       """
       Will parse a parameter line into the variables, self.parameter, self.value
       and self.units.
       """
       words = self.edit_line.split()

       if len(words) == 2:  self.parameter, self.value = words
       elif re.findall("\[[a-zA-Z]+\]", self.edit_line):
            unit_ind = [i for i, word in enumerate(words) if '[' in word and ']' in word]
            if len(unit_ind) == 1: self.unit = words[unit_ind[0]].strip('[]')
            else: self.unit = [words[i].strip('[]') for i in unit_ind]

            self.parameter, self.value = words[0], [words[i] for i in range(1, len(words)) if i not in unit_ind]

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
        * inp_file <list<inp_line>> => The text from the inp file with every line parsed into an inp_line object.
        * section <str>        => The name of the section to find the end of.
        * curr_line <int>      => The index of the line to start looking from.

    Outputs:
        An integer with the line in the inp file that contains the end section.
    """
    for line_num, parsed_line in enumerate(parsed_inp_lines[curr_line+1:]):
        if parsed_line.is_section and parsed_line.is_section_end \
           and parsed_line.section == section:

           return line_num+1+curr_line
    else:
      raise SystemExit("Can't find the end of the section '%s'" % section)


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
        all_lines = [inp_line(i) for i in inp_file]

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
               parse_inp_file(new_inp_file, new_inp_dict, all_lines[line_num+1:end_ind], new_full_data_dict)

               # Skip past the section as it is being parsed by the line above
               line_num += (end_ind - line_num)   # This also ignores the unecessary END ... line

        # Handle the parsing of paramters
        elif curr_line.is_parameter:
             inp_dict[curr_line.parameter.upper()] = curr_line.value
             full_data_dict[curr_line.parameter.upper()] = curr_line

        elif curr_line.is_include:
             inp_dict.setdefault('INCLUDE', []).append(curr_line.include_file)
             full_data_dict.setdefault('INCLUDE', []).append(curr_line)

        line_num += 1

    return {'params': inp_dict, 'full_data': full_data_dict, 'lines': all_lines}


def rm_comment_from_line(line):
    """
    Will remove any comments in a line for parsing.

    Inputs:
        * line   =>  line from input file

    Ouputs:
        The line with comments removed
    """
    poss_comment_strs = ['#','!']
    for s in poss_comment_strs:
        words = line.split(s)
        if len(words) >= 1:
            line = words[0]
            comment = s.join(words[1:])
            break
        else:
            line = ''.join(words[:-1])
            comment = ""

    return line, comment


def write_inp(nested_inp, s="", tab_depth=0, line_depth=0, tab="    ", prev="section"):
    """
    Will create a string containing the inp file.

    Inputs:
        * nested_inp => A nested dictionary containing all parameters (if order needs to be preserved use an ordered dictionary.

    Outputs:
        * A string containing the inp file.
    """
    good_types = (int, float, str)
    for key in nested_inp:
        # A little hack to put some nice whitespace between sections after a block of parameters
        curr = "param"
        if isinstance(nested_inp[key], (dict, type(OrderedDict()))) or key == "":
            curr = "section"
        if prev == "param" and curr == "section":
           s += "\n"
        prev = curr

        # Print a section and then call the function again recursively
        if isinstance(nested_inp[key], (dict, type(OrderedDict()))):
           # Write the extra parameters for the section
           if '' in nested_inp[key] and nested_inp[key]['']:
               section = "  ".join([key, '  '.join(nested_inp[key][''])])
               if type(nested_inp[key]['']) == str:
                  section = "  ".join([key, nested_inp[key]['']])
           else: section = key
           section = section.strip()

           # Write the section
           new_depth = max([len(i) + len(tab) for i in nested_inp[key]])
           tabs = tab * tab_depth
           if '@IF' in section.upper():
              s += "%s%s\n" % (tabs, section)
              s = write_inp(nested_inp[key], s, tab_depth+1, line_depth=new_depth, prev=prev)
              s = s.rstrip("\n")
              s += "\n%sEND%s\n\n" % (tabs, key.lstrip('@'))
           else:
              s += "%s&%s\n" % (tabs, section)
              s = write_inp(nested_inp[key], s, tab_depth+1, line_depth=new_depth, prev=prev)
              s = s.rstrip("\n")
              s += "\n%s&END %s\n\n" % (tabs, key)

        # Write any parameters within a section
        elif key != "":
           tabs = tab * tab_depth
           if isinstance(nested_inp[key], list):
              unit = ""
              if len(nested_inp[key]) == 2:
                 if nested_inp[key][1]: unit = "%s[%s]" % ("  ", nested_inp[key][1])
              keyStr = "%s%s" % (key.strip(), unit)
              keyStr = keyStr.ljust(line_depth)

              s += "%s%s%s%s\n" % (tabs, keyStr, tab, nested_inp[key][0])

           if isinstance(nested_inp[key], good_types):
              keyStr = key.ljust(line_depth)
              s += "%s%s%s%s\n" % (tabs, keyStr, tab, str(nested_inp[key]))

    return s


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
