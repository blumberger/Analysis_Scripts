#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 11:00:42 2018

A module to parse and edit and save inp files used with CP2K.

To load inp files use the function `parse_inp_file(filepath <str>)`.

To write them with correct indents etc use the function `write_inp(parsed_inp <dict>, filepath)`
"""
import re, os, copy
from collections import OrderedDict
import numpy as np

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
	def _parse_(self):
		self.data = INP_File(parse_inp_file(self.filepath))

		self.file_data = {gen_io.get_filename_from_filepath(self.filepath): copy.deepcopy(self.data)}

	def __str__(self):
		"""
		Overload the str magic function.
		"""
		return str(self.file_data)

	def write(self, filepath):
		count, orig_fpath = 0, filepath[:]
		while os.path.isfile(filepath) and count < 1000:
			fpath, ext = gen_io.remove_file_extension(orig_fpath)
			filepath = f"fpath_{count}.{ext}"
			count += 1
		
		self.data.write(filepath)

	def append(self, val):
		if type(val) == type(self):
			for line in val.data['lines']:
				self.data['lines'].append(line)

			self.data['params'].update(val.data['params'])
			self.data['full_data'].update(val.data['full_data'])

			self.file_data.update(val.file_data)
		else:
			raise SystemError(f"Currently can't append non {type(self)} to {type(self)}.")


class INP_File(dict):

	def __str__(self):
		return self.write()

	def __repr__(self):
		return self.__str__()

	def __init__(self, *args, **kwargs):
		self.update(*args, **kwargs)

		if not all(j in self.keys() for j in ('lines', 'params', 'full_data')):
			raise SystemError("Incorrect use of INP_File. Please input a dictionary with keys: 'lines', 'params', 'full_data'.")

	def check_section_path(self, section_path):
		"""
		Will check a section path exists with the nested dict structure
	  
		Inputs:
			* section_path <list> => The parameter path in the inp_dict

		Outputs:
			<bool> Whether the section is in the inp dict or not.
		"""
		is_sect, e = True, ""
		# Check if the path is specified correctly
		if len(section_path) == 1 and section_path[0] not in self['params']:
			e = f"Can't find the section '{section_path[0]}' from '{section_path}'"
			e += f"INP File keys: {self['params'].keys()}"
			is_sect = False
	  
		new_dict = copy.deepcopy(self['params'])
		for sect in section_path:
			if new_dict.get(sect) is None:
				e = f"Can't find the section '{sect}' from '{section_path}'"
				e += f"INP File keys: {self['params'].keys()}"
				is_sect = False
				break

			new_dict = new_dict[sect]

		return is_sect, e

	def add(self, section_path, new_val, throw_err=True):
		"""
		Will add a section or parameter in the inp file.

		N.B. Only parameters currently handled. I'll need to ammend line_num indexing
			 for sections

		Inputs:
			* parameter_path <list> => The parameter_path path in the inp_dict
			* new_val <*> => The new value of the paramter

		Output:
			Changes occur in-place
		"""		
		is_section = False
		self.__reset_all_line_nums__()

		new_dict = {i: self['full_data'][i] for i in self['full_data']}
		for val in section_path:
			if type(val) == list:
				raise SystemExit("Can't currently handle the repeated sections.")

			if val in new_dict:
				new_dict = new_dict[val]

			else:
				if val == section_path[-1]:
					add_line_num = new_dict[''].line_num + 1
					New_Line = INP_Line(f"{val}      {new_val}",
										section_path[-2], add_line_num)
					new_dict[val] = New_Line
				else:
					is_section = True
					raise SystemError("Can't currently handle adding sections")
				break

		else:
			# The parameter already exists
			if throw_err:
				raise SystemError(f"Parameter {val} already exists!")
			else:
				return

		self['lines'].insert(add_line_num, New_Line)
		self.__reset_all_line_nums__()

	def remove(self, section_path, throw_err=True):
		"""
		Will remove a section or parameter in the inp file.

		N.B. Only parameters currently handled. I'll need to properly remove lines
			 in self['line'] for sections.

		Inputs:
			* section_path <list> => The section path in the inp_dict

		Output:
			Everything is in-place
		"""
		is_sect, err = self.check_section_path(section_path)
		if not is_sect and throw_err:
			raise SystemError(err)
		elif not is_sect: return

		# Change the inp lines dict and full_data dict
		self.__reset_all_line_nums__()
		new_dict = {i: self['full_data'][i] for i in self['full_data']}
		for key in section_path[:-1]:
			new_dict = new_dict[key]
		line = new_dict[section_path[-1]]
		bad_ind = line.line_num

		if line.is_parameter:
			new_dict.pop(section_path[-1])

			del self['lines'][bad_ind]
			self.__reset_all_line_nums__()

			# Change the params dict too
			new_dict = {i: self['params'][i] for i in self['params']}
			for key in section_path[:-1]:
				new_dict = new_dict[key]
			new_dict.pop(section_path[-1])

		else:
			raise SystemExit("Still need to implement section removal!")
	
	def change_param(self, parameter_path, new_val, throw_err=True):
		"""
		Will change a parameter in the inp file.

		Inputs:
			* parameter_path <list> => The parameter_path path in the inp_dict
			* new_val <*> => The new value of the paramter

		Output:
			Changes occur in-place
		"""
		is_sect, err = self.check_section_path(parameter_path)
		if not is_sect and throw_err:
			raise SystemError(err)
		elif not is_sect: return

		# Change the inp lines dict and full_data dict
		new_dict = {i: self['full_data'][i] for i in self['full_data']}
		for key in parameter_path[:-1]:
			new_dict = new_dict[key]

		# Check if the thing is a paramter
		if not new_dict[parameter_path[-1]].is_parameter:
			raise SystemExit(f"{parameter_path} doesn't point to a paramter.")

		new_dict[parameter_path[-1]].set_param(new_val)
		self['lines'][new_dict[parameter_path[-1]].line_num].set_param(new_val)

		# Change parameter_path dict
		new_dict = {i: self['params'][i] for i in self['params']}
		for key in parameter_path[:-1]:
			new_dict = new_dict[key]
		new_dict[parameter_path[-1]] = new_val

	def write(self, filename=False):
		"""
		Will create a string containing the inp file.

		If a filepath is given this will write the file too.

		Inputs:
			* filename <str> OPTIONAL => A string specifying where to output the inp_file (if given).

		Outputs:
			* A string containing the inp file.
		"""
		inp_txt = ""
		indent_level = 0
		prev_line_type = ""

		# Loop over every line and add it to the
		param_len = 0
		for line_num, line in enumerate(self['lines']):
			indent = "   " * indent_level  # the indentation for each line
			line_txt = ""
			if not line.whitespace:
				# Handle writing the setions
				if line.is_section:
					if line.is_section_start:
						param_len = get_max_parameter_len_in_section(self['lines'], line_num)
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

	def __reset_all_line_nums__(self):
		"""Will reset all the line_numbers for each INP_Line object in the file."""
		for i in range(len(self['lines'])): self['lines'][i].line_num = i

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
	line_num = 0
	def  __init__(self, line, curr_section, line_num):
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
			elif 'if' in words[0].lower():
				self.parse_if_statement(words)
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

	def parse_if_statement(self, words):
		"""
		Will save various variables from the if control statement.

		Inputs:
			* words <list<str>> => The line split by ' ' into words.
		"""
		if 'endif' in words[0].lower():
			self.is_if = True
			self.is_end_if = True
		else:
			self.is_if = True
			self.is_end_if = False
			self.if_var_1 = words[1].strip("${} ")
			self.comparator = words[2].strip()
			self.if_var_2 = words[3].strip()

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
		elif re.findall(r"\[[a-zA-Z_]+\]", self.edit_line):
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

	def set_param(self, new_val):
		"""
		Will change the value of the line to something new.

		Inputs:
			* new_val <*> => The line's new value.

		Outputs:
			Everything happens inplace
		"""
		if not self.is_parameter:
			print(f"Warning tried to change non-parameter '{str(self)}' to '{new_val}'", end=" ")
			print(f"for line '{self.edit_line}'")

		self.value = new_val

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


def create_DECOMP_inp(mol_nums, ats_per_mol, filepath=False):
	"""
	Will write a DECOMP.inp file that includes the molecules that need activating in a
	surface hopping run.

	Inputs:
	  * mol_nums <list> => Molecule numbers to activate (zero-indexed -the python way)
	  * ats_per_mol <int> => Number of atoms per molecule
	  * filepath <str> OPTIONAL => Filepath of the DECOMP.inp file. If not given then 
	  							   a file won't be written.

	Outputs:
		<str> The DECOMP section
	"""
	if len(mol_nums) > 2000:
	  raise SystemExit("Can't write more than 2000 mols.")
	# Create the filetxt
	s = "&ENERGY_DECOMP\n\n    INDEX_MOL_DECOMP "
	s += ' '.join([str(i + 1) for i in mol_nums])
	s += "\n"+f"NUM_ACTIVE_ATOMS {len(mol_nums) * ats_per_mol}" + "\n\n"
	s += "&END ENERGY_DECOMP"

	if filepath:
		# Create the directory (if required)
		dir_ = gen_io.get_folder_from_filepath(filepath)
		if dir_ != "" and not os.path.isdir(dir_):
		  os.makedirs(dir_)

		with open(filepath, "w") as f:
		  f.write(s)

	return s


def create_AOM_include(nmol, active_mols, ats_per_mol, single_mol_AOM, filepath=False):
	"""
	Will write the AOM_COEFF.include file that lets SH run know which mols are active etc...

	Inputs:
	  * mol_mask <array> => A boolean array of len nmol, where True means active and False 
	  						is inactive.
	  * ats_per_mol <int> => Number of atoms per molecule
	  * single_mol_AOM <str> => The AOM coefficients for a single molecule.
	  * filepath <str> OPTIONAL => Filepath of the AOM_COEFF.include file (default is to not write)
	"""
	# Create inactive line
	inactive = "XX  1    0     0.00            0.00000000000\n" * ats_per_mol

	AOM_txt = ""
	for i in range(nmol):
		if i in active_mols:
			AOM_txt += single_mol_AOM
		else:
			AOM_txt += inactive


	# Write if we want to 
	if filepath is not False:
		with open(filepath, "w") as f:
			f.write(AOM_txt)

	return AOM_txt


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
		for line_num, i in enumerate(inp_file):
			parsed_line = INP_Line(i, curr_section, line_num=line_num)
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





# def create_section(section_path, inp_dict, section):
# 	"""
# 	Will create a section in the nested dictionary at a certain point.

# 	Inputs:
# 		* section_path => A string that points to the section to remove with '%' splitting each subsection e.g.
# 						MOTION%MD would remove the MD subsection within the MOTION section.
# 		* inp_dict => The nested dict containin the input parameters

# 	Output:
# 		* The dictionary is changed inplace so there is no output
# 	"""
# 	if section_path == "":
# 		for key in section:
# 			inp_dict[key] = section[key]
# 	else:
# 		section_path = section_path.upper()
# 		sections = check_section_path(section_path, inp_dict)
# 		if not sections: return

# 		# First check if we have a section or a parameter
# 		new_dict = copy.deepcopy(inp_dict)
# 		last_section = 0
# 		for sect in sections:
# 			if isinstance(new_dict[sect], (dict, type(OrderedDict()))):
# 				new_dict = new_dict[sect]
# 				last_section += 1
# 			else:
# 				break

# 		# Now get the section dict to change
# 		new_dict = inp_dict[sections[0]]
# 		for sect in sections[1:last_section]:
# 		 new_dict = new_dict[sect]

# 		# If we want the new section after a parameter
# 		for key in section:
# 			new_dict[key] = section[key]

 
  
# def remove_section(section_path, inp_dict):
# 	"""
# 	Will remove a section from the nested dictionary.
 
# 	Inputs:
# 	  * section_path => A string that points to the section to remove with '%' splitting each subsection e.g.
# 						MOTION%MD would remove the MD subsection within the MOTION section.
# 	  * inp_dict => The nested dict containin the input parameters
 
# 	Output:
# 	  * The dictionary is changed inplace so there is no output
# 	"""
# 	section_path = section_path.upper()
# 	sections = check_section_path(section_path, inp_dict)
# 	if not sections: return
 
# 	# Remove the section
# 	if len(sections) == 1:
# 		inp_dict.pop(sections[0])
# 	else:
# 		new_dict = inp_dict[sections[0]]
# 		for sect in sections[1:-1]:
# 			new_dict = inp_dict[sect]
# 		new_dict.pop(sections[-1])
