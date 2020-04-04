#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module containing methods relevant to general input and output operations.
"""

import os
import copy


class DataFileStorage(object):
	"""
	A class to act as a template for other classes that load and store data from files to use.

	The class takes the filepath as an input and calls a method named parse(). Any data stored
	under the name numeric_data will be manipulated via the operator overload functions.

	Inputs:
		* filepath <str> => The path to the file to be loaded.
	"""
	numeric_data = 0
	_poss_num_types_ = ('numeric_data', 'xyz_data', "csv_data")
	_numeric_data_types = []
	_defaults = {}
	metadata = {'file_type': 'txt'}
	def __init__(self, filepath):
		# Copy all dictionaries to make them unique for each instance
		all_vars = [i for i in dir(self) if i[0] != '_']
		for i in all_vars:
			if i[0] != '_':
				var = getattr(self, i)
				if not callable(var) and isinstance(var, (dict, list, tuple)):
					setattr(self, i, copy.deepcopy(var))

		self.filepath = filepath
		self.file_txt = open_read(self.filepath)

		# Set the default parameters
		for key in self._defaults:
			if key not in self.metadata:
				self.metadata[key] = self._defaults[key]

		self.parse()
		for var in dir(self):
			if var in self._poss_num_types_:
				self._numeric_data_types.append(var)


	# Dummy method to hold the place of an actual parser later
	def parse(self):
		"""
		Should be overridden in any child classes.
		"""
		pass

	def set_data(self):
		"""
		Should be overwritten, is a dummy function for deciding which data to calc with.
		"""
		pass

	# Overload adding
	def __add__(self, val):
		for i in self._numeric_data_types:
			if i in dir(self):
				att = getattr(self, i)
				att += float(val)
				setattr(self, i, att)
		return self

	# Overload multiplying
	def __mul__(self, val):
		for i in self._numeric_data_types:
			if i in dir(self):
				att = getattr(self, i)
				att *= float(val)
				setattr(self, i, att)
		return self

	# Overload subtracting
	def __sub__(self, val):
		for i in self._numeric_data_types:
			if i in dir(self):
				att = getattr(self, i)
				att -= float(val)
				setattr(self, i, att)
		return self

	# Overload division operator i.e. a / b
	def __truediv__(self, val):
		for i in self._numeric_data_types:
			if i in dir(self):
				att = getattr(self, i)
				att /= float(val)
				setattr(self, i, att)
		return self

	# Overload floor division operator i.e. a // b
	def __floordiv__(self, val):
		for i in self._numeric_data_types:
			if i in dir(self):
				att = getattr(self, i)
				att //= float(val)
				setattr(self, i, att)
		return self

	# Overload the power operator
	def __pow__(self, val):
		for i in self._numeric_data_types:
			if i in dir(self):
				att = getattr(self, i)
				att **= float(val)
				setattr(self, i, att)
		return self

	# str() would return filetxt by default
	def __str__(self):
		return self.file_txt


class Write_File(object):
	"""
	A parent class for other file writing classes.

	The general procedure for writing a file is first create the file text
	(saved as self.file_txt) and then write it as a file to the filepath given
	as an argument.

	The function create_file_str must be set and must return the file text.

	This class can also be used to write general files

	Inputs:
	   * Data_Class <class> => The class containing all the data to be written
	   * filepath <str>     => The path to the file to be written.
	   * extension <str> OPTIONAL   => The file extension. Default is False.
	"""
	def __init__(self, Data_Class, filepath, ext=False):
		self.filepath = filepath
		self.Data = Data_Class
		self.extension = ext

		# Get the str to write to a file
		self.file_txt = self.create_file_str()

		if type(self.file_txt) == list:
			if len(self.file_txt) == 1:
				self.write_single_file(self.file_txt[0], filepath, ext)

			else:
				for i, s in enumerate(self.file_txt):
					self.write_single_file(s, filepath, ext, i)

		else:
			self.write_single_file(self.file_txt, filepath, ext)


	def write_single_file(self, s, filepath, ext, file_num=False):
		"""
		Will write a single file from a string.

		Inputs:
			* s <str> => The text to write
		"""
		# Correct the extension
		if self.extension:
			filepath, user_ext = remove_file_extension(filepath)
			if ext is False: ext = user_ext

			if file_num:
				filepath_with_num = f"{filepath}_{count}.{ext}"
			
				count = 0
				while os.path.isfile(filepath_with_num):
					filepath_with_num = f"{filepath}_{count}.{ext}"
					count += 1
				filepath = filepath_with_num
			
			else:    
				filepath = f"{filepath}.{ext}"

		# Write the file with the correct extension
		with open(filepath, 'w') as f:
			f.write(s)


	def create_file_str(self):
	  """
	  To be overwritten by children to create the str to be written to a file.
  
	  By default we write the str(self.Data)
	  """
	  return str(self.Data)


# Reads a file and closes it
def open_read(filename, throw_error=True):
	"""
	A slightly unnecessary function to open a file and read it safely, with the choice to raise an error.

	Inputs:
		* filename <str> => The path to the file that needs opening
		* throw_error <bool> OPTIONAL => Choose whether to throw an error if the file can't be found. Default: True

	Outputs:
		<str> The file text.
	"""
	if os.path.isfile(filename):
		with open(filename, 'r') as f:
			txt = f.read()
		return txt
	else:
		if throw_error:
			raise SystemExit("The %s file doesn't exist!" % filename)
		return False


# Get folder from filepath
def get_folder_from_filepath(filepath):
	"""
	Will get the folderpath of the directory containing the file given in filepath.
 
	This will not raise an error if the filepath doesn't exist. If the filepath given
	is actually a folderpath the input will be returned.
 
	Inputs:
	  * filepath <str> => A str pointing towards a file.
	Outputs:
	  <str> folderpath
	"""
	if os.path.isdir(filepath) or '/' not in filepath:
		return ""

	folder = filepath[:filepath.rfind('/')]  # Splice up to the last '/'
	return folder

# Get filename from filepath
def get_filename_from_filepath(filepath):
	"""
	Will get the name of the file being pointed to in the filepath.
 
	This will not raise an error if the filepath doesn't exist. If the filepath given
	is actually a folderpath the top directory in that folderpath will be returned. If
	the filepath given doesn't point to any files it will simply be returned.
 
	Inputs:
	  * filepath <str> => A str pointing towards a file.
	Outputs:
	  <str> filename
	"""
	if '/' not in filepath:
		return filepath
 
	filename = filepath[filepath.rfind('/')+1:]  # Splice up to the last '/'
	return filename


def remove_file_extension(filepath):
	"""
	Will remove the file extension from a filepath.

	Inputs:
		* filepath <str> => The filepath which should have the extension removed.
	Outputs:
		(<str>, <str>) The filepath with it's extension removed.
	"""
	if '.' not in filepath:
		return filepath

	# The index of the extension
	ext_ind = filepath.rfind('.')

	# Splice the filepath into extension and filepath
	new_filepath = filepath[:ext_ind]
	ext = filepath[ext_ind+1:]
	return new_filepath, ext


def get_abs_path(filepath):
	"""
	Will check if a file/folder path exists and if it does return it's absolute path.
 
	If the filepath doesn't exists an error IOError will be raised.
 
	Inputs:
	  * filepath <str> => A str pointing towards a file.
	Outputs:
	  <str> filepath
	"""
	filepath = os.path.expanduser(filepath)
	filepath = os.path.abspath(filepath)
 
	if not os.path.isfile(filepath) and not os.path.isdir(filepath):
		raise IOError("Can't find file '%s'" % filepath)
	return filepath
