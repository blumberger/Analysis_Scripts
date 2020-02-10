#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module containing methods relevant to general input and output operations.
"""

import os

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
   if os.path.isfolder(filepath) or '/' not in filepath:
      return filepath

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
