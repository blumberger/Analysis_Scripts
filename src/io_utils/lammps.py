#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to store objects that perform I/O tasks on lammps files.
"""

import pandas as pd
import numpy as np
from io import StringIO
import re
import json

from src.io_utils import general_io as gen_io
from src.system import type_checking as type_check


class Lammps_Log_File(gen_io.DataFileStorage):
    """
    Will read the lammps log file.

    A class that will loop over all lines and parse the CSVs from the
    lammps log file.
    """
    data = []
    metadata = {'file_type': 'log_csv'}

    def __init__(self, filepath):
        super().__init__(filepath)

    def _parse(self):
        """
        Will loop over all steps and parse the csv from the lammps log file.
        """
        self.ltxt = self.file_txt.split("\n")
        self.__get_csv_lines__(50)
        self.__read_csv_lines__()


    def __read_csv_lines__(self):
        """
        Will read the csv lines (according to '__get_csv_lines__()') into DataFrames.

        The 2 lists: self.csv_starts and self.csv_ends tell the code where the DataFrame
        starts and ends.
        """
        for ifp, key in enumerate(self.csv_lines):
            # Create file-like object from a string
            fp = StringIO(self.csv_lines[key])
            self.data.append(pd.read_csv(fp, delim_whitespace=True))


    def __get_csv_lines__(self, same_line_tolerance=100):
        """
        Will determine which lines contain data in a csv format.

        This will loop over all lines deciding which ones are csv lines and find the start and the end
        of csv blocks.

        Inputs:
            * same_line_tolerance <int> OPTIONAL => How many lines can have the same length when split by whitespace
                                                    consequtively before being considered part of a csv. Any csv with
                                                    a number of rows below this will be ignored. DEFAULT: 100.
        """
        # Init some vars
        split_num, split_same = 0, 0
        is_csv_line = False, True
        self.csv_starts, self.csv_ends = [], []

        # Find where there are csv lines
        for line_num, line in enumerate(self.ltxt):

            # Check if the previous line has the same format as the current
            new_split_num = len(line.split())
            if split_num == new_split_num:
                split_same += 1
            else:
                split_same = 0

            # If there are more lines with a consistent format with the tolerance then
            #  declare it as a csv.
            if split_same > same_line_tolerance:
                # If the value of is_csv_line is False the previous line mustn't of been a csv line
                if is_csv_line is False:
                    self.csv_starts.append(line_num - 1 - same_line_tolerance)

                is_csv_line = True

            else:
                if is_csv_line is True:
                    self.csv_ends.append(line_num)

                is_csv_line = False

            split_num = new_split_num

        self.csv_line_nums = np.arange(len(self.ltxt))

        # Some error checking
        if len(self.csv_starts) != len(self.csv_ends):
            if len(self.csv_starts) == len(self.csv_ends) - 1:
                self.csv_ends.append(len(ltxt))
            else:
                err_msg = f"Number of csv block starts = {len(self.csv_starts)}"
                err_msg += "\n"
                err_msg += f"Number of csv block ends = {len(self.csv_ends)}"
                err_msg += "\n"
                err_msg += f"Can't find the end of the csv block in the file {self.filepath}"
                raise EOFError(err_msg)

        # Get the lines of the csvs
        self.csv_lines = {i: '\n'.join(self.ltxt[start: end])
                          for i, (start, end) in enumerate(zip(self.csv_starts, self.csv_ends))}


def write_lammps_log_CSVs(Lammps_Log_File, filepath):
    """
    Will write the different CSVs for each csv in the Lammps log file.

    Inputs:
        * Lammps_Log_File <Lammps_Log_File> => An instance of the Lammps_Log_File class.
    """
    for df_num, df in enumerate(Lammps_Log_File.data):
        filename, ext = gen_io.remove_file_extension(filepath)
        print(filename, ext)




class Lammps_Data_File(gen_io.DataFileStorage):
    """
    Will read the lammps log file.

    A class that will loop over all lines and parse the CSVs from the
    lammps log file.
    """
    data = {}
    metadata = {'file_type': 'log_csv'}
    sects = ('masses', 'atoms', 'bonds', 'angles', 'dihedrals',)

    def __init__(self, filepath):
        super().__init__(filepath)

    def _parse(self):
        """
        Will loop over all steps and parse the csv from the lammps log file.
        """
        self.ltxt = [i for i in self.file_txt.split("\n") if i]

        ftxt = self.file_txt.lower()

        # Should check the order of the sections here!

        # Get the sections in the file
        self.__divide_sections__()
        self.__parse_params_sect__()
        self.__parse_masses_sect__()

        headers = ("id", "mol_id", "at_type", "x", "y", "z", "ix", "iy", "iz",)
        self.data['atoms'] = self.__parse_numeric_section__(self.atoms_sect, headers)

        headers = ("id", "at_type", "at1", "at2")
        self.data['bonds'] = self.__parse_numeric_section__(self.bonds_sect, headers)

        headers = ("id", "at_type", "at1", "at2", "at3")
        self.data['angles'] = self.__parse_numeric_section__(self.angles_sect, headers)

        headers = ("id", "at_type", "at1", "at2", "at3", "at4")
        self.data['dihedrals'] = self.__parse_numeric_section__(self.dihedrals_sect, headers)

    def __divide_sections__(self):
        """
        Will divide the file up into sections -i.e. Atoms, Bonds, Masses etc...

        Section 1 gives number of atoms, bonds, dihedrals, angles, atom types etc...
        Section 2 gives masses of each atom type
        Section 3 gives atom coords
        Section 4 gives bonding info
        Section 5 gives angle info
        Section 6 gives dihedral info
        """
        # ftxt is the filetxt without whitespace
        ftxt = "\n".join(self.ltxt).lower()

        # Get first section
        divide = ftxt.lower().split('masses')
        self.__check_len_sect__('masses', divide)    # error checking
        ind = self.__search_in_list_of_str__(divide, lambda s: "xlo xhi" in s, 1)[0]
        self.params_sect = divide[ind]

        # Get masses section
        divide = divide[1-ind]
        divide = divide.split('atoms')
        self.__check_len_sect__('Masses', divide)     # error checking
        # Find the section with 2 columns
        self.masses_sect, divide = self.__get_numeric_section__(divide, 2)

        # Get atoms section
        divide = divide.split("bonds")
        self.__check_len_sect__('Atoms', divide)
        self.atoms_sect, divide = self.__get_numeric_section__(divide, 9)

        # Get bonds section
        divide = divide.split("angles")
        self.__check_len_sect__('Bonds', divide)
        self.bonds_sect, divide = self.__get_numeric_section__(divide, 4)

        # Get angles and dihedral section
        divide = divide.split("dihedrals")
        self.__check_len_sect__('Angles', divide)
        self.angles_sect, self.dihedrals_sect = self.__get_numeric_section__(divide, 5)

        self.dihedrals_sect = self.dihedrals_sect.replace("dihedrals", "")

    def __parse_params_sect__(self):
        """
        Will parse the parameters section and store the values in metadata
        """
        for line in self.params_sect.split("\n"):
            if not line:   continue  # remove whitespace

            words = line.split()
            if len(words) == 2:
                self.metadata[words[1].strip()] = type_check.eval_type(words[0])

            elif len(words) == 3:
                self.metadata["_".join(words[1:]).strip()] = type_check.eval_type(words[0])

            elif len(words) == 4:
                self.metadata[words[2].strip()] = type_check.eval_type(words[0])
                self.metadata[words[3].strip()] = type_check.eval_type(words[1])

    def __parse_masses_sect__(self):
        """
        Will parse the masses and store values in metadata
        """
        masses = {}
        for line in self.masses_sect.split("\n"):
            if not line: continue

            words = line.split()
            if len(words) == 2:
                masses[str(words[0])] = float(words[1])

        self.metadata['masses'] = masses

    def __parse_numeric_section__(self, sect_txt, headers):
        """
        Will parse a numeric data section and store values in the data dict
        """
        fp = StringIO(sect_txt.strip('\n '))
        return pd.read_csv(fp, delim_whitespace=True, names=headers)


    def __check_len_sect__(self, sect, div):
        """
        A quick function to check for errors in the input file section declarations

        Inputs:
          * sect <str> => The section to report as missing
          * div <list<str>> => The split of the filetxt for the section above
        """
        err_msg = f"Error in the file {self.filepath}."
        err_msg += "\n\n"
        err_msg += f"Mulitple declarations of the word {sect}. "
        err_msg += f"Only 1 is allowed."

        if len(div) != 2:
            raise SystemExit(err_msg)

    def __get_numeric_section__(self, sects, num_cols):
        """
        Will parse a numeric section from txt based on how many columns of data it has

        Inputs:
           * sects <list<str>> => 2 possible sections
           * num_cols <int> => The number of columns the section is allowed
        Outputs:
           <str> The section txt and the non-section txt
        """
        check = [i for i in sects[0].split("\n") if i]
        inds = self.__search_in_list_of_str__(check,
                                              lambda s: len(s.split()) == num_cols)
        if len(inds) > 0:
           ind = 0
        else:
           check = [i for i in sects[1].split("\n") if i]
           inds = self.__search_in_list_of_str__(check,
                                                 lambda s: len(s.split()) == num_cols)
           if len(inds) > 0:
               ind = 1
           else:
               raise SystemExit(f"Error in parsing file {self.filepath}. Please check the numerical data.")

        return  sects[ind], sects[1 - ind]



    def __search_in_list_of_str__(self, str_list, fnc, max_inds=False):
        """
        Will apply a function to a list of strings to search for a given pattern.

        The index where the pattern occurs is returned.
        The pattern is specified by the fnc

        Inputs:
           * str_list <list<str>> => A list to check
           * fnc <function> => Function to apply to check for pattern
           * max_inds <int> OPTIONAL => A maximum number of indices to find

        Outputs:
           <list<int>> A list of all indices where the pattern appears
        """
        inds = []
        for i, string in enumerate(str_list):
            if fnc(string):
                inds.append(i)

        if max_inds is not False and len(inds) > max_inds:
           raise SystemExit(f"Found too many inds found that fit the pattern in {fnc}")

        return inds

    def _set_xyz_data_(self):
        """
        Will set the xyz_data variable and cols and timesteps for the write_xyz function to use later.

        This doesn't affect the data it is only used for writing.
        """
        with open("src/data/period_table.json", 'r') as f:
            pt = json.load(f)
        at_cvt = {int(pt[i]['atomic_weight']): pt[i]['abbreviation'] for i in pt
                  if pt[i]['atomic_weight'] is not None}

        # Save the xyz data to allow the write_xyz function to find it and write it
        self.xyz_data = [self.data['atoms'][['x', 'y', 'z']].to_numpy()]
        self.xyz_data = np.array(self.xyz_data)

        self.cols = self.data['atoms']['at_type']

        for i in self.cols.unique():
            mass = self.metadata['masses'][str(i)]
            mass = int(mass)
            self.cols[self.cols == i] = at_cvt[mass]
        self.cols = np.array([self.cols], dtype=str)

        self.timesteps = np.array([0.0] * len(self.cols))
