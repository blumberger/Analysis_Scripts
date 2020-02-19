#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to store objects that perform I/O tasks on lammps files.
"""

import pandas as pd
import numpy as np
from io import StringIO
from collections import Counter
import re
import json
import copy

from src.io_utils import general_io as gen_io
from src.system import type_checking as type_check
from src.data import consts


class Lammps_Log_File(gen_io.DataFileStorage):
    """
    Will read the lammps log file.

    A class that will loop over all lines and parse the CSVs from the
    lammps log file.
    """
    csv_data = []
    metadata = {'file_type': 'log_csv', 'number_atoms': 0,
                'total_run_time': 0, 'run_times': []}
    name = "Lammps Log File"
    _write_types = ('csv', 'txt',)

    def parse(self):
        """
        Will loop over all steps and parse the csv from the lammps log file.
        """
        self._ltxt = self.file_txt.split("\n")
        self.get_csv_lines(50)
        self.read_csv_lines()
        self.get_metadata()

    def get_csv_lines(self, same_line_tolerance=100):
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
        self._csv_starts, self._csv_ends = [], []

        # Find where there are csv lines
        for line_num, line in enumerate(self._ltxt):

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
                    self._csv_starts.append(line_num - 1 - same_line_tolerance)

                is_csv_line = True

            else:
                if is_csv_line is True:
                    self._csv_ends.append(line_num)

                is_csv_line = False

            split_num = new_split_num

        # self.csv_line_nums = np.arange(len(self._ltxt))

        # Some error checking
        if len(self._csv_starts) != len(self._csv_ends):
            if len(self._csv_starts) == len(self._csv_ends) - 1:
                self._csv_ends.append(len(ltxt))
            else:
                err_msg = f"Number of csv block starts = {len(self._csv_starts)}"
                err_msg += "\n"
                err_msg += f"Number of csv block ends = {len(self._csv_ends)}"
                err_msg += "\n"
                err_msg += f"Can't find the end of the csv block in the file {self.filepath}"
                raise EOFError(err_msg)

        # Get the lines of the csvs
        self._csv_lines = {i: '\n'.join(self._ltxt[start: end])
                          for i, (start, end) in enumerate(zip(self._csv_starts, self._csv_ends))}

    def read_csv_lines(self):
        """
        Will read the csv lines (according to 'get_csv_lines()') into DataFrames.

        The 2 lists: self._csv_starts and self._csv_ends tell the code where the DataFrame
        starts and ends.
        """
        for ifp, key in enumerate(self._csv_lines):
            # Create file-like object from a string
            fp = StringIO(self._csv_lines[key])
            self.csv_data.append(pd.read_csv(fp, delim_whitespace=True))

    def get_metadata(self):
        """
        Will get extra metadata from a log file such as number of atoms etc...

        Will use regex to seach for key phrases to get metadata.
        """
        # First get the non-csv text
        if len(self._csv_starts) == 0:
            f_txt = self.file_txt

        else:
            f_txt = '\n'.join(self._ltxt[:self._csv_starts[0]])
            for i in range(len(self._csv_starts) - 1):
                start = self._csv_starts[i+1]
                end = self._csv_ends[i]
                f_txt += '\n'.join(self._ltxt[end:start])
            f_txt += '\n'.join(self._ltxt[self._csv_ends[-1]:])

        # Search the non-csv text for the strings
        natom = re.findall("with [0-9]+ atoms", f_txt)
        if len(natom) > 0:
                self.metadata['number_atoms'] = int(natom[0][5:-5])

        run_times = re.findall("run [0-9]+", f_txt)
        if len(run_times):
                self.metadata['run_times'] = [int(i[4:]) for i in run_times]
                self.metadata['total_run_time'] = sum(self.metadata['run_times'])

    def __repr__(self):
        return "Lammps Log File Class"


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
    Will read the lammps input data file (lammps.in).

    A class that will loop over all lines and parse the CSVs from the
    lammps log file.
    """
    _write_types = ('csv', 'xyz',)
    csv_data = {}
    metadata = {'file_type': 'log_csv'}
    sects = ('masses', 'atoms', 'bonds', 'angles', 'dihedrals',)
    name = "Lammps Input Data"

    def parse(self):
        """
        Will loop over all steps and parse the csv from the lammps log file.
        """
        self._ltxt = [i for i in self.file_txt.split("\n") if i]

        ftxt = self.file_txt.lower()

        # Should check the order of the sections here!

        # Get the sections in the file
        self.divide_sections()
        self.parse_params_sect()
        self.parse_masses_sect()

        headers = ("id", "mol_id", "at_type", "x", "y", "z", "ix", "iy", "iz",)
        self.csv_data['atoms'] = self.parse_numeric_section(self.atoms_sect, headers)

        headers = ("id", "at_type", "at1", "at2")
        self.csv_data['bonds'] = self.parse_numeric_section(self.bonds_sect, headers)

        headers = ("id", "at_type", "at1", "at2", "at3")
        self.csv_data['angles'] = self.parse_numeric_section(self.angles_sect, headers)

        headers = ("id", "at_type", "at1", "at2", "at3", "at4")
        self.csv_data['dihedrals'] = self.parse_numeric_section(self.dihedrals_sect, headers)

    def divide_sections(self):
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
        ftxt = "\n".join(self._ltxt).lower()

        # Get first section
        divide = ftxt.lower().split('masses')
        self.check_len_sect('masses', divide)    # error checking
        ind = self.search_in_list_of_str(divide, lambda s: "xlo xhi" in s, 1)[0]
        self.params_sect = divide[ind]

        # Get masses section
        divide = divide[1-ind]
        divide = divide.split('atoms')
        self.check_len_sect('Masses', divide)     # error checking
        # Find the section with 2 columns
        self.masses_sect, divide = self.get_numeric_section(divide, 2)

        # Get atoms section
        divide = divide.split("bonds")
        self.check_len_sect('Atoms', divide)
        self.atoms_sect, divide = self.get_numeric_section(divide, 9)

        # Get bonds section
        divide = divide.split("angles")
        self.check_len_sect('Bonds', divide)
        self.bonds_sect, divide = self.get_numeric_section(divide, 4)

        # Get angles and dihedral section
        divide = divide.split("dihedrals")
        self.check_len_sect('Angles', divide)
        self.angles_sect, self.dihedrals_sect = self.get_numeric_section(divide, 5)

        self.dihedrals_sect = self.dihedrals_sect.replace("dihedrals", "")

    def parse_params_sect(self):
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

    def parse_masses_sect(self):
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

    def parse_numeric_section(self, sect_txt, headers):
        """
        Will parse a numeric data section and store values in the data dict
        """
        fp = StringIO(sect_txt.strip('\n '))
        return pd.read_csv(fp, delim_whitespace=True, names=headers)


    def check_len_sect(self, sect, div):
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

    def get_numeric_section(self, sects, num_cols):
        """
        Will parse a numeric section from txt based on how many columns of data it has

        Inputs:
           * sects <list<str>> => 2 possible sections
           * num_cols <int> => The number of columns the section is allowed
        Outputs:
           <str> The section txt and the non-section txt
        """
        check = [i for i in sects[0].split("\n") if i]
        inds = self.search_in_list_of_str(check,
                                              lambda s: len(s.split()) == num_cols)
        if len(inds) > 0:
           ind = 0
        else:
           check = [i for i in sects[1].split("\n") if i]
           inds = self.search_in_list_of_str(check,
                                                 lambda s: len(s.split()) == num_cols)
           if len(inds) > 0:
               ind = 1
           else:
               raise SystemExit(f"Error in parsing file {self.filepath}. Please check the numerical data.")

        return  sects[ind], sects[1 - ind]



    def search_in_list_of_str(self, str_list, fnc, max_inds=False):
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

    def set_xyz_data(self):
        """
        Will set the xyz_data variable and cols and timesteps for the write_xyz function to use later.

        This doesn't affect the data it is only used for writing.
        """
        with open(consts.PT_FILEPATH, 'r') as f:
            pt = json.load(f)
        at_cvt = {int(pt[i]['atomic_weight']): pt[i]['abbreviation'] for i in pt
                  if pt[i]['atomic_weight'] is not None}

        # Save the xyz data to allow the write_xyz function to find it and write it
        self.xyz_data = [self.csv_data['atoms'][['x', 'y', 'z']].to_numpy()]
        self.xyz_data = np.array(self.xyz_data)

        # Get the cols
        self.cols = self.csv_data['atoms']['at_type']
        for i in self.cols.unique():
            mass = self.metadata['masses'][str(i)]
            mass = int(mass)
            self.cols[self.cols == i] = at_cvt[mass]
        self.cols = np.array([self.cols], dtype=str)

        # Get the timesteps
        self.timesteps = np.array([0.0] * len(self.cols))

class Lammps_Dump(gen_io.DataFileStorage):
    """
    Will read a Lammps snapshot (dump) file and store the data.

    This inherits from gen_io.DataFileStorage see this for more info.

    Inputs:
        * filepath <str> => The path to the file to be loaded.
    """
    _write_types = ('csv', 'xyz',)
    name = "Lammps Dump"
    _defaults = {'coordinate_wrapping': 'unwrapped'}
    
    def parse(self):
        """
        Will parse the file text and store the data.
        """
        # Split the file into sections
        self._all_items = [s.strip().split("\n")
                           for s in self.file_txt.split("ITEM: ") if s]

        # Get data names
        self._data_names = self.get_metadata_titles()

        # Get any extra headers for CSVs
        self._headers = [re.findall("[a-z]+", item[0])
                         for item in self._all_items]

        # Parse each item
        for item_num, item in enumerate(self._all_items):
            name = self._data_names[item_num]

            # Deal with simple key-value pairs
            if not self._headers[item_num] and len(item) == 2:
                self.metadata[name] = type_check.eval_type(item[1])

            # Parse the box bounds keyword
            elif self._data_names[item_num] == 'box_bounds':
                self.parse_box_bounds(item_num)

            # Parse the final dump
            else:
                self.parse_dump(item_num)

        # Fix periodic BCs
        self.wrapped_csv = copy.deepcopy(self.csv_data)
        self.fix_wrapping()

    def parse_dump(self, item_num):
        """
        Will parse the dump CSV data.

        Inputs:
            * item_num <int> => The index of the item in the self._all_items list.
        """
        item = self._all_items[item_num]
        fp = StringIO('\n'.join(item[1:]))
        self.csv_data = pd.read_csv(fp, names=self._headers[item_num],
                                    delim_whitespace=True)

    def parse_box_bounds(self, item_num):
        """
        Will parse the Box Bounds item from the snapshot file.

        Inputs:
            * item_num <int> => The index of the item in the self._all_items list.
        """
        item = self._all_items[item_num]
        if all(j in self._headers[item_num] for j in ('xy', 'xz', 'yz',)):
            xlo, xhi, xy = (float(i) for i in item[1].split())
            ylo, yhi, xz = (float(i) for i in item[2].split())
            zlo, zhi, yz = (float(i) for i in item[3].split())
            for name, val in (('xlo', xlo), ('xhi', xhi), ('xy', xy),
                              ('ylo', ylo), ('yhi', yhi), ('xz', xz),
                              ('zlo', zlo), ('zhi', zhi), ('yz', yz)):
                self.metadata[name] = val
            self.metadata['cell_type'] = "triclinic"
        else:
            xlo, xhi = (float(i) for i in item[1].split())
            ylo, yhi = (float(i) for i in item[2].split())
            zlo, zhi = (float(i) for i in item[3].split())
            for name, val in (('xlo', xlo), ('xhi', xhi),
                              ('ylo', ylo), ('yhi', yhi),
                              ('zlo', zlo), ('zhi', zhi)):
                self.metadata[name] = val
            self.metadata['cell_type'] = "orthogonal"
        self.get_unit_vecs()

    def get_unit_vecs(self):
        """
        Will get the unit vectors: a, b, c for the crystal system.
        """
        xlo, xhi = self.metadata['xlo'], self.metadata['xhi']
        ylo, yhi = self.metadata['ylo'], self.metadata['yhi']
        zlo, zhi = self.metadata['zlo'], self.metadata['zhi']

        if self.metadata['cell_type'] == 'triclinic':
            xy, xz, yz = self.metadata['xy'], self.metadata['xz'], self.metadata['yz']
            xhi = xhi - (xy + xz)

            min_y = yz if yz < 0 else 0
            max_y = yz if yz > 0 else 0
            ylo = ylo - min_y
            yhi = yhi - max_y

            self.metadata['a'] = np.array([xhi - xlo, 0, 0])
            self.metadata['b'] = np.array([xy, yhi - ylo, 0])
            self.metadata['c'] = np.array([xz, yz, zhi - zlo])
        else:
            self.metadata['a'] = np.array([xhi - xlo, 0, 0])
            self.metadata['b'] = np.array([0, yhi - ylo, 0])
            self.metadata['c'] = np.array([0, 0, zhi - zlo])

    def get_metadata_titles(self):
        """
        Will get the title of the metadata within the file.
        """
        words_to_remove = ('of',)
        item_names = (s[0] for s in self._all_items if s)
        item_names = (re.findall("[A-Z]+", i) for i in item_names)
        item_names = ([s.lower().strip() for s in i]
                       for i in item_names)
        item_names = ([s for s in i if s not in words_to_remove]
                      for i in item_names)
        return ['_'.join(i) for i in item_names]

    def fix_wrapping(self):
        """
        Will translate atom coords to fix wrapping of coords in periodic systems
        """
        # Check we can do the wrapping
        if not all(j in self.csv_data.columns for j in ('ix', 'iy', 'iz', 'x', 'y', 'z',)):
            return self.csv_data

        # Apply the wrapping
        unit_vectors = self.metadata['a'], self.metadata['b'], self.metadata['c']
        for unit_vec, wrap_dim in zip(unit_vectors, ('ix', 'iy', 'iz')):
            for idim, dim in enumerate('xyz'):
                if unit_vec[idim] != 0:
                    self.csv_data[dim] += self.csv_data[wrap_dim] * unit_vec[idim]

    def set_xyz_data(self):
        """
        Will set the xyz_data variable and cols and timesteps for the write_xyz function to use later.

        This doesn't affect the data it is only used for writing.
        """
        if self.metadata['coordinate_wrapping'] == 'unwrapped':
            df = self.csv_data
        else:
            df = self.wrapped_csv

        # Set the xyz data
        self.xyz_data = np.array([df[['x', 'y', 'z']].to_numpy()])

        # Set the xyz column data
        self.cols = np.array([df['type']]).astype(str)
        self.cols[0][self.cols[0] == '1'] = 'C'
        self.cols[0][self.cols[0] == '2'] = 'H'

        # Set the timesteps
        self.timesteps = np.array([0])
        if 'timestep' in self.metadata:
            self.timesteps = np.array([self.metadata['timestep']])
