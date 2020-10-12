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
import matplotlib.pyplot as plt

from src.io_utils import general_io as gen_io
from src.system import type_checking as type_check
from src.data import consts
from src.input_file import input_file_types as inp_types
from src.calc import molecule_utils as mol_utils
from src.parsing import general_parsing as gen_parse
from src.plot import general_plot as gen_plot


class Lammps_Log_File(gen_io.DataFileStorage):
    """
    Will read the lammps log file.

    A class that will loop over all lines and parse the CSVs from the
    lammps log file.

    Inputs:
        * filepath <str> => The path to the file to be loaded.
    """
    csv_data = []
    metadata = {'file_type': 'csv', 'number_atoms': 0,
                'total_run_time': 0, 'run_times': []}
    name = "Lammps Log File"
    _write_types = ('csv', 'txt',)

    def __init__(self, filepath):
        super().__init__(filepath)

    def _parse_(self):
        """
        Will loop over all steps and parse the csv from the lammps log file.
        """
        self._ltxt = self.file_txt.split("\n")
        self._set_csv_lines(50)
        self._read_csv_lines()
        self._get_metadata()
        self.append_csvs()

    def append_csvs(self):
        """
        Will append the csv files into multiple csvs with the same columns
        """
        # Get all unique headers
        col_heads = []
        for df in self.csv_data:
            heads = '|'.join(df.columns)
            if heads not in col_heads:
               col_heads.append(heads)

        # Collect similar dataframes
        self.collected_csv_data = [pd.DataFrame() for i in range(len(col_heads))]
        for df in self.csv_data:
            heads = '|'.join(df.columns)
            index = col_heads.index(heads)
            self.collected_csv_data[index] = self.collected_csv_data[index].append(df)

    def _set_csv_lines(self, same_line_tolerance=20):
        """
        Will determine which lines contain data in a csv format.

        This will loop over all lines deciding which ones are csv lines and find the start and the end
        of csv blocks.

        Will set the following attributes:
            * _csv_starts <list> => The start line of the csv data
            * _csv_ends <list> => The end line of the csv data
            * _csv_lines <dict> => The actual csv lines

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
            if len(self._csv_starts) == len(self._csv_ends) + 1:
                self._csv_ends.append(len(self._ltxt))

            else:
                err_msg = f"Number of csv block starts = {len(self._csv_starts)}"
                err_msg += "\n"
                err_msg += f"Number of csv block ends = {len(self._csv_ends)}" + "\n"
                err_msg += f"    * csv_starts = {self._csv_starts}" + "\n"
                err_msg += f"    * csv_ends = {self._csv_ends}" + "\n"
                err_msg += "\n"
                err_msg += f"Can't find the end of the csv block in the file {self.filepath}"
                raise EOFError(err_msg)

        # Get the lines of the csvs
        self._csv_lines = {i: '\n'.join(self._ltxt[start: end])
                          for i, (start, end) in enumerate(zip(self._csv_starts, self._csv_ends))}
        return self._csv_lines

    def _read_csv_lines(self):
        """
        Will read the csv lines (according to '_set_csv_lines()') into DataFrames.

        The 2 lists: self._csv_starts and self._csv_ends tell the code where the DataFrame
        starts and ends.
        """
        for ifp, key in enumerate(self._csv_lines):
            # Create file-like object from a string
            fp = StringIO(self._csv_lines[key])
            self.csv_data.append(pd.read_csv(fp, delim_whitespace=True))

    def _get_metadata(self):
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


class Lammps_Input_File(Lammps_Log_File):
    """
    Will read a input lammps data file.

    This is the file that defines the geometry of the system for lammps (e.g. where the atoms are
    which atom is bonded to which, the masses and types of atoms etc..)

    As with the other file data structures in this codebase the smaller data like ints, floats,
    small lists are stored in the self.metadata dict. The larger data is stored in separate 
    attributes. The important attributes in this class are:
        * Atoms <pd.DataFrame> => The atomic coords
        * Bonds <pd.DataFrame> => The bonding of the system
        * Angles <pd.DataFrame> => The angles
        * Dihedrals <pd.DataFrame> => The dihedrals

    This class inherits from the Lammps_Log_File for the get_all_csv_lines function. It
    probably isn't an ideal way to set things up as the lammps input file isn't really
    a logical child of the log file but it works...

    Inputs:
        * filepath <str> => The path to the file to be loaded.
    """
    csv_data = []
    metadata = {'file_type': 'lammps_input', 'number_atoms': 0}
    name = "Lammps Input File"
    _write_types = ('csv', 'txt',)
    _title_names = {
                    'atoms': (('index', 'mol_num', 'at_type', 'q', 'x', 'y', 'z', 'ix', 'iy', 'iz'),
                              ('index', 'mol_num', 'at_type', 'x', 'y', 'z', 'ix', 'iy', 'iz')),
                    'bonds': (('index', 'bond_type', 'at_1', 'at_2'),),
                    'angles': (('index', 'angle_type', 'at_1', 'at_2', 'at_3'),),
                    'dihedrals': (('index', 'dihedral_type', 'at_1', 'at_2', 'at_3', 'at_4'),),
                    'velocities': (('index', 'vx', 'vy', 'vz'),),
                    'pair coeffs': (('index', 'coeffA', 'coeffB'),),
                    'bond coeffs': (('index', 'coeffA', 'coeffB'),),
                    'angle coeffs': (('index', 'coeffA', 'coeffB'),),
                    'dihedral coeffs': (('index', 'coeffA', 'coeffB', 'coeffC', 'coeffD'),),
                    'masses': (('index', 'mass'),)
                   }

    def __init__(self, filepath):
        # Use the grandparent's init
        gen_io.DataFileStorage.__init__(self, filepath)

    def _parse_(self):
        """
        The main parsing function, called from the gen_io.DataFileStorage __init__.
        """
        self._ltxt = self.file_txt.split("\n")

        # Get the places that look like CSV lines
        self._set_csv_lines(same_line_tolerance=20)
        # Set the CSV names
        swap_keys = {}
        for i in self._csv_lines:
            key, _ = gen_parse.rm_comment_from_line(self._ltxt[self._csv_starts[i]-2])
            key = key.strip().lower()
            swap_keys[key] = i
        for i in swap_keys:
            self._csv_lines[i] = self._csv_lines[swap_keys[i]]
            self._csv_lines.pop(swap_keys[i])

        # Add any sections not caught by the set_csv lines function
        self._check_other_sects_()

        self.csv_data = self._read_csv_lines()

        self._set_metadata()

    def _check_other_sects_(self):
        """
        Will add any sections that haven't been added by the the _set_csv_lines() function.
        """
        txt = self.file_txt.lower()
        for sect_name in self._title_names:
            if sect_name == 'masses': continue
            if sect_name in self._csv_lines or sect_name not in txt:
                continue

            self._get_extra_lines_(sect_name)
            # print("\n\n" + f"Haven't parsed the '{sect_name.title()}' section.")
            # print(f"Please change the function Lammps_Input_File._check_other_sects_" +
            #       f" in order to implement this.")

    def _get_extra_lines_(self, sect_name):
        for i, line in enumerate(self._ltxt):
            if sect_name in line.lower():
                start_line = i

        num_entries_b = len(self._ltxt[start_line+2].split())
        for i, line in enumerate(self._ltxt[start_line+3:]):
            words = line.split()
            if num_entries_b == len(words):
                continue
            end_line = i + start_line + 3
            break

        self._csv_lines[sect_name] = '\n'.join(self._ltxt[start_line+2:end_line])

    def _read_csv_lines(self):
        """
        Will read the csv lines (according to '_set_csv_lines()') into DataFrames.

        The dict _csv_lines passes the lines for each csv file.

        The list: _csv_titles tells the code what the csv files are called.
        """
        csv_data = {}
        for ifp, key in enumerate(self._csv_lines):
            # Create file-like object from a string
            fp = StringIO(self._csv_lines[key])
            if key not in self._title_names:
                num_cols = len(self._ltxt[self._csv_starts[ifp]].split())
                raise SystemError("\n\n" + f"Can't read the {key} section. "
                             + "\nPlease input the titles for the columns in the 'Lammps_Input_File'"
                             + " class.\n\nYou can do this by creating a new section in the '_title_names'"
                             + f" attribute called '{key}' e.g:" + "\n\t" 
                             + "\tself._title_names = {\n\t                             ...\n"
                             + f"\t                             '{key}': (" + "..., " * num_cols
                             + ")\n\t                           })")

            headers = self._get_csv_headers_(self._csv_lines[key], key)
            csv_data[key] = pd.read_csv(fp, delim_whitespace=True, names=headers)
        return csv_data

    def _get_csv_headers_(self, lines, key):
        """
        Will find the correct headers for the section in the csv lines.

        Inputs:
            * lines <str> => The csv lines
            * key <str>   => The title fo the section

        Outputs:
            headers <list>
        """
        first_line = lines[:500].splitlines()[0]
        words = first_line.split()
        num_words = len(words)
        
        poss_vals = []
        for i, vals in enumerate(self._title_names[key]):
            if len(vals) == num_words:
                poss_vals.append(i)


        if len(poss_vals) > 1:
            raise SystemExit(f"Can't determine the csv headers for the section: {key}" +
                             "\n\nPlease check your file, change the self._title_names"+
                             " variable in src/io_utils/lammps.py or add a new test "  +
                             "within the function Lammps_Input_File._get_csv_headers_().")

        elif len(poss_vals) == 0:
            raise SystemExit(f"No valid csv headers for the section: {key}" +
                             "\n\nPlease check your file or change the self._title_names"+
                             " variable in src/io_utils/lammps.py.")

        else:
            return self._title_names[key][poss_vals[0]]


    def _set_metadata(self):
        """
        Will loop over lines and get the metadata from the header of the file.
        """
        header_ltxt = self._ltxt[:self._csv_starts[0]-2]
        self.metadata['header_txt'] = '\n'.join(header_ltxt)

        self.metadata['title_line'] = self._ltxt[0]

        # First parse the masses and remove them
        mass_start, del_lines = 0, []
        for line_num, line in enumerate(header_ltxt):
            if line == "Masses":
                for j, line2 in enumerate(header_ltxt[line_num:]):
                    del_lines.append(j + line_num)
                    if len(line2.split()) == 2:
                        mass_start = line_num + j
                        break
                else: raise SystemError("Error parsing the 'Masses' section.")
                break

        self.metadata['atom_types'] = {}
        for i, line in enumerate(header_ltxt[mass_start:]):
            splitter = line.split()
            if len(splitter) == 2:
                atom_type = type_check.eval_type(splitter[0])
                mass = type_check.eval_type(splitter[1])
                self.metadata['atom_types'].setdefault(atom_type, {})['mass'] = mass
                at_type = mol_utils.get_atom_type(mass)
                self.metadata['atom_types'][atom_type]['details'] = at_type
                del_lines.append(i + mass_start)
            else:  break

        # Tidy up
        for i in reversed(sorted(del_lines)): del header_ltxt[i]


        # Parse the rest of the header file
        for line in header_ltxt:
            splitter = line.split()
            if len(splitter) % 2 == 0:
                for i in range(len(splitter) // 2):
                    # name = type_check.eval_type(splitter[1])
                    name = splitter[len(splitter) // 2 + i].strip()
                    val = type_check.eval_type(splitter[i])
                    self.metadata[name] = val

            if len(splitter) == 3:
                if splitter[2] == "types":
                    self.metadata[f"number {splitter[1]} types"] = type_check.eval_type(splitter[0])


class Write_Lammps_Input(gen_io.Write_File):
    """
    Will write a lammps input file.

    The main method in this class is creating a string that can be written to a file. This
    is done in the create_file_str method. This string is then handled by the parent class
    and written to a file.

    See gen_io.Write_File for more info.
    
    Inputs:
       * Data_Class <class> => The class containing all the data to be written
       * filepath <str>     => The path to the file to be written.
       * extension <str> OPTIONAL   => The file extension. Default is False.
    """
    ordered_sects = ('pair coeffs', 'bond coeffs', 'angle coeffs', 'dihedral coeffs', 
                     'atoms', 'bonds', 'angles', 'dihedrals', 'velocities', )

    def __set_data__(self):
        """
        Will check we have the correct info for writing the file.
        """
        if type(self.Data) == inp_types.Vars:
            if 'lammps_input' in self.Data:
                self.Data = self.Data['lammps_input']
            else:
                raise SystemError(f"Can't write data of types: {self.Data.keys()} as a Lammps input_file.")
        
        elif type(self.Data) == Lammps_Input_File:
            pass

        else:
            raise SystemError(f"Can't write data of types: {type(self.Data)} as a Lammps input_file.")

    def _create_header_str_(self):
        """Will create the header string for the input file."""
        s = self.Data.metadata['title_line'] + "\n\n"

        s += f"{self.Data.metadata['atoms']} atoms" + "\n"
        s += f"{self.Data.metadata['bonds']} bonds" + "\n"
        s += f"{self.Data.metadata['angles']} angles" + "\n"
        s += f"{self.Data.metadata['dihedrals']} dihedrals" + "\n\n"

        s += f"{self.Data.metadata['number atom types']} atom types" + "\n"
        s += f"{self.Data.metadata['number bond types']} bond types" + "\n"
        s += f"{self.Data.metadata['number angle types']} angle types" + "\n"
        s += f"{self.Data.metadata['number dihedral types']} dihedral types" + "\n\n"

        s += f"{self.Data.metadata['xlo']} {self.Data.metadata['xhi']} xlo xhi" + "\n"
        s += f"{self.Data.metadata['ylo']} {self.Data.metadata['yhi']} ylo yhi" + "\n"
        s += f"{self.Data.metadata['zlo']} {self.Data.metadata['zhi']} zlo zhi" + "\n"

        if 'xy' in self.Data.metadata and 'xz' in self.Data.metadata and 'yz' in self.Data.metadata:
              s += f"{self.Data.metadata['xy']} {self.Data.metadata['xz']} {self.Data.metadata['yz']} xy xz yz" + "\n\n"
        else: s += "\n"

        s += f"Masses\n"
        for i in self.Data.metadata['atom_types']:
            s += "\n" + f"{i} {self.Data.metadata['atom_types'][i]['mass']}"

        return s

    def _get_csv_strs_(self):
        """Will return the csv data as a string."""
        s = ""
        for key in self.ordered_sects:
            if key not in self.Data.csv_data: continue
            s += "\n\n" + key.title() + "\n\n"
            s += self.Data.csv_data[key].to_string(index=False, header=False, justify="left")

        return s

    def create_file_str(self):
        """
        Will create the text that can be written to a file.
        """
        self.__set_data__()
        s = self._create_header_str_()
        s += self._get_csv_strs_()

        return s


class Lammps_Data_File(gen_io.DataFileStorage):
    """
    Will read the lammps input data file (lammps.in).

    """
    _write_types = ('csv', 'xyz',)
    csv_data = {}
    metadata = {'file_type': 'log_csv'}
    sects = ('masses', 'atoms', 'bonds', 'angles', 'dihedrals',)
    name = "Lammps Input Data"

    def __init__(self, filepath):
        super().__init__(filepath)
        self.metadata = copy.deepcopy(self.metadata)
        self.csv_data = copy.deepcopy(self.csv_data)

    def _parse_(self):
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

        headers = { 9: ("id", "mol_id", "at_type", "x", "y", "z", "ix", "iy", "iz",),
                   10: ("id", "mol_id", "at_type", "q", "x", "y", "z", "ix", "iy", "iz",),
                  }

        ncol = len([i for i in self.atoms_sect[0:1000].split("\n") if i][0].split())
        if ncol not in headers:
            raise SystemExit(f"I don't have some header inputs for an atom section with {ncol} columns." +
                              "\nI" + f" can handle ({', '.join(map(str, headers.keys()))}) numbers " +
                              "of columns.")
        at_headers = headers[ncol]
        self.csv_data['atoms'] = self.parse_numeric_section(self.atoms_sect, at_headers)

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
        self.masses_sect, divide = self.get_numeric_section(divide, 2, 'masses')

        # Get atoms section
        divide = divide.split("bonds")
        self.check_len_sect('Atoms', divide)
        try:
            self.atoms_sect, divide = self.get_numeric_section(divide, 9, 'atoms')
        except SystemExit as e:
            self.atoms_sect, divide = self.get_numeric_section(divide, 10, "atoms")

        # Get bonds section
        divide = divide.split("angles")
        self.check_len_sect('Bonds', divide)
        self.bonds_sect, divide = self.get_numeric_section(divide, 4, 'bonds')

        # Get angles and dihedral section
        divide = divide.split("dihedrals")
        self.check_len_sect('Angles', divide)
        self.angles_sect, self.dihedrals_sect = self.get_numeric_section(divide, 5, 'angles')

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
        at_types = {int(i): {'mass': i, 
                             'details': mol_utils.get_atom_type(masses[i])} for i in masses
                   }
        self.metadata['atom_types'] = at_types


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

    def get_numeric_section(self, sects, num_cols, sect_name):
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
                inds = self.search_in_list_of_str(check,
                                                 lambda s: len(s.split("\t")) == num_cols)
                if len(inds) == 0:
                    msg = "Error in parsing a lammps input data file."
                    msg += f"The {sect_name.title()} section cannot "
                    msg += f"be parsed into {num_cols} columns." + "\n\n"
                    raise SystemExit(msg+f"Parsing Error {self.filepath}.")

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

    # def set_xyz_data(self):
    #     """
    #     Will set the xyz_data variable and cols and timesteps for the write_xyz function to use later.

    #     This doesn't affect the data it is only used for writing.
    #     """
    #     with open(consts.PT_FILEPATH, 'r') as f:
    #         pt = json.load(f)
    #     at_cvt = {int(pt[i]['atomic_weight']): pt[i]['abbreviation'] for i in pt
    #               if pt[i]['atomic_weight'] is not None}

    #     # Save the xyz data to allow the write_xyz function to find it and write it
    #     self.xyz_data = [self.csv_data['atoms'][['x', 'y', 'z']].to_numpy()]
    #     self.xyz_data = np.array(self.xyz_data)

    #     # Get the cols
    #     self.cols = self.csv_data['atoms']['at_type']
    #     for i in self.cols.unique():
    #         mass = self.metadata['masses'][str(i)]
    #         mass = int(mass)
    #         self.cols[self.cols == i] = at_cvt[mass]
    #     self.cols = np.array([self.cols], dtype=str)

    #     # Get the timesteps
    #     self.timesteps = np.array([0.0] * len(self.cols))

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
    def _parse_(self):
        """
        Will parse the file text and store the data.
        """
        # Split the file into sections
        self._all_items = [s.strip().split("\n")
                           for s in self.file_txt.split("ITEM: ") if s]

        # Get data names
        self._data_names = self._get_metadata_titles()

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

        # Add some useful metadata
        unique_mols = self.csv_data['mol'].unique()
        self.metadata['nmol'] = len(unique_mols)
        self.metadata['natom'] = len(self.csv_data)
        if 'atoms_per_molecule' not in self.metadata:
            self.metadata['atoms_per_molecule'] = np.sum(self.csv_data['mol'] == unique_mols[0])

        # Fix periodic BCs
        self.csv_data['timestep'] = self.metadata['timestep']
        self.wrapped_csv = copy.deepcopy(self.csv_data)
        self.unwrap_coords()

    def set_data(self):
        """
        Will set the wrapped_csv or unwrapped_csv to be the csv_data.

        This depends on the metadata inputted.
        """
        if self.metadata['coordinate_wrapping'] == 'wrapped' or not self.unwrapped_avail:
            self.csv_data = self.wrapped_csv
        else:
            self.csv_data = self.unwrapped_csv

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
            
            for name, val in (('xlo_bound', xlo), ('xhi_bound', xhi), ('xy', xy),
                              ('ylo_bound', ylo), ('yhi_bound', yhi), ('xz', xz),
                              ('zlo_bound', zlo), ('zhi_bound', zhi), ('yz', yz)):
                self.metadata[name] = val

            min_x = min([0, xy, xz, xy+xz])
            max_x = max([0, xy, xz, xy+xz])
            self.metadata['xlo'] = xlo - min_x
            self.metadata['xhi'] = xhi - max_x

            min_y, max_y = min((yz, 0)), max((yz, 0))
            self.metadata['ylo'] = ylo - min_y
            self.metadata['yhi'] = yhi - max_y

            self.metadata['zlo'], self.metadata['zhi'] = zlo, zhi

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

            self.metadata['a'] = np.array([xhi - xlo, 0.0, 0.0])
            self.metadata['b'] = np.array([xy, yhi - ylo, 0.0])
            self.metadata['c'] = np.array([xz, yz, zhi - zlo])

        else:
            self.metadata['a'] = np.array([xhi - xlo, 0, 0])
            self.metadata['b'] = np.array([0, yhi - ylo, 0])
            self.metadata['c'] = np.array([0, 0, zhi - zlo])

    def _get_metadata_titles(self):
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

    def proj(self, vec1, vec2):
        return np.dot(vec1, vec2) / np.linalg.norm(vec1)

    def unwrap_coords(self):
        """
        Will translate atom coords to fix wrapping of coords in periodic systems
        """
        self.unwrapped_avail = True
        self.unwrapped_csv = copy.deepcopy(self.csv_data)

        # Check we can do the wrapping
        if not all(j in self.unwrapped_csv.columns for j in ('ix', 'iy', 'iz', 'x', 'y', 'z',)):
            self.unwrapped_avail = False
            return self.unwrapped_csv

        unit_vectors = self.metadata['a'], self.metadata['b'], self.metadata['c']
        for unit_vec, wrap_dim in zip(unit_vectors, ('ix', 'iy', 'iz')):
            for idim, dim in enumerate('xyz'):
                if unit_vec[idim] != 0:
                    self.unwrapped_csv[dim] += self.unwrapped_csv[wrap_dim] * unit_vec[idim]

    def remove_split_mols(self):
        """
        Will remove the split molecules from the system (keeping only the whole mols in the center).
        """
        self.unwrapped_avail = True
        self.rm_split_mols = copy.deepcopy(self.csv_data)

        unit_vectors = self.metadata['a'], self.metadata['b'], self.metadata['c']

        # Check we can do the wrapping
        if not all(j in self.rm_split_mols.columns for j in ('ix', 'iy', 'iz', 'x', 'y', 'z', 'mol')):
            self.unwrapped_avail = False
            return self.rm_split_mols

        # Find the split mols
        mol_crds = self.rm_split_mols[['mol', 'ix', 'iy', 'iz']].groupby("mol", axis=0).std()
        mask = (mol_crds['ix'] != 0) | (mol_crds['iy'] != 0) | (mol_crds['iz'] != 0)
        split_mols = mol_crds[mask].index
        self.rm_split_mols['split'] = self.rm_split_mols['mol'].apply(lambda x: x in split_mols)

        # Get a map which maps the molecules that are still there to the original molecules
        split_ind = 0
        self.split_to_mol = {}
        ats_per_mol = self.metadata['atoms_per_molecule']
        nmol = len(self.rm_split_mols) // ats_per_mol
        for mol_ind in range(nmol):
            all_splits = self.rm_split_mols['split'][mol_ind*ats_per_mol : (mol_ind+1)*ats_per_mol]
            split = not any(all_splits)
            if split:
                self.split_to_mol[split_ind] = mol_ind
                split_ind += 1

        # Remove the split molecules together.
        self.rm_split_mols = self.rm_split_mols[self.rm_split_mols['split'] == False]


    def unwrap_add_images(self):
        """
        Will unwrap the molecules but for every molecule that has been split an image molecule
        will be put in every place the atoms of that molecule are.
        """
        self.add_img_csv = copy.deepcopy(self.csv_data)
        unit_vectors = self.metadata['a'], self.metadata['b'], self.metadata['c']

        # Check we can do the wrapping
        if not all(j in self.add_img_csv.columns for j in ('ix', 'iy', 'iz', 'x', 'y', 'z', 'mol')):
            self.unwrapped_avail = False
            return self.add_img_csv

        # Find the split mols
        mol_crds = self.add_img_csv[['mol', 'ix', 'iy', 'iz']].groupby("mol", axis=0).std()
        mask = (mol_crds['ix'] != 0) | (mol_crds['iy'] != 0) | (mol_crds['iz'] != 0)
        split_mol_inds = mol_crds[mask].index
        self.add_img_csv.loc[:, 'split'] = self.add_img_csv.loc[:, 'mol'].apply(lambda x: x in split_mol_inds)

        for i in self.add_img_csv['mol'].unique():
            mol_df = self.add_img_csv[self.add_img_csv['mol'] == i]
            print(mol_df)

            for unit_vec, wrap_dim in zip(unit_vectors, ('ix', 'iy', 'iz')):
                for idim, dim in enumerate('xyz'):
                    if unit_vec[idim] != 0:
                        mol_df[dim] += mol_df[wrap_dim] * unit_vec[idim]

            break
        raise SystemExit


    def unwrap_split_mols(self):
        """
        Will wrap the molecules back as close to the unit cell as they go.
        """
        self.unwrapped_avail = True
        self.unwrapped_csv = copy.deepcopy(self.csv_data)

        unit_vectors = self.metadata['a'], self.metadata['b'], self.metadata['c']

        # Check we can do the wrapping
        if not all(j in self.unwrapped_csv.columns for j in ('ix', 'iy', 'iz', 'x', 'y', 'z', 'mol')):
            self.unwrapped_avail = False
            return self.unwrapped_csv

        # Find the split mols
        mol_crds = self.unwrapped_csv[['mol', 'ix', 'iy', 'iz']].groupby("mol", axis=0).std()
        mask = (mol_crds['ix'] != 0) | (mol_crds['iy'] != 0) | (mol_crds['iz'] != 0)
        split_mol_inds = mol_crds[mask].index
        self.unwrapped_csv.loc[:, 'split'] = self.unwrapped_csv.loc[:, 'mol'].apply(lambda x: x in split_mol_inds)

        # Wrap the split molecules together.
        unsplit = self.unwrapped_csv[self.unwrapped_csv['split'] == False]

        mask = self.unwrapped_csv['split'] == True
        tmp = self.unwrapped_csv.loc[mask, :]
        unit_vectors = self.metadata['a'], self.metadata['b'], self.metadata['c']
        for unit_vec, wrap_dim in zip(unit_vectors, ('ix', 'iy', 'iz')):
            for idim, dim in enumerate('xyz'):
                if unit_vec[idim] != 0:
                    self.unwrapped_csv.loc[mask, dim] += tmp.loc[:, wrap_dim] * unit_vec[idim]

    def append(self, val):
        """
        Will append a value to the csv_data.
        """
        if type(val) == type(self):
            for attr in ('csv_data', 'wrapped_csv', 'unwrapped_csv'):
                csv_data = getattr(self, attr)
                csv_data = csv_data.append(val.csv_data)
                csv_data.index = range(len(csv_data))
                setattr(self, attr, csv_data)
