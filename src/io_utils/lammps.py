#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to store objects that perform I/O tasks on lammps files.
"""

import pandas as pd
import numpy as np
from io import StringIO

from src.io_utils import general_io as gen_io


class Lammps_File(gen_io.DataFileStorage):
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
        

def write_lammps_log_CSVs(Lammps_File, filepath):
    """
    Will write the different CSVs for each csv in the Lammps log file.

    Inputs:
        * Lammps_File <Lammps_File> => An instance of the Lammps_File class.
    """
    for df_num, df in enumerate(Lammps_File.data):
        filename, ext = gen_io.remove_file_extension(filepath)
        print(filename, ext)

