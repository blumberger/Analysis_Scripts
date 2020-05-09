#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to read and store the Pseudo_Hamiltonian
"""

import pandas as pd
import numpy as np

from src.io_utils import general_io as gen_io

from src.system import type_checking as type_check


class Pseudo_Ham(gen_io.DataFileStorage):
    """
    Will read the Pseudo Hamiltonian.

    A class that will loop over all lines and parse the CSVs from the
    lammps log file.
    """
    csv_data = []
    metadata = {}
    name = "Pseudo Hamiltonian"
    _write_types = ('txt',)

    def __init__(self, filepath):
        super().__init__(filepath)


    def _parse_(self):
        """
        Will parse the pseudo hamiltonian file.
        """
        ltxt = self.file_txt.split("\n")
        # self.get_step_data()

        self.data = []
        counts = {}

        data = {'title_lines': [], 'couplings': [], 'site_energies': []}
        start = False
        for line in ltxt:

            words = line.split()
            if len(words) == 3:
                mol1, mol2, Hval = words
                mol1, mol2, Hval = int(mol1), int(mol2), float(Hval)
                start = True

                # Minus 1 for tranlating from 0 indexing (python) to 1 indexing (FORTRAN)
                data.setdefault(mol1 - 1, {})[mol2 - 1] = Hval

                if mol1 != mol2:
                    data['couplings'].append(Hval)
                else:
                    data['site_energies'].append(Hval)

            else:
                if start is False:
                    data['title_lines'].append(line)

                else:
                    start = False
                    self.data.append(data)
                    data = {'title_lines': [line], 'couplings': [], 'site_energies': []}

    def get_data(self):
        """
        Will just return the data to make so the user doesn't need to know the name of it in the class.
        """
        return self.data

    def append(self, val):
        if type(val) == type(self):
            for i in range(len(val.data)):
                if i < len(self.data):
                    new_dict = val.data[i]
                    curr_dict = self.data[i]
                    for key in new_dict:
                        if key not in curr_dict:
                            curr_dict[key] = new_dict[key]
                        else:
                            print(f"Duplicate Entries in both dicts: {key}")
            
        else:
            raise TypeError(f"\n\nCan't append object of type {type(val)} to {self.name}.")
