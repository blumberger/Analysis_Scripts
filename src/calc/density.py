#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate nearest neighbour lists
"""
import numpy as np

from src.calc import general_types as gen_type
from src.system import type_checking as type_check

class Density(gen_type.Calc_Type):
    """
    Will calculate the Density from the data contained within a data file.

    Inputs:
        * Variable <Variable> => An instance of the Variable class

    Important Attributes:
        * required_metadata <tuple> => Any keys that are required in the metadata dictionary.
        * required_calc <tuple> => Any values that need calculating to calculate this value.
        * data <*> => The data that has been calculated.
    """
    required_metadata = ()

    # Need these 3 attributes to create a new variable type
    data = {}
    metadata = {'file_type': 'json'}
    name = "Nearest Neighbour Calculator"

    def calc(self):
        """
        Will calculate the nearest neighbour list from an xyz file.

        This function works by looping over all steps and then over atoms. The
        distances between the atom and all atoms with indices higher than itself
        are calculated (in order to avoid pointless calculations).

        The distances and indices of the atoms at a certain distance are then
        stored in the self.data dict with the first key indicating the step and
        the second key being either 'distances' or 'atom_indices' the third key
        will be the atom index.
        """
        data = self.Var.data.data
        did_dens_calc = False
        data_count = 0

        if type(data) == list:
            for df in data:
                if all(j in df.columns for j in ('Lx', 'Ly', 'Lz',)):
                    self.data[data_count] = self.__calc_dens__(df)
                    did_dens_calc = True
                    data_count += 1

        elif type(data) == pd.DataFrame:
            if all(j in df.columns for j in ('Lx', 'Ly', 'Lz',)):
                self.data[data_count] = self.__calc_dens__(data)
                did_dens_calc = True

        if did_dens_calc is False:
            raise SystemExit("Can't find the required DataFrame headers to calculate the density")

    def __calc_dens__(self, df):
        """
        Will calculate the density from a single dataframe
        """
        vol = df['Volume']
        print(vol)

    def __str__(self):
        """
        Overload the string function
        """
        for step in self.data:
            for key in self.data[step]:
                self.data[step][key] = self.data[step][key].tolist()

        return str(self.data)
