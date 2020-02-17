#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate nearest neighbour lists
"""
import numpy as np
import pandas as pd

from src.calc import general_types as gen_type
from src.calc import molecule_utils as mol_utils
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
    required_metadata = ('molecular_mass', "atoms_per_molecule", "number_atoms")
    required_data_names = ('csv', )

    # Need these 3 attributes to create a new variable type
    csv_data = {}
    metadata = {'file_type': 'csv'}
    name = "Densities"

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
        data = self.Var.data.csv_data
        did_dens_calc = False
        data_count = 0

        self.__get_tot_mass__()

        if type(data) == list:
            for df in data:
                if all(j in df.columns for j in ('Lx', 'Ly', 'Lz',)):
                    self.csv_data[data_count] = self.__calc_dens__(df)
                    did_dens_calc = True
                    data_count += 1

        elif type(data) == pd.DataFrame:
            if all(j in df.columns for j in ('Lx', 'Ly', 'Lz',)):
                self.csv_data[data_count] = self.__calc_dens__(data)
                did_dens_calc = True

        if did_dens_calc is False:
            raise SystemExit("Can't find the required DataFrame headers to calculate the density")

    def __get_tot_mass__(self):
        """
        Will get the total mass of the system in kg.

        First we find the number of molecules then the total mass.
        """
        # Get num mol
        nmol = mol_utils.get_nmol(self.Var.metadata['number_atoms'],
                                  self.Var.metadata['atoms_per_molecule'])
        self.tot_mass = nmol * self.Var.metadata['molecular_mass']

    def __calc_dens__(self, df):
        """
        Will calculate the density from a single dataframe

        We assume the mass is in kg and the volume is in A^3

        This will calculate densities in g/cm^3
        """
        new_df = pd.DataFrame({})
        # Convert to grams
        tot_mass = self.tot_mass * 1e3
        if 'Volume' in df.columns:
            vol = df['Volume'] * 1e-24
        else:
            vol = df['Lx'] * df['Ly'] * df['Lz'] * 1e-24

        if 'Step' in df:  new_df['Step'] = df['Step']
        new_df['Volumes'] = vol
        new_df['Density'] = tot_mass / vol

        return new_df

    def __str__(self):
        """
        Overload the string function
        """
        s = ""
        for key in self.csv_data:
            s += f"Dataframe {key}:"
            s += "\n\n" + str(self.csv_data[key])

        return str(self.csv_data)
