#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module that contains general types that can be subclassed to create new calculation types.
"""
import copy
import numpy as np

from collections import Counter

class Calc_Type(object):
    """
    A skeletal type that can be inheritted when creating a new calc type.

    Inputs:
        * Variable <Variable> => An instance of the Variable class.

    Important Attributes:
        * required_metadata <tuple> => Any keys that are required in the metadata dictionary.
        * required_calc <tuple> => Any values that need calculating to calculate this value.
        * data <*> => The data that has been calculated.

    Public Methods:
        * calc => To be overridden to calculate the property in question.
    """
    _write_types = ('txt',)
    required_metadata = ()
    required_calc = ()
    required_data_names = ()
    _defaults = {}

    # Require these 3 objects for the formation of a new variable type
    name = "General Calc Type"
    metadata = {}
    data = 0

    def __init__(self, Variable):
        """
        Just check the required metadata is there and call the calc function.
        """
        # Copy all data to make them unique for each instance
        all_vars = [i for i in dir(self) if i[0] != '_']
        for i in all_vars:
            if i[0] != '_':
                var = getattr(self, i)
                if not callable(var) and isinstance(var, (dict, list, tuple)):
                    setattr(self, i, copy.deepcopy(var))

        # Check we have all the data we need to calculate the property
        self.Var = Variable
        for key in self.required_metadata:
            if key not in self.Var.metadata:
               raise KeyError(f"Please load the data '{key}' into the variable '{self.Var.name}'")
            else:
               self.metadata[key] = self.Var[key]

        # Set the default parameters
        for key in self._defaults:
            if key not in self.metadata:
                if key not in self.Var.metadata:
                    self.metadata[key] = self._defaults[key]
                else:
                    self.metadata[key] = self.Var.metadata[key]

        for name in self.required_data_names:
            all_attrs = dir(self.Var.data)
            if name not in all_attrs and f"{name}_data" not in all_attrs:
                var_name = self.Var.name
                raise AttributeError("\n\n"
                                     + f"I can't calculate {self.name} from the"
                                     + f" variable '{var_name}'."
                                     + "\n\nYou need to give me a"
                                     + f" variable with some '{name}' data in.")

    def get_xyz_data(self):
        """
        Will get the data to use from the inputted class.
        """
        if 'xyz_data' in dir(self.Var.data):
            self.compute_data = self.Var.data.xyz_data
        elif 'csv_data' in dir(self.Var.data):
            self.Var.data.set_data()
            self.compute_data = self.Var.data.csv_data[['x', 'y', 'z']].to_numpy()
            self.compute_data = np.array([self.compute_data])

    def get_cols_from_CSV(self, number_each_atom, ats_per_mol):
        """
        Will get the element type columns from a csv file.
        """
        df = self.Var.data.csv_data

        # Error check
        if 'type' not in df.columns:
            raise SystemError("\n\n\nI can't calculate the COMs of the mols. I "
                              + "can't find the atoms types in the data."
                              + "I don't what atom types are in the mol")

        # Error check -if there is a mol with equal nums of atom of x and y type.
        num_elm_in_mol = [number_each_atom[i] for i in number_each_atom]
        if len(set(num_elm_in_mol)) != len(num_elm_in_mol):
            raise SystemError("\n\n\nCan't compute center of masses for RDF."
                            + " I don't know what types the atoms are and can't"
                            + " work them out as the molecule has 2 or more "
                            + "elements with the same number of atoms.\n\n"
                            + " There are no other RDF methods implemented.")

        # Compare how many atoms of each type are in each molecule and how many there should be.
        num_at_types = Counter(df.loc[:ats_per_mol-1, 'type'])
        at_types = {i: number_each_atom[i] for i in number_each_atom}
        cvt_type = {}
        for i in num_at_types:
            for elm in at_types:
                if at_types[elm] == num_at_types[i]:
                    print(f"Atom type '{i}' is element {elm}.")
                    cvt_type[str(i)] = elm
                    break
            else:
                raise SystemError("Could determine what element each atom type is.")
            at_types.pop(elm)

        # Create the cols array
        self.cols = df['type'].to_numpy().astype(str)
        self.natom = len(self.cols)
        for i in np.unique(self.cols):
            self.cols[self.cols == i] = cvt_type[i]
        self.cols = np.array([self.cols])

    def get_cols(self, number_each_atom, ats_per_mol):
        """
        Will try to get the atom element types from the number of each type in a
        molecule.
        """

        if 'cols' not in dir(self.Var.data):
            if 'csv_data' in dir(self.Var.data):
                self.get_cols_from_CSV(number_each_atom, ats_per_mol)
            else:
                raise SystemError("Can't compute center of masses for RDF."
                               + " There are no other RDF methods implemented.")
        else:
            self.cols = self.Var.data.cols

    def calc(self):
        """
        A placeholder method. In child classes this should be replaced with the appropriate calc
        function.
        """
        print("Please override the 'calc' method in the class '{self.__name__}'")

    # Overload the type convertors
    def __str__(self):
        """return str(data)"""
        return str(self.data)
    def __float__(self):
        """return float(data)"""
        return float(self.data)
    def __int__(self):
        """return int(data)"""
        return int(self.data)

    # Overload mathematical operators
    def __add__(self, val):
        """Ammend data attribute and return self"""
        self.data += val
        return self
    def __sub__(self, val):
        """Ammend data attribute and return self"""
        self.data -= val
        return self
    def __mul__(self, val):
        """Ammend data attribute and return self"""
        self.data *= val
        return self
    def __truediv__(self, val):
        """Ammend data attribute and return self"""
        self.data /= val
        return self
    def __floordiv__(self, val):
        """Ammend data attribute and return self"""
        self.data //= val
        return self
    def __pow__(self, val):
        """Ammend data attribute and return self"""
        self.data **= val
        return self
