#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module that contains general types that can be subclassed to create new calculation types.
"""
import copy
import numpy as np

from collections import Counter

from src.calc import molecule_utils as mol_utils




class Calc_Type(object):
    """
    A skeletal type that can be inheritted when creating a new calc type.

    Inputs:
        * Variable <Variable> => An instance of the Variable class.

    Important Attributes:
        * required_metadata <tuple> => Any keys that are required in the metadata dictionary.
        * required_calc <tuple> => Any values that need calculating to calculate this value.
        * data <*> => The data that has been calculated.
        * required_data_types <tuple> => What classes are required to calculate the property.
                                         At the moment this only supports logical 'and's by simply
                                         adding new str entries into the tuple. The 'str's available
                                         are the keys to the load_fncs type in the input_file class.

    Public Methods:
        * calc => To be overridden to calculate the property in question.
    """
    _write_types = ('txt',)
    required_metadata = ()
    required_calc = ()
    required_data_types = ()
    _defaults = {}

    # Require these 3 objects for the formation of a new variable type
    name = "General Calc Type"
    metadata = {}
    data = 0

    # The shortcuts for more all data types that contain a more general type i.e. xyz and lammps_dump types contain position data.
    type_shortcut_dict = {'pos': ('xyz', 'lammps_dump',)}

    def __init__(self, Variable):
        """
        Just check we have the required properties for calculating the quantitiy.
        """
        self.Var = Variable

        # Copy all data to make them unique for each instance
        all_vars = [i for i in dir(self) if i[0] != '_']
        for i in all_vars:
            if i[0] != '_':
                var = getattr(self, i)
                if not callable(var) and isinstance(var, (dict, list, tuple)):
                    setattr(self, i, copy.deepcopy(var))

        # Add the required calc types to the class
        self.required_data_types = LogicalTuple([self.type_shortcut_dict.setdefault(i.lower(), i) for i in self.required_data_types])
        data_types_loaded = self.Var['data_loaded']
        if not self.required_data_types.all_in(data_types_loaded):
            err_msg = f"Can't calculate property '{self.name}' as the data"
            err_msg += f" '{str(self.required_data_types)}' haven't been loaded."
            raise SystemExit(err_msg)

        # Else just set the required data to the previous data variable.
        else:
            self.Var.required_data  = self.Var.data

        # Set the default parameters
        for key in self._defaults:
            if key not in self.Var.metadata:
                self.Var.metadata[key] = self._defaults[key]

            self.metadata[key] = self.Var.metadata[key]

        # Check we have all the data we need to calculate the property
        for key in self.required_metadata:
            if key not in self.Var.metadata:
               raise KeyError(f"Please load the data '{key}' into the variable '{self.Var.name}'")
            else:
               self.metadata[key] = self.Var[key]


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
        else:
            pass

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

    def _get_mol_crds(self):
        """
        A function to get rearrange some data to make it easier and more efficient to loop over mols.
        """
        ats_per_mol = self.metadata['atoms_per_molecule']
        if 'cols' not in dir(self):
            self.get_cols(self.metadata['number_each_atom'],
                          ats_per_mol)

        self.natom = len(self.compute_data[0])
        all_mol_crds = mol_utils.atoms_to_mols(self.compute_data, ats_per_mol)
        self.nmol = all_mol_crds.shape[1]
        mol_col = np.reshape(self.cols[0], (self.nmol, ats_per_mol))
        return all_mol_crds, mol_col

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


class LogicalTuple(tuple):
    """
    A subclass of the tuple object to make it slightly easier to handle certain objects.

    This can also be used to check if 2 tuples are logicallly equivalent assuming all entries
    along axis 0 of the container are 'and' conditions and all values along axis 1 are 'or' values.
    E.g:
        ('a', 'b', ('c', 'd')) is equivalent to saying 'a' and 'b' and ('c' or 'd'). Here this is
        logically equivalent to the tuples ('a', 'b', 'c') and ('a', 'b', 'd') (where or is inclusive).
    """

    def __new__(cls, data):
        try:
            return super().__new__(cls, tuple(data))
        except TypeError:
            return super().__new__(cls, (data, ))

    def __init__(self, data):
        """
        Run some checks on the tuple to see if the input is OK.

        Should also check it's dimensionality too... It can be 1D or 2D.
        """
        # Set all lists to tuples
        try:
            data = list(data)
        except TypeError:
            print("\n\nPlease enter an iterable to the class LogicalTuple().")
            print(f"You entered {data}, did you forget to use a comma in your tuple?")
            print(f"E.g. ({data}) instead of ({data}, )..."+"\n\n")
            raise TypeError(f"'{type(data)}' not iterable.")

        for i, val in enumerate(data):
            if type(val) == list:
                data[i] = tuple(val)
            elif not isinstance(val, (tuple, int, float, str)):
                raise TypeError("Incorrect use of LogicalTuple type.")

        self.data = tuple(data)
        self.all_tuples = tuple(i for i in self.data if type(i) == tuple)
        self.all_single_vals = tuple(i for i in self.data if type(i) != tuple)

        # Create the string from the data
        self.str_data = ""
        for ind, i in enumerate(self.data):
            if type(i) != tuple:
                self.str_data += f"{i}"
            else:
                joined_str = " and/or ".join([f"'{str(j)}'" for j in i])
                self.str_data += f"({joined_str})"

            if ind < len(self.data)-1:
                self.str_data += " and "

    def __contains__(self, key):
        """
        Overload 'in' operator to check whether something is in any of the tuple (it can be 2D).

        This will check if the input 'key' is found within the tuple anywhere
        E.g. if the tuple in question is (1, 2, 3, (4, 5)) then any tuple
        containing any length of combination of 1, 2, 3, 4, 5 will be classed
        as being contained within the tuple so (1, 2, 3) or (1, 4, 5) or (1, 1, 2) etc...
        """
        try:
            '1' in key
            for i in key:
                if i in self.data:
                    continue
                if any(i in j for j in self.all_tuples):
                    continue
                else:
                    return False

        except TypeError:
            for i in self.__iter__():
                if isinstance(i, (tuple, list)):
                    if key in i:
                        return True
                elif key == i:
                    return True
            return False

        return True

    def all_in(self, compare):
        """
        Check all of the tuple is found in the input.

        Will return True if all the tuple saved in this class can be found
        within the input. However, all of the input doesn't need to be
        found in the tuple.

        For example, if the tuple saved in this class is (1, 2, 3, (4, 5)) then
        any inputs containing all of 1 and 2 and 3 and (4 and/or 5) and any other
        object will pass the test.

        This differs from the equality operator '=' in the fact it allows the
        input to contain extra objects.
        """
        # Check all the ANDs
        for i in self.all_single_vals:
            if i in compare: continue
            else: return False

        # Check all the ORs
        for i in self.all_tuples:
            if not any(j in compare for j in i):
                return False

        return True

    def __eq__(self, comp):
        """
        Overload the = operator to check if 2 tuples are "logically equivalent".

        This is a more stringent requirement than used for the 'in' operator.
        If all values in tuple are found within the comparator and all values
        within the comparator are found within the tuple then this will be
        satisfied assuming the and/or rules mentioned in the main docstring.

        E.g. if the tuple saved in the class is (1, 2, 3, (4, 5)) then any
        tuple containing 1 and 2 and 3 and (4 and/or 5) will be allowed.
        """
        test_list = list(comp)

        for i in self.data:

            # Check the ORs
            if type(i) == tuple:
                for j in i:
                    if j in test_list:
                        for k in i:
                            if k in test_list:
                                test_list.remove(k)
                        break
                else:
                    return False

            # Check the ANDs
            elif i in test_list:
                test_list.remove(i)
                continue

            else:
                print(f"NO AND {i}")
                return False

        return not len(test_list)

    def __ne__(self, comp):
        """Will just return the inverse of the equal operator."""
        return not self.__eq__(comp)

    def __str__(self):
        """Will overload the str() function."""
        return self.str_data
