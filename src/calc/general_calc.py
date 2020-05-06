#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module that contains general types that can be subclassed to create new calculation types.
"""
import copy
import matplotlib.pyplot as plt
import numpy as np

from collections import Counter

from src.calc import molecule_utils as mol_utils

from src.input_file import function_dicts as f_dicts


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
        if type(data_types_loaded) == bool or not self.required_data_types.all_in(data_types_loaded):
            err_msg = f"Can't calculate property '{self.name}' as the data"
            err_msg += f" '{str(self.required_data_types)}' haven't been loaded."
            err_msg += "\n\n" + f"Data types loaded = {data_types_loaded}" + "\n"
            raise SystemError(err_msg)

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

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

    def _calc_(self):
        """
        A function to override with the function to calculate the required property
        """
        print("\n\n" +f"Please overload the _calc_ function in {self.name}" + "\n\n")
        raise SystemExit("General Calc _calc_ function not overloaded!")

    def calc(self):
        """
        A placeholder method. In child classes this should be replaced with the appropriate calc
        function.
        """
        # Handle the metadata used to calculate things
        curr_metadata = self.Var.metadata.copy()

        for key in self.Var.metadata:
            if key not in self.metadata:
                self.metadata[key] = self.Var.metadata[key]

        # Calculate and prerequisites and save them in the new object
        for calc in self.required_calc:
            print(f"Calculating {calc}")           
            Calculated_Object = f_dicts.calc_fncs[calc](self.Var)
            Calculated_Object.calc()
            setattr(self, calc, Calculated_Object)

        # Preserve the metadata from before the calculation.
        self.Var.metadata_update(curr_metadata)
        self._calc_()

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

    def _plot_xyz_data(self, ats, a=False, args={'ls': 'None', 'marker': '.'}):
        if a is False:
            f = plt.figure()
            a = f.add_subplot(111, projection="3d")

            a.set_xlabel("X")
            a.set_ylabel("Y")
            a.set_zlabel("Z")

            a.view_init(elev=-3, azim=64)

            a.set_xticks([])
            a.set_yticks([])
            a.set_zticks([])

        a.plot(ats[:, 0], ats[:, 1], ats[:, 2], **args)

        return a

    def _rotate_plot_(self, ax, elev_increment=False, azim_increment=1, init_elev=30, init_azim=0, step_func=False, args=()):
        """
        Will rotate a 3D plot.

        Inputs:
            * ax <plt.axis> => The axis to rotate
            * elev_increment <float> => The value to increment the elevation axis by
            * azim_increment <float> => The value to increment the azimuthal axis by
            * init_elev <float> => The initial elevation of the axis
            * init_azim <float> => The initial azimuth angle of the axis
            * step_func <function> => A function to carry out at each step
        """
        count = 0
        if elev_increment and azim_increment:
            for elev_ang in np.arange(init_elev, 360 + init_elev, elev_increment):
                for azim_ang in np.arange(init_azim, 360 + init_azim, azim_increment):
                    count += 1
                    ax.view_init(azim=azim_ang, elev=elev_ang)
                    if step_func:
                        step_func(count, *args)

        elif elev_increment and not azim_increment:
            for elev_ang in np.arange(init_elev, 360 + init_elev, elev_increment):
                count += 1
                ax.view_init(azim=init_azim, elev=elev_ang)
                if step_func:
                    step_func(count, *args)

        elif azim_increment and not elev_increment:
            for azim_ang in np.arange(init_azim, 360 + init_azim, azim_increment):
                count += 1
                ax.view_init(azim=azim_ang, elev=init_elev)
                if step_func:
                    step_func(count, *args)


    ########################################################
    ## For file writing
    def get_csv_data(self):
        """
        Just a dummy function to return a useful error message to the user.

        When implemented in children it will return the csv data as a pandas
        DataFrame so the Write_CSV class can write a csv.
        """
        if 'csv_data' not in dir(self):
            raise SystemError(f"No function 'get_csv_data()' in the '{self.name}' class.")

        return self.csv_data

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

    def __setitem__(self, key, val):
        """Will check if there is a data variable if not return error"""
        try:
            self.data[key] = val
        except AttributeError:
            raise SystemError(f"No variable 'data' implemented in '{self.name}'")
        except:
            raise SystemError(f"Can't index '{self.name}.data' ({type(self.data)}) with index '{key}' ({type(key)})")

    def __getitem__(self, key):
        """Will check if there is a data variable if not return error"""
        try:
            return self.data[key]
        except AttributeError:
            raise SystemError(f"No variable 'data' implemented in '{self.name}'")
        except:
            raise SystemError(f"Can't index '{self.name}.data' ({type(self.data)}) with index '{key}' ({type(key)})")


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
