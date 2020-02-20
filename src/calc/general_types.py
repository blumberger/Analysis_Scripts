#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module that contains general types that can be subclassed to create new calculation types.
"""

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
        # Check we have all the data we need to calculate the property
        self.Var = Variable
        for key in self.required_metadata:
            if key not in self.Var.metadata:
               raise KeyError(f"Please load the data '{key}' into the variable '{self.Var.name}'")

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
