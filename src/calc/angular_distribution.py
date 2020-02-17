#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate nearest neighbour lists
"""
import numpy as np

from src.calc import general_types as gen_type
from src.system import type_checking as type_check

class Angular_Dist(gen_type.Calc_Type):
    """
    Will calculate the angular distributions of the molecules in a system.

    Inputs:
        * Variable <Variable> => An instance of the Variable class

    Important Attributes:
        * required_metadata <tuple> => Any keys that are required in the metadata dictionary.
        * required_calc <tuple> => Any values that need calculating to calculate this value.
        * data <*> => The data that has been calculated.
    """
    required_metadata = ('long_axis_atoms', 'short_axis_atoms')

    # Need these 3 attributes to create a new variable type
    data = {}
    metadata = {'file_type': 'json'}
    name = "Angular Distributions"

    def calc(self):
        """
        Will calculate the angular distribution of the molecular system.

        This will loop over each molecule, find the rotation (wrt x,y,z) of the
        long and short axes of the molecule then create a histogram of this
        data.
        """
        pass
