#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Some useful utilities for the C wrapper functions.
"""

import numpy as np

def to_list(arr):
	"""
	If the input is a numpy array then convert it to a list.

	Inputs:
		* arr <array | list> => The input to convert
	Ouptuts:
		<list> The list of the input
	"""
	if type(arr) == type(np.array(1)):
		return arr.tolist()
	elif type(arr) == list:
		return arr
	else:
		return list(arr)


def is_int(num, msg=False):
    """
    Will check whether a parameter is a float or not.

    If a msg is supplied a TypeError will be thrown if the float isn't an int.

    Inputs:
        * num <float> => Number to check
        * msg <str> => The msg to report if it isn't an int
    Outputs:
        <bool>
    """
    if not num.is_integer():
        if msg is not False:
            raise TypeError(msg)
        else:
            return False
    else:
        return True


def get_nmol(natom, natom_in_mol):
    """
    Will divide a number of atoms into molecules and raise an error if it not.

    Inputs:
        * natom <int> => The number of atoms in the system
        * natom_in_mol <int> => The number of atoms in 1 molecule.
    Outputs:
        <int> Number of molecules.
    """
    nmol = natom / natom_in_mol
    err_msg = "\n\n\nPlease double check that you have entered the correct "
    err_msg += "number of atoms per molecule! \nThe number of atoms per "
    err_msg += "molecule does give an integer number of molecules!\n"
    err_msg += f"Number atoms = {natom}"
    err_msg += "\n"
    err_msg += f"Number atoms per molecule = {natom_in_mol}"

    if is_int(nmol, err_msg):
        nmol = int(nmol)

    return nmol