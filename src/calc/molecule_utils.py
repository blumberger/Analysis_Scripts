#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module holds some utilities regarding molecular system manipulation such as
reshaping atomic coords to molecular coords etc...
"""

import numpy as np

from src.system import type_checking as type_check


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

    if type_check.is_int(nmol, err_msg):
        nmol = int(nmol)

    return nmol


def atoms_to_mols(crds, num_ats_in_mol, cart_dims=3):
    """
    Will reshape array to divide up the atoms into arrays with molecules.

    i.e. reshape array from (num_ats, cart_dims) to (num_mols, ats_in_mol, cart_dims)

    Inputs:
        * crds <np.NDArray> => coordinates in shape (num_ats, cart_dims)
        * num_ats_in_mol <int> => how many atoms in a single molecule

    Outputs:
        * <np.NDArray> array of shape (num_mols, ats_in_mol, cart_dims)
    """
    nmol = get_nmol(len(crds), ats_in_mol)
    return np.reshape(crds, (nmol, num_ats_in_mol, cart_dims))
