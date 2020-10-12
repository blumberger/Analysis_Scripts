#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module holds some utilities regarding molecular system manipulation such as
reshaping atomic coords to molecular coords etc...
"""
import re
import copy
import numpy as np
from sklearn.cluster import DBSCAN
from collections import Counter

from src.system import type_checking as type_check
from src.io_utils import json_files as json

PT = json.read_json("src/data/periodic_table.json")
PT_abbrv = {PT[i]['abbreviation']: {**PT[i], **{'full_name': i}} for i in PT}


def get_atom_type(mass):
    """
    Get the atom type from the atomic mass.
    """
    for i in PT:
        mass_i = PT[i]['atomic_weight']
        if isinstance(mass_i, (float, int)):
            if np.isclose(mass_i, mass, atol=0.1):
                new_dict = copy.deepcopy(PT[i])
                new_dict['name'] = i
                return new_dict

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

def mols_to_atoms(mol_crds):
    """
    Will reshape an array of molecules to an array of atoms.

    Inputs:
        * mol_crds <array> => An array of shape (nstep, nmol, nat_per_mol, ndim)
                              or (nmol, nat_per_mol, ndim).

    Outputs:
        np.array:
            An array of shape (nstep, nat, ndim)
                           or (nat, ndim).
    """
    if len(np.shape(mol_crds)) == 3:
        nmol, nat_per_mol, ndim = mol_crds.shape
        nat = nmol * nat_per_mol
        return np.reshape(mol_crds, (nat, ndim))

    elif len(np.shape(mol_crds)) == 4:
        nstep, nmol, nat_per_mol, ndim = mol_crds.shape
        nat = nmol * nat_per_mol
        return np.reshape(mol_crds, (nstep, nat, ndim))

    else:
        raise SystemExit("Incorrect shape array should be (nmol, nat_per_mol, ndim) " +
                         "or (nstep, nmol, nat_per_mol, ndim)")

def mols_to_cols(mol_cols):
    """
    Will reshape an array of molecules to an array of atoms.

    Inputs:
        * mol_cols <array> => An array of shape (nstep, nmol, nat_per_mol)
                              or (nmol, nat_per_mol).

    Outputs:
        np.array:
            An array of shape (nstep, nat)
                           or (nat).
    """
    if len(np.shape(mol_cols)) == 2:
        nmol, nat_per_mol = mol_cols.shape
        nat = nmol * nat_per_mol
        return np.reshape(mol_cols, (nat,))

    elif len(np.shape(mol_cols)) == 3:
        nstep, nmol, nat_per_mol = mol_cols.shape
        nat = nmol * nat_per_mol
        return np.reshape(mol_cols, (nstep, nat,))

    else:
        raise SystemExit("Incorrect shape array should be (nmol, nat_per_mol, ndim) " +
                         "or (nstep, nmol, nat_per_mol, ndim)")
    

def atoms_to_mols(crds, num_ats_in_mol, cart_dims=3, nstep=1):
    """
    Will reshape array to divide up the atoms into arrays with molecules.

    i.e. reshape array from (num_ats, cart_dims) to (num_mols, ats_in_mol, cart_dims)

    Inputs:
        * crds <np.NDArray> => coordinates in shape (num_ats, cart_dims)
        * num_ats_in_mol <int> => how many atoms in a single molecule

    Outputs:
        * <np.NDArray> array of shape (num_mols, ats_in_mol, cart_dims)
    """
    if type(crds) == list:
        crds = np.array(crds)

    if len(crds.shape) == 3:
        nstep = len(crds)
        nmol = get_nmol(len(crds[0]), num_ats_in_mol)
        return np.reshape(crds, (nstep, nmol, num_ats_in_mol, cart_dims))

    elif len(crds.shape) == 2:
        nmol = get_nmol(len(crds), num_ats_in_mol)
        return np.reshape(crds, (nmol, num_ats_in_mol, cart_dims))

    else:
        raise SystemError("\n\n\nWrong shape for atomic coordinate array!\n\n"
                          + "It should either be (nstep, natom, 3) or"
                          + " (natom, 3)\n\n"
                          + f"Current shape = {crds.shape}")

def cols_to_mols(cols, num_ats_in_mol, nstep=1):
    """
    Will reshape array to divide up the atoms into arrays with molecules.

    i.e. reshape array from (num_ats, cart_dims) to (num_mols, ats_in_mol, cart_dims)

    Inputs:
        * crds <np.NDArray> => coordinates in shape (num_ats, cart_dims)
        * num_ats_in_mol <int> => how many atoms in a single molecule

    Outputs:
        * <np.NDArray> array of shape (num_mols, ats_in_mol, cart_dims)
    """
    if type(cols) == list:
        cols = np.array(cols)


    if len(cols.shape) == 2:
        nstep = len(cols)
        nmol = get_nmol(len(cols[0]), num_ats_in_mol)
        return np.reshape(cols, (nstep, nmol, num_ats_in_mol))

    elif len(cols.shape) == 1:
        nmol = get_nmol(len(cols), num_ats_in_mol)
        return np.reshape(cols, (nmol, num_ats_in_mol))

    else:
        raise SystemError("\n\n\nWrong shape for atomic coordinate array!\n\n"
                          + "It should either be (nstep, natom) or"
                          + " (natom, )\n\n"
                          + f"Current shape = {cols.shape}")


def get_K_nearest_neighbours(at, crds, K, dist=False):
    """
    Will get the K nearest neighbours for at.

    Inputs:
        * at <array 3> => The point to get the neighbours for
        * crds <array> => The other points to check.
        * K <int> => How many nearest neighbours to get
        * dist <array> OPTIONAL => The distances of the atom in question with
                                    other coords.
    """
    if dist is False:
        dist = np.linalg.norm(crds - at, axis=1)

    sorting = sorted(zip(dist, np.arange(len(dist))))[1: K+1]

    return [i[1] for i in sorting]


def substring_is_in(substring, string):
    """
    A function to check whether a subtring (containing regex) can be found in another string.

    Inputs:
        * string <str> => The string to check
        * substring <str> => The string to look for
    Outputs:
        <bool> Whether the substring is contained
    """
    matches = [i for i in re.findall(substring, string) if i]
    if len(matches): return True
    else: return False


def get_bonding_info(all_mol_crds, bond_info, cols, types, NN=False, cutoff=5):
    """
    Will determine the bonding between crds in a molecule.

    Inputs:
        * al_mol_crds <array> => All the molecular coordinates (nmol, nat_per, 3)
        * bond_info <dict> => Which atoms can bond to which key value pairs
        * cols <array> => String with what type of atom it is
        * types <array> => What element the atom is
        * NN <array<tuple>> OPTIONAL => matrix of tuples with tuples containing
                                        all (distances, at_nums, atom_type)
        * cutoff <float> OPTIONAL => The max dist that 2 atoms can bond
    Outputs:
        <list<dict>> key = atom, value = all atoms it bonds with. Duplicates are included
    """
    all_bonds = []

    # Get element types from the columns and types
    elm = np.array([[types[j] for j in i] for i in cols])

    # Loop over all molecules
    for imol, mol_crds in enumerate(all_mol_crds):
        at_bonds = {}

        # Loop over all atoms in a molecule
        for iat in range(len(mol_crds)):

            # # This atom won't bond with anything
            # if elm[imol, iat] not in bond_info:  continue

            # Get the element of the atom
            at_type = elm[imol, iat]
            at1_name = cols[imol, iat]

            # Get which closest atoms (the atoms is most likely to be bonded to).
            num_bonds = PT_abbrv[at_type]['number_bonds'][0]
            sorting = NN[iat][1:]

            inds = []
            for dist, ind, a_type in sorting:
                # If we are full up then stop the loop
                if len(inds) == num_bonds: break

                # If the elements are allowed to be bonded we may add them to the list
                if float(dist) < cutoff:

                    # If the first atom doesn't match the regex then skip it
                    for at1_regex in bond_info:
                        if substring_is_in(at1_regex, at1_name):
                            break
                    else:  continue

                    # If the second atom doesn't match any regex then skip that
                    for at2_regex in bond_info[at1_regex]:
                        if substring_is_in(at2_regex, a_type):
                            inds.append(int(ind))
                            break


            # Add the atoms to the bonding dict
            at_bonds[iat+1] = inds

        all_bonds.append(at_bonds)

    return all_bonds


def get_all_atom_chains(mol_crds, original_iat, chain, types, bond_info,
                        chain_ind=0, inds=[], new_iat=False,
                        all_inds=[]):
    """
    Will get all chains of bonded atoms that satisfy the pattern in 'chain'

    Inputs:
        * mol_crds <array> => All the coordinates of the molecule in question
        * original_iat <int> => The atom to check for the chain
        * chain <list> => A list of strings with what types of atoms are chained
                          with others
        * types <array> => The atom types for each atom in molecule
        * bond_info <dict> => All bonds each atom can form
    """
    # Initialise
    if new_iat is False:
        # Lists are sticky so reset them
        all_inds, inds = [], []
        new_iat = original_iat
    inds.append(new_iat)

    if len(inds) == len(chain):
        all_inds.append(inds)
        return all_inds

    # terminate if we violate the chain rule
    # if chain_ind == len(chain):   return all_inds
    if not substring_is_in(chain[chain_ind], types[new_iat-1]):
        return all_inds


    # Loop over each atom that the current atom is bonded to
    for iter_iat in bond_info[new_iat]:
        # If the atoms follow the chain rule and haven't been visited before
        if iter_iat in inds:  continue
        if not substring_is_in(chain[chain_ind+1], types[iter_iat-1]): continue

        # Recursively call self to iterate through all atoms
        all_inds = get_all_atom_chains(mol_crds, original_iat, chain, types,
                                      bond_info, chain_ind+1, inds[:],
                                      iter_iat, all_inds)
    return all_inds


def get_atom_masses(at_labels):
    """
    Will get atomic masses from atomic labels

    Inputs:
        * at_labels <arr> => The label of the atoms can be any shape
    Outputs:
        * <arr> The masses as floats for each label
    """
    unique_types = np.unique(at_labels)
    for elm in unique_types:
        if elm in PT_abbrv:
            at_labels[at_labels == elm] = PT_abbrv[elm]['atomic_weight']
        else:
            raise SystemError(f"I don't understand the atom type '{elm}'")

    return at_labels.astype(float)

def get_COM(all_mol_crds, mol_col):
    """
    Will calculate the Center of Mass of a list of molecules.

    N.B if the array only has 3 dimensions then it will assume that nsteps is missing

    Inputs:
        * all_mol_crds <array> => The molecular coordinate in shape (nstep, nmol, nat_per_mol, 3)
        * mol_col <array> => The type of each molecule as a string in shape (nmol, nat_per_mol)
    Outputs:
        <array> The center of masses in shape (nstep, nmol, 3)
    """
    # If we don't have the number of steps info
    add_ax = False
    if len(np.shape(all_mol_crds)) == 3:
        add_ax = True
        all_mol_crds = [all_mol_crds]

    # Get atom masses
    masses = get_atom_masses(mol_col)

    # Loop over all mols and get COM
    nstep = len(all_mol_crds)
    mass_1_mol = np.sum(masses[0])
    COMs = [[crds[:, :, 0] * masses, crds[:, :, 1] * masses, crds[:, :, 2] * masses]
             for crds in all_mol_crds]
    COMs = np.sum(COMs, axis=3) / mass_1_mol
    COMs = np.swapaxes(COMs, 1, 2)

    if add_ax:
        COMs = COMs[0]

    return COMs

def get_split_mols(all_mol_crds):
    """
    Will use stddev to detemine if a molecule is 'split' or not.
    """
    if len(np.shape(all_mol_crds)) != 4:
        raise SystemError(f"Only works for 4D data, this is {len(np.shape(all_mol_crds))}D  data. It's shape is {np.shape(all_mol_crds)}.")

    std_ats = np.std(all_mol_crds, axis=2)
    mask = np.any(std_ats > 7, axis=2)

    split_mols = np.array([step_crds[mask[i]] for i, step_crds in enumerate(all_mol_crds)])
    non_split_mols = np.array([step_crds[~mask[i]] for i, step_crds in enumerate(all_mol_crds)])
    return split_mols, non_split_mols, mask

def get_COM_split_mols(all_mol_crds, mol_col):
    """
    Will correct for the split molecules in periodically wrapped systems.

    This works by only calculating the COM of molecules that aren't split (have a stddev below 10).
    For those molecules that are split the centroid of the largest fragment is calculated.

    Inputs:
        * all_mol_crds <array> => The molecular coordinate in shape (nstep, nmol, nat_per_mol, 3)
        * mol_col <array> => The type of each molecule as a string in shape (nmol, nat_per_mol)
    Outputs:
        <array> The center of masses in shape (nstep, nmol, 3)
    """
    # If we don't have the number of steps info
    add_ax = False
    if len(np.shape(all_mol_crds)) == 3:
        add_ax = True
        all_mol_crds = np.array([all_mol_crds])

    # Get atom masses

    # Loop over all mols and get COM
    split_mols, non_split_mols, mask = get_split_mols(all_mol_crds)
    non_split_col = mol_col[~mask[0]]

    # Handle the non-split molecules
    masses = get_atom_masses(non_split_col)
    mass_1_mol = np.sum(masses[0])
    COMs = np.zeros(list(np.shape(mask)) + [3])
    tmp = [[crds[:, :, 0] * masses, crds[:, :, 1] * masses, crds[:, :, 2] * masses]
             for crds in non_split_mols]
    tmp = np.sum(tmp, axis=3) / mass_1_mol

    # If there aren't any split mols just use the normal COM
    if np.shape(split_mols)[1] == 0:
        COMs = np.swapaxes(tmp, 1, 2)

    else:
        # COMs for mols not split
        COMs[~mask] = np.swapaxes(tmp, 1, 2)

        # Now deal with the ones that have been split by the wrapping.
        COMs[mask] = [[get_largest_mol_fragment_centroid(mol_data) for mol_data in step_data] for step_data in split_mols]

    if add_ax:
        COMs = COMs[0]

    return COMs


def get_largest_mol_fragment(single_mol):
    """
    Will use a clustering algorithm (DBSCAN) to determine the largest molecular fragment.

    Inputs:
        * single_mol <array> => The coordinates of a single molecule.
    """
    db = DBSCAN(eps=2.5, min_samples=2).fit(single_mol)
    labels_mode = max(set(db.labels_), key=list(db.labels_).count)
    return single_mol[db.labels_ == labels_mode]


def get_largest_mol_fragment_centroid(single_mol):
    """
    Will use a clustering algorithm (DBSCAN) to determine the largest molecular fragment.

    Inputs:
        * single_mol <array> => The coordinates of a single molecule.
    """
    db = DBSCAN(eps=2.5, min_samples=2).fit(single_mol)
    labels_mode = max(set(db.labels_), key=list(db.labels_).count)
    return np.mean(single_mol[db.labels_ == labels_mode], axis=0)


def get_topo_info(all_mol_crds, angle_info, all_bonds, cols, types, NN=False):
    """
    Will determine the bonding between crds in a molecule.

    Inputs:
        * al_mol_crds <array> => All the molecular coordinates (nmol, nat_per, 3)
        * angle_info <dict> => Which atoms can make angles to which
        * all_bonds <dict> => Full bonding structure from 'get_b'
        * cols <array> => String with what type of atom it is
        * types <array> => What element the atom is
        * NN <array<tuple>> OPTIONAL => matrix of tuples with tuples containing
                                        all (distances, at_nums, atom_type)
        * cutoff <float> OPTIONAL => The max dist that 2 atoms can bond
    """
    all_angs = []

    # Get element types from the columns and types
    elm = np.array([[types[j] for j in i] for i in cols])

    # Loop over all molecules
    for imol, mol_crds in enumerate(all_mol_crds):
        at_angs = {iat: [] for iat in range(1, len(mol_crds)+1)}

        # Loop over all atoms in a molecule
        for iat, at in enumerate(mol_crds):
            for iang in angle_info:
                all_ang = get_all_atom_chains(mol_crds, iat+1, angle_info[iang],
                                              cols[imol], all_bonds[imol])

                for ang in all_ang:
                    at_angs[iat+1].append(ang)

        all_angs.append(at_angs)

    return all_angs
