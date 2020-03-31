#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate nearest neighbour lists
"""
import numpy as np
import time

from src.calc import general_calc as gen_calc
from src.calc import molecule_utils as mol_utils

from src.system import type_checking as type_check

class NN(gen_calc.Calc_Type):
    """
    Will calculate the Nearest neighbour list from the data contained within a data file.

    The nearest neighbour list will be stored in a dictionary named 'data'. The
    structure of this dictionary is:
        {istep:
               {
                'distances': np.NDArray(natom, natom),
                'closest_atom_indices': np.NDArray(natom, natom),
               }
         }

    If the number of atoms per molecule is specified then 2 more keys will
    appear at each step:
        {istep:
            {
             'distances': np.NDArray(natom, natom),
             'closest_atom_indices': np.NDArray(natom, natom),
             'closest_atoms_mol_grouped': np.NDArray(nmol, natom, natom-1),
             'distances_mol_grouped': np.NDArray(nmol, natom, natom-1),
            }
        }

    Inputs:
        * Variable <Variable> => An instance of the Variable class

    Important Attributes:
        * required_metadata <tuple> => Any keys that are required in the metadata dictionary.
        * required_calc <tuple> => Any values that need calculating to calculate this value.
        * data <*> => The data that has been calculated.
    """
    required_metadata = ()
    required_data_types = ('pos',)

    # Need these 3 attributes to create a new variable type
    metadata = {'file_type': 'json'}
    name = "Nearest Neighbour Lists"

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
        self.data = []
        all_xyz_data = self.Var.data.get_xyz_data()
        all_cols = self.Var.data.get_xyz_cols()

        # Loop over all the xyz data and cols we have
        for xyz_data, cols in zip(all_xyz_data, all_cols):

            at_crds = np.array([i[cols[0] != 'Ne'] for i in xyz_data])
            self.natom = len(at_crds[0])
            self.nstep = len(at_crds)
            self.step_data = {}

            # Calculate the nearest neighbour lists for each step
            for step in range(self.nstep):
                self.step_data[step] = {}

                # Get coords
                crds = at_crds[step]

                # Get distances between neighbours
                self.get_distances(crds)

                # Get a sorted list of atom indices by distance
                self.get_nearest_atom_inds()

                # If we have some molecule metadata
                if 'atoms_per_molecule' in self.Var.metadata:
                    self.at_per_mol = self.Var.metadata['atoms_per_molecule']
                    self.nmol = mol_utils.get_nmol(self.natom, self.at_per_mol)
                    self.reshape_at_dist()
                    self.get_nearest_atom_inds_per_mol()
                    self.step_data[step]['closest_atoms_mol_grouped'] = self.closest_at_per_mol
                    self.step_data[step]['distances_mol_grouped'] = self.all_dist_per_mol

                # Save data in dict
                self.step_data[step]['distances'] = self.all_dist
                self.step_data[step]['closest_atom_indices'] = self.closest_ats

            self.data.append(self.step_data)

        return self.data

    def get_distances(self, crds):
        """
        Will get all distances between atoms and store them in a 2D numpy array.

        A numpy array of size (natom, natom) is created and indexing can be done by:
            all_dist[i, j]
        This will give the distance between atom i and j.

        The data is stored as self.all_dist.

        Inputs:
            * crds <np.NDArray> => The atomic coordinates to find distances of
        """
        self.all_dist = np.zeros((self.natom, self.natom))
        # Loop over upper triangle of atom pairs
        for iat in range(self.natom-1):
            # Get the atom indices
            at_inds = np.arange(len(crds))

            # Calc distances between atoms (only upper triangle though)
            at_msk = at_inds > iat
            all_ut_dist = crds[at_msk] - crds[iat]
            all_ut_dist = np.linalg.norm(all_ut_dist, axis=1)

            self.all_dist[iat, iat+1:] = all_ut_dist

        # Get lower triangle indices
        self.all_dist = self.all_dist + self.all_dist.T

    def get_nearest_atom_inds(self):
        """
        Will sort the atom indices by the distance from the current atom.

        This should be optimised!

        This will create a numpy array of size (natom, natom). Indexing as:
             closest_ats[i]
        will give a list of atom indices sorted by distance from atom i.

        The data is stored as self.closest_ats
        """
        # Create empty data structure
        self.closest_ats = np.zeros((self.natom, self.natom-1), dtype=int)

        # Get and sort distances
        all_at_inds = np.arange(self.natom)
        for iat in range(self.natom):
            at_inds = all_at_inds[all_at_inds != iat]
            dist = self.all_dist[iat, at_inds]

            at_inds = [i[1] for i in sorted(zip(dist, at_inds))]
            self.closest_ats[iat] = at_inds


    def reshape_at_dist(self):
        """
        Will reshape the all_dist array if we have the number of atoms per molecule metadata.

        Will reshape all_dist into:
             (nmol, nat_per_mol, nat_per_mol)
        """
        self.all_dist_per_mol = np.zeros((self.nmol, self.at_per_mol, self.at_per_mol))
        for imol in range(self.nmol):
           start, end = self.at_per_mol*imol, (imol+1)*self.at_per_mol
           self.all_dist_per_mol[imol] = self.all_dist[start:end,
                                                               start:end]

    def get_nearest_atom_inds_per_mol(self):
        """
        Will sort the atom indices by the distance from the current atom for each molecule.

        Will create array of size:
            (nmol, nat_per_mol, nat_per_mol)

        indexing the array as: arr[i][j] will return a list of atom indices for molecule i,
        with elements sorted by distance from atom j.
        """
        self.closest_at_per_mol = np.zeros((self.nmol,
                                            self.at_per_mol,
                                            self.at_per_mol-1), dtype=int)

        # Get and sort distances
        all_at_inds = np.arange(self.at_per_mol)
        for imol in range(self.nmol):
           for iat in range(self.at_per_mol):
               at_inds = all_at_inds[all_at_inds != iat]
               dist = self.all_dist_per_mol[imol, iat, at_inds]

               at_inds = [i[1] for i in sorted(zip(dist, at_inds))]
               self.closest_at_per_mol[imol, iat] = at_inds

    def json_data(self):
        """
        Will return data in a form that the json writer can write.

        This function is just used for writing
        """
        return_data = []

        # Loop over num files
        for ifile in range(len(self.data)):

            # Loop over steps
            file_return_data = {}
            for istep in range(len(self.data[ifile])):
                file_return_data[istep] = {}
                for key in self.data[ifile][istep]:
                    file_return_data[istep][key] = self.data[ifile][istep][key].tolist()

            return_data.append(file_return_data)

        return return_data

    def __str__(self):
        """
        Overload the string function
        """
        for step in self.data:
            for key in self.data[step]:
                self.data[step][key] = self.data[step][key]

        return str(self.data)
