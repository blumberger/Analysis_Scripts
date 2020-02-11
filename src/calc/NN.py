#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate nearest neighbour lists
"""
import numpy as np

from src.calc import general_types as gen_type
from src.system import type_checking as type_check

class NN(gen_type.Calc_Type):
    """
    Will calculate the Nearest neighbour list from the data contained within a data file.

    Inputs:
        * Variable <Variable> => A ass or a class that has been derived from Data_File

    Important Attributes:
        * required_metadata <tuple> => Any keys that are required in the metadata dictionary.
        * required_calc <tuple> => Any values that need calculating to calculate this value.
        * data <*> => The data that has been calculated.
    """
    required_metadata = ()

    # Need these 3 attributes to create a new variable type
    data = {}
    metadata = {'file_type': 'json'}
    name = "Nearest Neighbour Calculator"

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
        XYZFile = self.Var.data
        cols = XYZFile.cols
        at_crds = np.array([i[cols[0] != 'Ne'] for i in XYZFile.numeric_data])
        self.natom = len(at_crds[0])
        self.nstep = len(at_crds)
       
        # Calculate the nearest neighbour lists for each step
        for step in range(self.nstep):
            self.data[step] = {}

            # Get coords
            crds = at_crds[step]

            # Get distances between neighbours
            self.__get_distances(crds)

            # Get a sorted list of atom indices by distance
            self.__get_nearest_atom_inds()

            # If we have some molecule metadata
            if 'num_at_per_mol' in self.Var.metadata:
                self.__get_nmol()
                self.__reshape_at_dist()
                self.__get_nearest_atom_inds_per_mol()
                self.data[step]['closest_atoms_mol_grouped'] = self.closest_at_per_mol.tolist()
                self.data[step]['distances_mol_grouped'] = self.all_dist_per_mol.tolist()
                 
            
            # Save data in dict (use tolist for writing later)
            self.data[step]['distances'] = self.all_dist.tolist()
            self.data[step]['closest_atom_indices'] = self.closest_ats.tolist()
            
        return self.data

    def __get_distances(self, crds):
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

    def __get_nearest_atom_inds(self):
        """
        Will sort the atom indices by the distance from the current atom.

        This will create a numpy array of size (natom, natom). Indexing as:
             closest_ats[i]
        will give a list of atom indices sorted by distance from atom i.

        The data is stored as self.closest_ats
        """
        # Create empty data structure
        self.closest_ats = np.zeros((self.natom, self.natom-1))

        # Get and sort distances
        all_at_inds = np.arange(self.natom)
        for iat in range(self.natom):
            at_inds = all_at_inds[all_at_inds != iat]
            dist = self.all_dist[iat, at_inds]
            
            at_inds = [i[1] for i in sorted(zip(dist, at_inds))]
            self.closest_ats[iat] = at_inds

    def __get_nmol(self):
        """
        Will get the number of molecules based on the sizes of various arrays.
        """
        self.at_per_mol = self.Var.metadata['num_at_per_mol']
        self.nmol = self.natom / self.at_per_mol
 
        # Error checking for self.nmol
        err_msg = "Number of atoms per molecule doesn't neatly divide up the atoms in each xyz step!"
        if type_check.is_int(self.nmol, err_msg):  self.nmol = int(self.nmol)
 
 
    def __reshape_at_dist(self):
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
 
    def __get_nearest_atom_inds_per_mol(self):
        """
        Will sort the atom indices by the distance from the current atom for each molecule.
 
        Will create array of size:
            (nmol, nat_per_mol, nat_per_mol)
 
        indexing the array as: arr[i][j] will return a list of atom indices for molecule i,
        with elements sorted by distance from atom j.
        """
        self.closest_at_per_mol = np.zeros((self.nmol,
                                            self.at_per_mol,
                                            self.at_per_mol-1))

        # Get and sort distances
        all_at_inds = np.arange(self.at_per_mol)
        for imol in range(self.nmol):
           for iat in range(self.at_per_mol):
               at_inds = all_at_inds[all_at_inds != iat]
               dist = self.all_dist_per_mol[imol, iat, at_inds]
            
               at_inds = [i[1] for i in sorted(zip(dist, at_inds))]
               self.closest_at_per_mol[imol, iat] = at_inds
 
