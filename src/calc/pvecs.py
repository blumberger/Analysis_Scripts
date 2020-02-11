#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate pvecs from some xyz coordinates.
"""
import numpy as np

from src.calc import general_types as gen_type
from src.system import type_checking as type_check

class PVecs(gen_type.Calc_Type):
    """
    Will calculate the pvecs from the data contained within a data file.

    Inputs:
        * Variable <Variable> => A ass or a class that has been derived from Data_File

    Important Attributes:
        * required_metadata <tuple> => Any keys that are required in the metadata dictionary.
        * required_calc <tuple> => Any values that need calculating to calculate this value.
        * data <*> => The data that has been calculated.
    """
    required_metadata = ('num_at_per_mol',)
    required_calc = ('NN', )

    # Need these 3 attribute to create a new variable type
    metadata = {'file_type': 'xyz'}
    name = "P-Vec Calculator"
    data = []

    def calc(self):
        """
        Will calculate the pvecs from self.Data_File.numeric_data
        """
        XYZFile = self.Var.data
        at_crds = XYZFile.numeric_data
        self.nstep = XYZFile.nstep
        self.at_per_mol = self.Var.metadata['num_at_per_mol']
        
        # First remove 'Ne' atoms
        self.cols = XYZFile.cols
        at_crds = np.array([i[self.cols[0] != 'Ne'] for i in XYZFile.numeric_data])
        self.cols = np.array([c[self.cols[0] != 'Ne'] for c in self.cols])
        natom = len(at_crds[1])
        
        # Reshape the at_crds array to get mol_crds array
        mol_crds, ats_per_mol = self.__reshape_at_crds(at_crds)
        self.__get_pvec_ats()

        self.data = np.zeros((self.nstep,
                               self.nmol * len(self.pvec_ats[0]),
                               3))
        for step in range(self.nstep):
           at_count = 0
           for imol in range(self.nmol):
               crds = mol_crds[step, imol]
               for iats in self.pvec_ats[imol]:
                   disp1 = crds[iats[0]] - crds[iats[1]]
                   disp2 = crds[iats[2]] - crds[iats[1]]

                   pvec = np.cross(disp1, disp2)
                   pvec /= np.linalg.norm(pvec)

                   self.data[step, at_count] = pvec
                   at_count += 1

    def __get_pvec_ats(self):
        """
        Will get for each carbon atom the 3 closest atoms at step 0.

        This function will determine which atoms to calculate the pvecs with by
        finding the closest atom to each carbon atom at step 0. We only calculate
        them at step 0 to keep the sign of the pvecs consistent throughout the
        steps.

        These atoms are stored in a dictionary -> self.pvec_ats. The first key is
        the mol num and then indexing this list will give the 3 atoms used to
        calculate the pvecs for that atom.
        """
        # Get some initial data
        self.cols = np.reshape(self.cols, (self.nstep, self.nmol, self.at_per_mol))
        self.C_ats = [np.arange(len(c))[c == 'C'] for c in self.cols[0]][0]

        step_NN = self.NN[0]
        at_inds = np.array(step_NN['closest_atoms_mol_grouped'])
        at_inds = at_inds.astype(int)
        self.pvec_ats = {}

        # Loop over all mols and C atoms
        for imol in range(self.nmol):
           self.pvec_ats.setdefault(imol, [])
           for C_at in self.C_ats:
              closest_3_ats = at_inds[imol][C_at, :3]
              self.pvec_ats[imol].append(closest_3_ats)

    def __reshape_at_crds(self, at_crds):
        """
        Will use the self.at_per_mol information to reshape the at_crds array from (self.nstep, natom, 3) to (self.nstep, nmol, self.at_per_mol, 3).

        Inputs:
            * at_crds <np.NDArray> => The atomic coordinates
        Ouputs:
            (<np.NDArray>, <np.NDArray>) The atomic coordinates organised by molecule and the atomic indices for each mol.
        """
        err_msg = "Number of atoms per molecule doesn't neatly divide up the atoms in each xyz step!"

        # Define some consts
        natom = len(at_crds[1])
        self.nmol = natom / self.at_per_mol

        # Error checking for self.nmol
        if type_check.is_int(self.nmol, err_msg):
           self.nmol = int(self.nmol)

        ats_per_mol = np.reshape(np.arange(natom), (self.nmol, self.at_per_mol))
        mol_crds = np.reshape(at_crds, (self.nstep, self.nmol, self.at_per_mol, 3))

        return mol_crds, ats_per_mol

    def __str__(self):
        """
        Overload the string function to display the data in an xyz format
        """
        # Get the atom numbers
        ats = [j for imol in range(self.nmol) for j in self.C_ats + (imol*self.at_per_mol)]

        raise SystemExit("BREAK")
