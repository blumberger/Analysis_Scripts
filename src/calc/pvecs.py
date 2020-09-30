#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate pvecs from some xyz coordinates.
"""
import numpy as np

from src.calc import general_calc as gen_calc
from src.system import type_checking as type_check

class PVecs(gen_calc.Calc_Type):
    """
    Will calculate the pvecs from the data contained within a data file.

    Inputs:
        * Variable <Variable> => A ass or a class that has been derived from Data_File

    Important Attributes:
        * required_metadata <tuple> => Any keys that are required in the metadata dictionary.
        * required_calc <tuple> => Any values that need calculating to calculate this value.
        * data <*> => The data that has been calculated.
    """
    required_metadata = ('atoms_per_molecule',)
    required_calc = ('NN', )
    required_data_types = ('pos',)

    # Need these 3 attribute to create a new variable type
    metadata = {'file_type': 'xyz'}
    name = "P-Vecs"

    def _calc_(self):
        """
        Will calculate the pvecs from self.xyz_data_File.xyz_data
        """
        self.data = []
        all_xyz_data = self.Var.data.get_xyz_data()
        all_cols = self.Var.data.get_xyz_cols()


        # First remove 'Ne' atoms
        all_xyz_data = np.array([
                         [step_xyz[step_cols != 'Ne']
                         for step_cols, step_xyz in zip(file_cols, file_xyz)]
                        for file_cols, file_xyz in zip(all_cols, all_xyz_data)])
        all_cols = np.array([[step_cols[step_cols != 'Ne'] for step_cols in file_cols]
                             for file_cols in all_cols])


        # Loop over the xyz data of each file that is loaded.
        for ifile, (at_crds, self.cols) in enumerate(zip(all_xyz_data, all_cols)):
            self.nstep = len(at_crds)
            self.at_per_mol = self.Var.metadata['atoms_per_molecule']
            self.natom = len(at_crds[0])

            # Reshape the at_crds array to get mol_crds array
            mol_crds, ats_per_mol = self.__reshape_at_crds(at_crds)
            self.__get_pvec_ats(ifile)

            self.pvecs = np.zeros((self.nstep,
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

                       self.pvecs[step, at_count] = pvec
                       at_count += 1


    def __get_pvec_ats(self, ifile):
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

        step_NN = self.NN[ifile][0]
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
        Will use the self.at_per_mol information to reshape the at_crds array from (nstep, natom, 3) to (nstep, nmol, at_per_mol, 3).

        Inputs:
            * at_crds <np.NDArray> => The atomic coordinates
        Ouputs:
            (<np.NDArray>, <np.NDArray>) The atomic coordinates organised by molecule and the atomic indices for each mol.
        """
        err_msg = "Number of atoms per molecule doesn't neatly divide up the atoms in each xyz step!"

        # Define some consts
        self.natom = len(at_crds[1])
        self.nmol = self.natom / self.at_per_mol

        # Error checking for self.nmol
        if type_check.is_int(self.nmol, err_msg):
           self.nmol = int(self.nmol)

        ats_per_mol = np.reshape(np.arange(self.natom), (self.nmol, self.at_per_mol))
        mol_crds = np.reshape(at_crds, (self.nstep, self.nmol, self.at_per_mol, 3))

        return mol_crds, ats_per_mol

    def get_xyz_data(self):
        """Will return the xyz_data that has been created"""
        return self.pvecs

    def get_xyz_cols(self):
        """
        Set the xyz data attributes required for writing an xyz file.

        These are xyz_data, cols and timesteps.
        """
        mols = [str(imol) + "    " for j in self.C_ats for imol in range(self.nmol)]
        ats = [str(j) + "     " for imol in range(self.nmol) for j in self.C_ats + (imol*self.at_per_mol)]
        self.cols = np.char.add(mols, ats)
        return np.array([self.cols] * self.nstep)

    def get_xyz_timesteps(self):
        """Will return 0 for all timesteps."""
        return np.array([0.0] * self.nstep)

    def __str__(self):
        """
        Overload the string function to display the data in an xyz format
        """
        # Get the atom numbers
        mols = [str(imol) + "    " for j in self.C_ats for imol in range(self.nmol)]
        ats = [str(j) + "     " for imol in range(self.nmol) for j in self.C_ats + (imol*self.at_per_mol)]
        cols = np.char.add(mols, ats)

        head_str = f'{len(ats)}\nPvecs. Step:  '
        self.pvecs = self.xyz_data.astype(str)
        
        xyz = (['    '.join(line) for line in step_data] for step_data in self.xyz_data)
        s = (head_str + "%s\n"%step + '\n'.join(np.char.add(cols, step_data)) + "\n"
             for step, step_data in enumerate(xyz))

        return ''.join(s)
