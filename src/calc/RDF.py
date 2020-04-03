#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate radial distribution functions of mols
"""



import numpy as np
import json
from collections import Counter

# Own C modules
from src.wrappers import RDF_wrap as rdf

# Own Python Modules
from src.data import consts

from src.calc import general_calc as gen_calc
from src.calc import molecule_utils as mol_utils
from src.calc import geometry as geom



class RDF(gen_calc.Calc_Type):
    """
    Will calculate the radial distribution functions of the molecules in a system.

    The calc function is the function that is called to calculate the RDF, see
    the Calc_Type.

    Inputs:
        * Variable <Variable> => An instance of the Variable class

    Important Attributes:
        * required_metadata <tuple> => Any keys that are required in the metadata dictionary.
        * required_calc <tuple> => Any values that need calculating to calculate this value.
        * data <*> => The data that has been calculated.
    """
    _write_types = ('json', )
    required_metadata = ('atoms_per_molecule', 'number_each_atom')
    _defaults = {'rdf_type': 'intermolecular',
                 'max_dist': 1.0, 'number_bins': 400}
    required_data_types = ('lammps_dump',)

    # Need these 3 attributes to create a new variable type
    metadata = {'file_type': 'json'}
    name = "Radial Distribution Function"
    with open(consts.PT_FILEPATH) as f: PT = json.load(f)

    def _calc_(self):
        """
        Will calculate the radial distribution function for the system.
        """
        self.ats_per_mol = self.metadata['atoms_per_molecule']
        self.Var['coordinate_wrapping'] = 'wrapped'
        self.Var.data.get_xyz_data()
        self.Var.data.get_xyz_cols(self.metadata['number_each_atom'],
                      self.ats_per_mol)

        # Set cell vecs in the required format
        ABC = [ [self.Var['xlo'], self.Var['xhi'], self.Var['xy']],
                [self.Var['ylo'], self.Var['yhi'], self.Var['xz']],
                [self.Var['zlo'], self.Var['zhi'], self.Var['yz']] ]

        # Set the atomic coords
        at_crds = self.compute_data[0]

        # Get the atomic types
        at_types = self.cols[0]

        
        self.radii, self.rdf = rdf.calc_RDF(self.compute_data[0], self.ats_per_mol, self.cols[0],
                                            ABC, ['C'], ['C'], dr=0.01, cutoff=12.5)

        raise SystemExit("RDF calculator is still in progress!")

    def get_vitals(self, xyz):
        """
        Will get some vital properties for calculating the RDF.

        The properties that are calculated are:
            * self.x_len <int> => The length in x
            * self.y_len <int> => The length in y
            * self.z_len <int> => The length in z
            * self.cutoff <float> => The maximum dist we use for calculating RDF
            * self.dr <float> => The spacing between bins
            * self.V <float> => The volume of the full system.
            * self.radii <array> => The bin edges
            * self.RDF <array> => The RDF array.

        Inputs:
            * xyz <array> => The position array of shape (nstep, natom, 3)
        """
        self.x_len = np.max(xyz[:, :, 0]) - np.min(xyz[:, :, 0])
        self.y_len = np.max(xyz[:, :, 1]) - np.min(xyz[:, :, 1])
        self.z_len = np.max(xyz[:, :, 2]) - np.min(xyz[:, :, 2])
        if type(self.metadata['max_dist']) == float:
            if self.metadata['max_dist'] <= 1:
                self.cutoff = min([self.x_len, self.y_len, self.z_len])
                self.cutoff *= self.metadata['max_dist']
            else:
                self.cutoff = self.metadata['max_dist']
        elif type(self.metadata['max_dist']) == int:
            if self.metadata['max_dist'] == 1:
                self.cutoff = min([self.x_len, self.y_len, self.z_len])
            self.cutoff = self.metadata['max_dist']
        else:
            raise SystemExit("I don't know how to handle the 'max_dist' parameter."
                             + " It should be an int or a float.")

        self.dr = self.cutoff / self.metadata['number_bins']

        # Assume a cube.
        self.V = self.x_len * self.y_len * self.z_len

        self.radii = np.arange(0, self.metadata['number_bins']+1, dtype=float)
        self.radii += 0.5
        self.radii *= self.dr

        self.RDF = np.zeros(len(self.radii))
        self.vols = (self.radii[:-1]+(self.dr/2))**2 * 4 * np.pi * self.dr


    def calc_RDF(self, pos_data, mol_inds, rdf_type="intermolecular"):
        """
        Will calculate RDF for the provided position data.

        Inputs:
            * pos_data <array> => The data to calculate the RDF from
            * mol_inds <array> => Same length as pos_data giving the molecular index
        """
        if rdf_type == "intermolecular":
            make_mask = lambda mol_ind: mol_inds != mol_ind
        elif rdf_type == "intramolecular":
            make_mask = lambda mol_ind: mol_inds == mol_ind
        else:
            raise SystemError("\n\n\nI don't understand the type of RDF you want."
                 + "\n\nChoose from:\n\t* 'intermolecular'\n\t* 'intramolecular'")

        mol_ind = mol_inds[0]
        mask = mol_inds != mol_ind
        pos = pos_data[mask]

        # Loop over all atoms and get the RDF contribution from each one
        self.RDF = np.zeros(self.metadata['number_bins'])
        self.rho = self.V / self.N
        self.norm = 1/(self.vols * self.rho)
        for at_num, xyz in enumerate(pos_data):
            print(f"\r{at_num}/{len(pos_data)}", end="\r")
            if mol_inds[at_num] != mol_ind:
                mol_ind = mol_inds[at_num]
                mask = mol_inds != mol_ind
                pos = pos_data[:at_num][mask[:at_num]]

            # We only need to do up to the atom number as the NN matrix is symmetric
            dist = np.linalg.norm(pos - xyz, axis=1)
            dist = dist[dist < self.cutoff]

            # Get the histogram of the data
            C, bin_edges = np.histogram(dist, bins=self.radii)
            self.RDF += C *self.norm

        #print(self.radii)
        #print("\r                                               ", end="\r")
        #import matplotlib.pyplot as plt
        #plt.plot(self.radii[:-1], self.RDF)
        #plt.show()



    def calc_shell_volumes(self):
        """
        Will calculate the volumes of all the spherical shells.
        """
        self.dr = self.max_dist / self.nbins
        self.radii = np.linspace(0, self.nbins * self.dr, self.nbins)
        self.shell_volumes = np.zeros(len(self.radii))
        for i, r in enumerate(self.radii):
            v1 = geom.volume_sphere(r)
            v2 = geom.volume_sphere(r + self.dr)
            self.shell_volumes[i] = v2 - v1

    def get_max_dist(self, pos):
        """
        Will get the maximum distance to go up to calculating the RDF
        """
        dist = np.max(pos, axis=0) - np.min(pos, axis=0)
        self.max_dist = np.linalg.norm(dist)

    def get_COMs(self, mol_crds, elm_names):
        """
        Will get the center of masses of each molecule in an array.

        Inputs:
            * mol_crds <np.NDArray> => Array must be of shape:
                                        (nstep, nmol, nat_per_mol, 3)
            * elm_name <list|array> => Elemental symbol for each atom on 1
                                       molecule. Must be of shape (nat_per_mol)
        """
        # First get the masses from the periodic table
        unique_elm_names = np.unique(elm_names)
        masses = {n: self.PT[i]['atomic_weight'] for n in unique_elm_names
                      for i in self.PT if self.PT[i]['abbreviation'] == n}
        for name in masses:
            elm_names[elm_names == name] = masses[name]
        masses = elm_names.astype(float)

        # Now multiply coords by masses
        tot_mass = sum(masses)
        COMs = [[crds[:, 0] * masses, crds[:, 1] * masses, crds[:, 2] * masses]
                for crds in mol_crds]
        COMs = np.sum(COMs, axis=2) / tot_mass

        return COMs
