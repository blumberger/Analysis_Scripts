#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate radial distribution functions of mols
"""



import numpy as np
import json
from collections import Counter

from src.data import consts

from src.calc import general_types as gen_type
from src.calc import molecule_utils as mol_utils
from src.calc import geometry as geom


class RDF(gen_type.Calc_Type):
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
                 'max_dist': False, 'number_bins': False}
    required_calc = ()

    # Need these 3 attributes to create a new variable type
    metadata = {'file_type': 'json'}
    name = "Radial Distribution Function"
    with open(consts.PT_FILEPATH) as f: PT = json.load(f)

    def calc(self):
        """
        Will calculate the radial distribution function for the system.
        """
        self.dr = 0.01
        self.cutoff = 12.5
        self.ats_per_mol = self.metadata['atoms_per_molecule']

        self.get_xyz_data()
        self.get_cols(self.metadata['number_each_atom'],
                      self.ats_per_mol)

        xyz = self.compute_data
        mol_xyz = mol_utils.atoms_to_mols(xyz, self.ats_per_mol, nstep=len(xyz))
        nmol = mol_xyz.shape[1]
        mol_col = np.reshape(self.cols[0], (nmol, self.ats_per_mol))
        mask = mol_col == "C"

        x_len = np.max(xyz[:, :, 0]) - np.min(xyz[:, :, 0])
        y_len = np.max(xyz[:, :, 1]) - np.min(xyz[:, :, 1])
        z_len = np.max(xyz[:, :, 2]) - np.min(xyz[:, :, 2])
        self.V = x_len * y_len * z_len

        self.bins = np.arange(0, self.cutoff+self.dr, self.dr)
        self.RDF = np.zeros(len(self.bins))

        # Loop over all available steps
        for step_xyz in mol_xyz:
            just_carbons = step_xyz[mask]
            self.N = len(just_carbons)
            mol_inds, _ = np.mgrid[0:nmol,0:self.ats_per_mol]
            mol_inds = mol_inds[mask]
            self.calc_RDF(just_carbons, mol_inds,
                          rdf_type=self.metadata['rdf_type'])

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

        N = len(pos_data)
        mol_ind = mol_inds[0]
        mask = mol_inds != mol_ind

        # Loop over all atoms and get the RDF contribution from each one
        for at_num, xyz in enumerate(pos_data):

            if mol_inds[at_num] != mol_ind:
                mol_ind = mol_inds[at_num]
                mask = mol_inds != mol_ind

            # Counting all distances which means we are doing double to calculations
            # This can be optimised later.
            dist = np.linalg.norm(pos_data[mask] - xyz, axis=1)
            dist = dist[dist < self.cutoff]


            bin_index = (dist // self.dr).astype(int)
            for i, bin_ind in enumerate(bin_index):
                R1 = dist[i]
                R2 = R1 + self.dr
                vol_sect = geom.volume_concentric_spheres(R1, R2)

                self.RDF[bin_ind] += self.N / (vol_sect * self.V * 2)

        import matplotlib.pyplot as plt
        plt.plot(self.bins, self.RDF)
        plt.show()



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

    def get_num_bins(self):
        """
        Will set the number of bins paramter.
        """
        self.nbins = self.metadata['number_bins']
        if type(self.nbins) != int:
            self.nbins = 70

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
