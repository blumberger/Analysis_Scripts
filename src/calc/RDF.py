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
                 'max_dist': 1.0, 'number_bins': 400}
    required_calc = ()

    # Need these 3 attributes to create a new variable type
    metadata = {'file_type': 'json'}
    name = "Radial Distribution Function"
    with open(consts.PT_FILEPATH) as f: PT = json.load(f)

    def calc(self):
        """
        Will calculate the radial distribution function for the system.
        """
        self.ats_per_mol = self.metadata['atoms_per_molecule']

        # Grab the relevant data.
        self.get_xyz_data()
        self.get_cols(self.metadata['number_each_atom'],
                      self.ats_per_mol)

        # Reshape the data a bit to sort into different molecules
        xyz = self.compute_data
        mol_xyz = mol_utils.atoms_to_mols(xyz, self.ats_per_mol, nstep=len(xyz))
        nmol = mol_xyz.shape[1]
        mol_col = np.reshape(self.cols[0], (nmol, self.ats_per_mol))
        mask = mol_col == "C"

        # Will get the properties required to calculate the RDF
        self.get_vitals(xyz)

        # Loop over all available steps and calc RDF for each
        for step_xyz in mol_xyz:
            just_carbons = step_xyz[mask]
            self.N = len(just_carbons)   # Num all carbon atoms

            mol_inds, _ = np.mgrid[0:nmol,0:self.ats_per_mol]
            mol_inds = mol_inds[mask]
            self.calc_RDF(just_carbons, mol_inds,
                          rdf_type=self.metadata['rdf_type'])

        self.RDF /= len(mol_xyz)

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

        print(self.radii)
        print("\r                                               ", end="\r")
        import matplotlib.pyplot as plt
        plt.plot(self.radii[:-1], self.RDF)
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
