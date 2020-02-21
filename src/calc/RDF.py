#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate radial distribution functions of mols
"""
import numpy as np
import json

from src.data import consts

from src.calc import general_types as gen_type
from src.calc import molecule_utils as mol_utils
from src.calc import geometry as geom


class RDF(gen_type.Calc_Type):
    """
    Will calculate the radial distribution functions of the molecules in a system.

    Inputs:
        * Variable <Variable> => An instance of the Variable class

    Important Attributes:
        * required_metadata <tuple> => Any keys that are required in the metadata dictionary.
        * required_calc <tuple> => Any values that need calculating to calculate this value.
        * data <*> => The data that has been calculated.
    """
    _write_types = ('json', )
    required_metadata = ('atoms_per_molecule', )
    _defaults = {'rdf_type': 'intermolecular',
                 'max_dist': False, 'number_bins': False}
    required_calc = ()

    # Need these 3 attributes to create a new variable type
    data = {'rdf': [], 'r': []}
    metadata = {'file_type': 'json'}
    name = "Radial Distribution Function"
    with open(consts.PT_FILEPATH) as f: PT = json.load(f)

    def get_data(self):
        """
        Will get the data to use from the inputted class.
        """
        if 'xyz_data' in dir(self.Var.data):
            self.compute_data = self.Var.data.xyz_data
        elif 'csv_data' in dir(self.Var.data):
            self.compute_data = self.Var.data.csv_data[['x', 'y', 'z']].to_numpy()
            self.compute_data = np.array([self.compute_data])

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
        COMS = np.sum(COMs, axis=2)/tot_mass
        return COMS

    def calc(self):
        """
        Will calculate the angular distribution of the molecular system.

        This will loop over each molecule, find the rotation (wrt long ax, short
        ax of central molecule) of the long and short axes of the molecule then
        create a histogram of this data.
        """
        ats_per_mol = self.Var.metadata['atoms_per_molecule']
        self.get_data()
        all_at_crds = self.compute_data
        self.get_num_bins()

        self.RDF = np.zeros(self.nbins)
        self.vols = np.zeros(self.nbins)

        # Loop over all steps
        for at_crds in self.compute_data:
            mol_crds = mol_utils.atoms_to_mols(at_crds, ats_per_mol)
            self.COMs = self.get_COMs(mol_crds, self.Var.data.cols[0, :36])

            self.get_max_dist(self.COMs)
            self.tot_volume = geom.volume_sphere(self.max_dist)
            self.calc_shell_volumes()
            self.vols += self.shell_volumes * len(self.COMs)

            # Loop over all mols
            for i, mol1 in enumerate(self.COMs):
                # self.vols += self.shell_volumes

                # Loop over all mol pairs
                all_dist = np.linalg.norm(self.COMs[i:] - mol1, axis=1)
                for dist in all_dist:
                    index = int(dist // self.dr)
                    if 0 < index < self.nbins:
                        self.RDF[index] += 2.0

            # Now normalise
            vol_per_n = self.tot_volume / len(self.COMs)
            for i, value in enumerate(self.RDF):
                self.RDF[i] = value * vol_per_n / self.vols[i]
