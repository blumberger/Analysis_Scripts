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
    data = {'rdf': [], 'r': []}
    metadata = {'file_type': 'json'}
    name = "Radial Distribution Function"
    with open(consts.PT_FILEPATH) as f: PT = json.load(f)

    def calc(self):
        """
        Will calculate the radial distribution function for the system.
        """
        self.ats_per_mol = self.Var.metadata['atoms_per_molecule']
        self.at_types = self.Var.metadata['number_each_atom']

        self.get_data()
        all_at_crds = self.compute_data
        self.get_num_bins()

        self.RDF = np.zeros(self.nbins)
        self.vols = np.zeros(self.nbins)

        # Loop over all steps
        for at_crds in self.compute_data:
            mol_crds = mol_utils.atoms_to_mols(at_crds, self.ats_per_mol)
            self.COMs = self.get_COMs(mol_crds, self.cols[0, :36])

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
        rho = len(self.COMs) / self.tot_volume
        norm = rho * self.vols
        self.RDF /= norm

    def get_data(self):
        """
        Will get the data to use from the inputted class.
        """
        if 'xyz_data' in dir(self.Var.data):
            self.compute_data = self.Var.data.xyz_data
        elif 'csv_data' in dir(self.Var.data):
            self.compute_data = self.Var.data.csv_data[['x', 'y', 'z']].to_numpy()
            self.compute_data = np.array([self.compute_data])

        # Only for COMs
        self.get_cols()

    def get_cols_from_CSV(self):
        """
        Will get the element type columns from a csv file.
        """
        df = self.Var.data.csv_data
        # Error check
        if 'type' not in df.columns:
            raise SystemError("\n\n\nI can't calculate the COMs of the mols. I "
                              + "can't find the atoms types in the data."
                              + "I don't what atom types are in the mol")

        # Error check -if there is a mol with equal nums of atom of x and y type.
        num_elm_in_mol = [self.at_types[i] for i in self.at_types]
        if len(set(num_elm_in_mol)) != len(num_elm_in_mol):
            raise SystemError("\n\n\nCan't compute center of masses for RDF."
                            + " I don't know what types the atoms are and can't"
                            + " work them out as the molecule has 2 or more "
                            + "elements with the same number of atoms.\n\n"
                            + " There are no other RDF methods implemented.")

        # Compare how many atoms of each type are in each molecule and how many there should be.
        num_at_types = Counter(df.loc[:self.ats_per_mol-1, 'type'])
        at_types = {i: self.at_types[i] for i in self.at_types}
        cvt_type = {}
        for i in num_at_types:
            for elm in at_types:
                if at_types[elm] == num_at_types[i]:
                    print(f"Atom type '{i}' is element {elm}.")
                    cvt_type[str(i)] = elm
                    break
            else:
                raise SystemError("Could determine what element each atom type is.")
            at_types.pop(elm)

        # Create the cols array
        self.cols = df['type'].to_numpy().astype(str)
        self.natom = len(self.cols)
        for i in np.unique(self.cols):
            self.cols[self.cols == i] = cvt_type[i]
        self.cols = np.array([self.cols])

    def get_cols(self):
        """
        Will try to get the atom element types from the number of each type in a
        molecule.
        """

        if 'cols' not in dir(self.Var.data):
            if 'csv_data' in dir(self.Var.data):
                self.get_cols_from_CSV()
            else:
                raise SystemError("Can't compute center of masses for RDF."
                               + " There are no other RDF methods implemented.")
        else:
            self.cols = self.Var.data.cols

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
