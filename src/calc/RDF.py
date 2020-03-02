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
        self.data = {'rdf': [], 'r': []}
        self.ats_per_mol = self.Var.metadata['atoms_per_molecule']
        self.at_types = self.Var.metadata['number_each_atom']

        self.get_data()
        self.get_cols()
        all_at_crds = self.compute_data
        print(np.shape(self.compute_data))
        self.get_num_bins()

        self.RDF = np.zeros(self.nbins)
        self.vols = np.zeros(self.nbins)

        # Loop over all steps
        for at_crds in self.compute_data:
            mol_crds = mol_utils.atoms_to_mols(at_crds, self.ats_per_mol)
            self.COMs = self.get_COMs(mol_crds, self.cols[0, :36])

            N = len(self.COMs)
            # print(self.COMs[:, 0] - 39.788086)

            self.get_max_dist(self.COMs)
            self.tot_volume = geom.volume_sphere(self.max_dist)
            self.calc_shell_volumes()
            self.vols += self.shell_volumes * len(self.COMs)

            # Loop over all mols
            for i, mol1 in enumerate(self.COMs):
                # self.vols += self.shell_volumes

                # Loop over all mol pairs
                all_dist = np.linalg.norm(self.COMs[i+1:] - mol1, axis=1)
                for dist in all_dist:
                    index = int(dist // self.dr)
                    if 0 < index < self.nbins:
                        self.RDF[index] += 2.0


        # Now normalise
        rho = (N*(N-1))  / self.tot_volume
        norm = rho * self.vols
        self.RDF /= norm

        # import matplotlib.pyplot as plt
        # plt.plot(self.radii, self.RDF)
        # plt.xlim([0, 15])
        # plt.show()

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

        # import matplotlib.pyplot as plt
        # from mpl_toolkits.mplot3d import Axes3D
        #
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection="3d")
        # ax.plot(COMs[:, 0], COMs[:, 1], COMs[:, 2], 'ro')
        #
        # # Catoms = self.compute_data[0][self.cols[0] == "C"]
        # # Hatoms = self.compute_data[0][self.cols[0] == "H"]
        # # ax.plot(Catoms[:, 0], Catoms[:, 1], Catoms[:, 2], 'k.')
        # # ax.plot(Hatoms[:, 0], Hatoms[:, 1], Hatoms[:, 2], 'y.')
        #
        # ax.set_xlim([0, 80])
        # ax.set_ylim([0, 80])
        # ax.set_zlim([0, 80])
        #
        # print(self.Var.data.filepath)
        # plt.show()
        # raise SystemExit("RBREAK")


        # s = f"{len(COMs)}"+"\n"+"time: 0.0 fs"+"\n"
        # for i in COMs:
        #     s += f"C {i[0]} {i[1]} {i[2]}" + "\n"
        # print(s)
        # print(COMs)
        #
        # fp = self.Var.data.filepath
        # if '50ps' in fp:
        #     with open("COMs_50ps.xyz", "w") as f:
        #         f.write(s)
        # elif 'crystal' in fp:
        #     with open("COMs_crystal.xyz", "w") as f:
        #         f.write(s)
        # else:
        #     with open("COMs_1ns.xyz", "w") as f:
        #         f.write(s)

        return COMs
