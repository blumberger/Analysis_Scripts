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
    required_metadata = ()
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
        self.volumes = np.zeros(len(self.radii))
        for i, r in enumerate(self.radii):
            v1 = geom.volume_sphere(r)
            v2 = geom.volume_sphere(r + dr)
            self.volumes[i] = v2 - v1

    def get_num_bins(self):
        """
        Will set the number of bins paramter.
        """
        self.nbins = self.metdata['number_bins']
        if type(num_bins) != int:
            num_bins = 70

    def get_max_dist(self):
        """
        Will get the maximum distance to go up to calculating the RDF
        """
        dist = np.max(self.compute_data[0], axis=0)                           \
               - np.min(self.compute_data[0], axis=0)
        self.max_dist = np.linalg.norm(dist) * 0.6

    def calc(self):
        """
        Will calculate the angular distribution of the molecular system.

        This will loop over each molecule, find the rotation (wrt long ax, short
        ax of central molecule) of the long and short axes of the molecule then
        create a histogram of this data.
        """
        self.get_data()
        all_at_crds = self.compute_data

        self.get_num_bins()
        self.get_max_dist()




    def get_dr(self, max_dist):
        """
        Will determine the spacing between bins.

        Inputs:
            * max_dist <float> => The maximum distance from the center.
        """
        num_bins = self.metadata['number_bins']
        if type(num_bins) != int:
            # 2 * sqrt[ num of atoms ]
            num_bins = 2 * np.sqrt(len(self.compute_data[0]))

        return max_dist / num_bins
