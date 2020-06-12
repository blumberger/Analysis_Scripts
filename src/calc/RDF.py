#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate radial distribution functions of mols
"""
import numpy as np
import json
import pandas as pd
import matplotlib.pyplot as plt

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
    _write_types = ('csv', )
    required_metadata = ('atoms_per_molecule', 'number_each_atom')
    _defaults = {'rdf_type': 'intermolecular', 'plot_RDF': False,
                 'max_dist': 1.0, 'number_bins': 400, 'atom_list_1': ['C'],
                 'atom_list_2': ['C'], 'dr': False, 'cutoff': False}
    required_data_types = ('pos',)

    # Need these 3 attributes to create a new variable type
    metadata = {'file_type': 'csv'}
    name = "Radial Distribution Function"
    with open(consts.PT_FILEPATH) as f: PT = json.load(f)

    def _calc_(self):
        """
        Will calculate the radial distribution function for the system.
        """
        self.ats_per_mol = self.metadata['atoms_per_molecule']
        self.Var['coordinate_wrapping'] = 'wrapped'
        all_xyz_data = self.Var.data.get_xyz_data()
        all_cols = self.Var.data.get_xyz_cols(self.metadata['number_each_atom'])
        nfiles, nstep, natom, ndim = np.shape(all_xyz_data)

        # Set cell vecs in the required format
        ABC = [ [self.Var['xlo'], self.Var['xhi'], self.Var['xy']],
                [self.Var['ylo'], self.Var['yhi'], self.Var['xz']],
                [self.Var['zlo'], self.Var['zhi'], self.Var['yz']] ]

        # Set the atomic coords
        # for ifile in range(len(all_xyz_data)):
        file_data = all_xyz_data[0]
        file_cols = all_cols[0]

        system_info = geom.get_system_size_info(all_xyz_data)
        
        cutoff = self.metadata['cutoff']
        if cutoff is False or cutoff == 'auto':
            cutoff = min([system_info['xlen'], system_info['ylen'], system_info['zlen']]) / 2.

        dr = self.metadata['dr']
        if dr is False or dr == "auto": dr = cutoff / (2 * np.sqrt(nstep * natom))

        self.rdf = []
        for step_xyz, step_cols in zip(file_data, file_cols):
            self.radii, rdf_vals = rdf.calc_RDF(step_xyz, self.ats_per_mol, step_cols,
                                                ABC, self.metadata['atom_list_1'][:],
                                                self.metadata['atom_list_2'][:],
                                                dr=dr, cutoff=cutoff)

            self.rdf.append(rdf_vals)

        self.rdf = np.mean(self.rdf, axis=0)

        if self.metadata['plot_RDF']:
            self._plot_()
            plt.show()

    def _plot_(self, a=False):
        """
        Will make a quick plot of the rdf vs radius
        """
        if a is False:
            f, a = plt.subplots()
        a.plot(self.radii, self.rdf)
        a.set_ylabel("RDF")
        a.set_xlabel(r"R [$\AA$]")


    def get_csv_data(self):
        """
        Will create the string to write.
        """
        return pd.DataFrame({'rdf': self.rdf, 'radii': self.radii})