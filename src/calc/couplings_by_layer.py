#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate couplings from hamiltonian data
"""
import numpy as np
import pandas as pd	

from src.calc import general_calc as gen_calc
from src.calc import molecule_utils as mol_utils


class Layer_Couplings(gen_calc.Calc_Type):
	"""
	Will get the distribution of couplings from hamiltonian data for each predefined layer. 
	"""
	required_data_types = ('pos', 'pseudo_ham',)
	_defaults = {'zmin': 'all', 'zmax': 'all', 'ymin': 'all', 'ymax': 'all', 'xmin': 'all', 'xmax': 'all'}
	required_metadata = ('xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax', "atoms_per_molecule")
	name = "Layer Couplings"

	def calc(self):

		self.ats_per_mol = self.metadata['atoms_per_molecule']

		# Get the xyz data
		self.xyz_data = self.Var.required_data['xyz'].xyz_data
		self.compute_data = self.xyz_data
		self.cols = self.Var.required_data['xyz'].cols
		self.all_mol_crds, self.mol_col = self._get_mol_crds()

		# Get the Hamiltonian data
		self.ham_data = self.Var.required_data['pseudo_ham'].get_data()

		# Reshape the data and 
		all_COMs = mol_utils.get_COM(self.all_mol_crds, self.mol_col)

		mol_nums = self.get_mol_nums_in_layer(all_COMs[0])

		# import matplotlib.pyplot as plt
		# from mpl_toolkits.mplot3d import Axes3D

		# f = plt.figure()
		# a = f.add_subplot(111, projection="3d")

		# # print(self.all_mol_crds[0, mol_nums, :, 0].flatten())
		# # a.plot(self.all_mol_crds[0, mol_nums, :, 0].flatten(), self.all_mol_crds[0, mol_nums, :, 1].flatten(), self.all_mol_crds[0, mol_nums, :, 2].flatten(), 'ko')
		# a.plot(all_COMs[0, mol_nums, 0].flatten(), all_COMs[0, mol_nums, 1].flatten(), all_COMs[0, mol_nums, 2].flatten(), 'ko')
		# # a.plot(self.all_mol_crds[0, :, :, 0].flatten(), self.all_mol_crds[0, :, :, 1].flatten(), self.all_mol_crds[0, :, :, 2].flatten(), 'r.', alpha=0.2)

		# a.set_xlabel("X")
		# a.set_ylabel("Y")
		# a.set_zlabel("Z")

		# plt.show()

	def get_mol_nums_in_layer(self, mol_crds):
		"""
		Will get the molecule numbers of molecules in a layer.
		"""
		xmin, xmax = self.metadata['xmin'], self.metadata['xmax']
		ymin, ymax = self.metadata['ymin'], self.metadata['ymax']
		zmin, zmax = self.metadata['zmin'], self.metadata['zmax']

		nmol = np.shape(mol_crds)[0]
		mask = np.ones(nmol, dtype=bool)
		if xmax != 'all': mask = (mask) & (mol_crds[:, 0] > xmin)
		if xmax != 'all': mask = (mask) & (mol_crds[:, 0] < xmax)
		if ymin != 'all': mask = (mask) & (mol_crds[:, 1] > ymin)
		if ymin != 'all': mask = (mask) & (mol_crds[:, 1] < ymax)
		if zmin != 'all': mask = (mask) & (mol_crds[:, 2] > zmin)
		if zmin != 'all': mask = (mask) & (mol_crds[:, 2] < zmax)

		mol_nums = np.arange(nmol)
		mols_in_layer = mol_nums[mask]

		return mols_in_layer