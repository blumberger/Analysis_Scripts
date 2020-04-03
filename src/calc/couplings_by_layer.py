#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate couplings from hamiltonian data
"""
import numpy as np
import pandas as pd	

from src.calc import general_calc as gen_calc
from src.calc import molecule_utils as mol_utils
from src.calc import geometry as geom

from src.io_utils import CP2K_inp_files as cp2k_inp

from src.data import consts


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Layer_Couplings(gen_calc.Calc_Type):
	"""
	Will get the distribution of couplings from hamiltonian data for each predefined layer. 
	"""
	required_data_types = ('pos', 'pseudo_ham',)
	required_calc = ("long_ax_rotation", "mol_layers")
	_defaults = {'zmin': 'all', 'zmax': 'all', 'ymin': 'all', 'ymax': 'all', 'xmin': 'all', 'xmax': 'all'}
	required_metadata = ('xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax', "atoms_per_molecule")
	name = "Layer Couplings"

	def _calc_(self):
		"""
		Will calculate the couplings in each layer.

		This works by first 
		"""
		self.ats_per_mol = self.metadata['atoms_per_molecule']

		# Get the xyz data
		self.Var['coordinate_wrapping'] = 'wrapped'
		xyz_data = self.Var.data.get_xyz_data()
		cols = self.Var.data.get_xyz_cols()
		couplings = self.Var.data['pseudo_ham'].data[0]

		for ifile in range(len(xyz_data)):
			for istep in range(len(xyz_data[ifile])):
				step_data = xyz_data[ifile][istep]
				cols = cols[ifile][istep]
				
				mol_crds = mol_utils.atoms_to_mols(step_data, self.metadata['atoms_per_molecule'])
				mol_col = mol_utils.cols_to_mols(cols, self.metadata['atoms_per_molecule'])

				COM = mol_utils.get_COM_split_mols(mol_crds, mol_col)
				rotated_COM = geom.rotate_crds(COM, self.long_ax_rotation.xy_rotation_matrix)

				rotated_mol_crds = geom.rotate_crds(mol_crds, self.long_ax_rotation.xy_rotation_matrix)

				

				for ilayer in range(len(layer_starts)-1):
					start, end = layer_starts[ilayer], layer_starts[ilayer+1]
					mask = (rotated_COM[:, 2] < end) & (rotated_COM[:, 2] > start)

					at_data = rotated_mol_crds[mask]

					data = rotated_COM[mask]
					mol_nums = np.arange(len(rotated_COM))[mask]
					
					coupl = {mol_nums[i-1]: {mol_nums[j-1]: couplings[i][j] for j in couplings[i]}
							 for i in couplings if type(i) == int}

					# x, y, z, u, v, w, c = [], [], [], [], [], [], []
					# for imol1 in mol_nums:
					# 	if len(coupl[imol1]) <= 1: continue

					# 	for j in coupl[imol1]:
					# 		if j == imol1: continue

					# 		x.append(rotated_COM[imol1, 0])		
					# 		y.append(rotated_COM[imol1, 1])		
					# 		z.append(rotated_COM[imol1, 2])

					# 		u.append(rotated_COM[j, 0] - rotated_COM[imol1, 0])
					# 		v.append(rotated_COM[j, 1] - rotated_COM[imol1, 1])
					# 		w.append(rotated_COM[j, 2] - rotated_COM[imol1, 2])
					# 		c.append(np.abs(coupl[imol1][j])*50)

					# a = self._plot_xyz_data(rotated_COM, alpha=0.5)
					# self._plot_xyz_data(data, a, fmt='ro')
					a = self._plot_xyz_data(at_data[0], fmt="k.", alpha=0.2)
					for i in at_data[1:]:
						self._plot_xyz_data(i, a, fmt="b.", alpha=0.2)
					for ix, iy, iz, iu, iv, iw, il in zip(x, y, z, u, v, w, c):
						a.quiver(ix, iy, iz, iu, iv, iw, length=il, lw=3, color='k')
					
					plt.tight_layout()

					plt.show()

					# plt.savefig(f"Layer_{ilayer}.png")
					# plt.close()
					break








