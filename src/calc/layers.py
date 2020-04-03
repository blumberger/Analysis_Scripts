#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to find any layers in a molecular structure
"""
import numpy as np
import pandas as pd	
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from src.calc import general_calc as gen_calc
from src.calc import molecule_utils as mol_utils
from src.calc import geometry as geom

class Molecular_Layers(gen_calc.Calc_Type):
	"""
	Will get the distribution of couplings from hamiltonian data for each predefined layer. 
	"""
	required_data_types = ('pos', )
	required_calc = ("long_ax_rotation", )
	_defaults = {'allow_boundaries': False, "plot_layers": True}
	required_metadata = ("atoms_per_molecule", )
	name = "Molecular Layers"

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


		for ifile in range(len(xyz_data)):
			for istep in range(len(xyz_data[ifile])):
				step_data = xyz_data[ifile][istep]
				cols = cols[ifile][istep]
				
				mol_crds = mol_utils.atoms_to_mols(step_data, self.metadata['atoms_per_molecule'])
				mol_col = mol_utils.cols_to_mols(cols, self.metadata['atoms_per_molecule'])

				COM = mol_utils.get_COM_split_mols(mol_crds, mol_col)
				rotated_COM = geom.rotate_crds(COM, self.long_ax_rotation.xy_rotation_matrix)
				self.sys_info = geom.get_system_size_info(rotated_COM)

				self.layer_starts = self.get_layers(rotated_COM)
				self.layer_mols = self.get_layer_mols(rotated_COM)

				if self.metadata['plot_layers'] is True:
					self.plot_layers()
					plt.show()

	def get_layer_mols(self, mols):
		"""
		Will get the molecules in each layer.

		Inputs:
			* mols <array> => The molecules in the system.
		Outputs:
			* <list> The molecules arranged into layers
		"""
		layer_mols = []
		for ilayer in range(len(self.layer_starts)-1):
			start, end = self.layer_starts[ilayer], self.layer_starts[ilayer+1]
			mask = (mols[:, 2] < end) & (mols[:, 2] > start)

			layer_mols.append(mols[mask])

		return layer_mols


	def get_layers(self, xyz):
		"""
		Will find the layers in the system.

		Inputs:
			* xyz <arr> => The pos array of shape (N, 3)
		
		Outputs:
			<list> The start position of each layer.
		"""
		# Scan across the system and find the layer structure.
		self.z, self.num = [], []
		for zmin in np.linspace(self.sys_info['zmin'], self.sys_info['zmax'], 500):
			mask = (xyz[:, 2] > zmin) & (xyz[:, 2] < zmin + 5)
			self.z.append(zmin)
			self.num.append(sum(mask))

		# Smooth the data
		self.smoothed_df = pd.DataFrame({'num': self.num, 'z': self.z})
		self.smoothed_df = self.smoothed_df.rolling(10, center=True).mean().dropna()
		# self.smoothed_df['gradient'] = np.abs(self.smoothed_df['num'] - np.roll(self.smoothed_df['num'], 10))
		self.smoothed_df.index = np.arange(len(self.smoothed_df))

		# Get some starting points for the minima
		max_data, min_data = np.max(self.smoothed_df['num']), np.min(self.smoothed_df['num'])
		data_range = max_data - min_data

		# Use a steepest descent-eqsue algorithm to get true local minima
		true_min = []
		for i in np.linspace(self.sys_info['zmin'], self.sys_info['zmax'], 100):
			z_ind = self.smoothed_df.index[np.argmin(np.abs(self.smoothed_df['z'] - i))]

			new_ind, _ = geom.find_local_min(self.smoothed_df['num'], z_ind, self.metadata['allow_boundaries'])

			if new_ind is not False:
				true_min.append(self.smoothed_df.loc[new_ind, 'z'])

		_, true_min = geom.cluster_points(true_min, 0.03)

		return true_min

	def plot_layers(self):
		"""
		"""
		f = plt.figure()
		a2D = f.add_subplot(121)
		a3D = f.add_subplot(122, projection="3d")
		
		a3D.set_xlabel("X"); a3D.set_ylabel("Y"); a3D.set_zlabel("Z");
		a3D.set_xticks([]);  a3D.set_yticks([]);  a3D.set_zticks([]);
		for direction in ["xzero", "yzero"]:
			a3D.axis[direction].set_axisline_style("-|>")
			a3D.axis[direction].set_visible(True)

		a3D.view_init(elev=-3, azim=64)
		for i in self.layer_mols:
			self._plot_xyz_data(i, a3D)

		a2D.plot(self.z, self.num)
		a2D.plot(self.smoothed_df['z'], self.smoothed_df['num'], 'k--', alpha=0.5)

		for i in self.layer_starts:
			a2D.axvline(i, ls='--', color='k', lw=1)
		
		a2D.set_ylabel("Mol Count")
		a2D.set_xlabel(r"Z [$\AA$]")
		
		plt.tight_layout()
		return a2D
		


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