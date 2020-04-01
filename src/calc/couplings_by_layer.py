#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate couplings from hamiltonian data
"""
import numpy as np
import pandas as pd	
# from sklearn.cluster import DBSCAN 
# from sklearn.preprocessing import StandardScaler, normalize 

from src.calc import general_calc as gen_calc
from src.calc import molecule_utils as mol_utils
from src.calc import geometry as geom

from src.data import consts


class Layer_Couplings(gen_calc.Calc_Type):
	"""
	Will get the distribution of couplings from hamiltonian data for each predefined layer. 
	"""
	required_data_types = ('pos', 'pseudo_ham',)
	required_calc = ("long_ax_rotation",)
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
		xyz_data = self.Var.data.get_xyz_data()
		cols = self.Var.data.get_xyz_cols()
		couplings = self.Var.data['pseudo_ham'].data[0]

		for ifile in range(len(xyz_data)):
			for istep in range(len(xyz_data[ifile])):
				step_data = xyz_data[ifile][istep]
				cols = cols[ifile][istep]
				
				mol_crds = mol_utils.atoms_to_mols(step_data, self.metadata['atoms_per_molecule'])
				mol_col = mol_utils.cols_to_mols(cols, self.metadata['atoms_per_molecule'])

				# print(np.std(mol_crds))
				import matplotlib.pyplot as plt
				plt.hist(np.std(mol_crds, axis=1), bins="auto")
				plt.show()
				
				COM = mol_utils.get_COM(mol_crds, mol_col)
				# rotated_COM = geom.rotate_crds(COM, self.long_ax_rotation.xy_rotation_matrix)
				# rotated_crds = geom.rotate_crds(step_data, self.long_ax_rotation.xy_rotation_matrix)

				# self.sys_info = geom.get_system_size_info(rotated_COM)
				# layer_starts = self.get_layers(rotated_COM)
				# count = 0
				# for start, end in zip(layer_starts[:-1], layer_starts[1:]):
				# 	mask = (rotated_COM[:, 2] < end) & (rotated_COM[:, 2] > start)
				# 	data = rotated_COM[mask]
				# 	mol_nums = np.arange(len(rotated_COM))[mask]

				# 	mask = (rotated_crds[:, 2] < end) & (rotated_crds[:, 2] > start)
				# 	at_layer = rotated_crds[mask]

				# 	layer_couplings = []
				# 	x, y, z = [], [], []
				# 	u, v, w, s = [], [], [], []
					# for mol_id, mol_num in enumerate(mol_nums):
					# 	if mol_num not in couplings: continue

					# 	mol_couplings = couplings[mol_num]
					# 	if len(mol_couplings) == 1: continue


					# 	for j in mol_couplings:
					# 		if j == mol_num: continue
							
					# 		layer_couplings.append(mol_couplings[j])


					# 		x.append(data[mol_id, 0])
					# 		y.append(data[mol_id, 1])
					# 		z.append(data[mol_id, 2])

					# 		iu, iv, iw = rotated_COM[j, 0]-data[mol_id, 0], rotated_COM[j, 1]-data[mol_id, 1], rotated_COM[j, 2]-data[mol_id, 2]
					# 		# iu, iv, iw = 1, 1, 1
					# 		u.append(iu)
					# 		v.append(iv)
					# 		w.append(iw)

					# 		s.append(mol_couplings[j]*10)
					# 	break

					# layer_couplings = np.array(layer_couplings) * consts.Ha_to_meV




					# import matplotlib.pyplot as plt
					# from mpl_toolkits.mplot3d import Axes3D

					# f = plt.figure()
					# a = f.add_subplot(111, projection="3d")
					# a.set_xlabel("X")
					# a.set_ylabel("Y")
					# a.set_zlabel("Z")

					# a.view_init(elev=90, azim=64)

					# # a.plot(x, y, z, 'bo')
					# # a.plot(at_layer[:, 0], at_layer[:, 1], at_layer[:, 2], 'k.', alpha=0.2)
					# a.plot(rotated_COM[:, 0], rotated_COM[:, 1], rotated_COM[:, 2], 'k.', alpha=0.2)
					# v = rotated_COM[1] - rotated_COM[0]
					# for i in mol_nums:
					# 	if i not in couplings: continue

					# 	for j in couplings[i]:
					# 		pos1, pos2 = rotated_COM[i], rotated_COM[j]
					# 		v = pos2 - pos1
					# 		a.quiver([pos1[0]], [pos1[1]], [pos1[2]], [v[0]], [v[1]], [v[2]])
					# 	break
					# for ix, iy, iz, iu, iv, iw, il in zip(x, y, z, u, v, w, s):
					# 	a.quiver([ix], [iy], [iz], [iu], [iv], [iw])



					# a.set_xlim([0, 100])
					# a.set_ylim([0, 100])
					# # a.set_zlim([0, 100])
					# plt.tight_layout()
					# # plt.savefig(f"Quiver_{count}.png")
					# plt.show()
					# count += 1
					# break
					# raise SystemExit("BREAK")

					# plt.hist(layer_couplings, bins="auto")
					# plt.show()


	def get_layers(self, xyz):
		"""
		Will find the layers in the system.

		Inputs:
			* xyz <arr> => The pos array of shape (N, 3)
		
		Outputs:
			<list> The start position of each layer.
		"""
		# Scan across the system and find the layer structure.
		z, num = [], []
		for zmin in np.linspace(self.sys_info['zmin'], self.sys_info['zmax'], 500):
			mask = (xyz[:, 2] > zmin) & (xyz[:, 2] < zmin + 5)
			z.append(zmin)
			num.append(sum(mask))

		# Smooth the data
		df = pd.DataFrame({'num': num, 'z': z})
		df = df.rolling(int(len(df['z']) / 50), center=True).mean().dropna()
		df['gradient'] = np.abs(df['num'] - np.roll(df['num'], 10))
		df.index = np.arange(len(df))

		# Get some starting points for the minima
		max_data, min_data = np.max(df['num']), np.min(df['num'])
		data_range = max_data - min_data
		min_grad, max_grad = np.min(df['gradient']), np.max(df['gradient'])
		grad_range = max_grad - min_grad
		minima = df['z'] [ (df['gradient'] < (grad_range * 0.2)) & (df['num'] < min_data + (0.7 * data_range))]

		# Cluster the minima to remove noise
		_, clustered_minima = geom.cluster_points(minima)

		# Use a steepest descent-eqsue algorithm to get true local minima
		true_min = []
		for i in clustered_minima:
		    z_ind = df.index[np.argmin(np.abs(df['z'] - i))]

		    new_ind, _ = geom.find_local_min(df['num'], z_ind)

		    if new_ind is not False:
		        true_min.append(df.loc[new_ind, 'z'])

		_, true_min = geom.cluster_points(true_min, 0.1)

		# import matplotlib.pyplot as plt
		# plt.plot(z, num)
		# # for i in minima:
		# # 	plt.axvline(i, alpha=0.5, color='r')
		# for i in true_min:
		# 	plt.axvline(i, ls='--', color='k', lw=1)
		# plt.ylabel("Mol Count")
		# plt.xlabel(r"Z [$\AA$]")
		# plt.show()

		return true_min


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