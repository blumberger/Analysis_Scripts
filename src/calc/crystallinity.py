#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate the degree of crystallinity of a structure.

The degree of crystallinity is calculated by seeing how far (as a percentage)
a structure's density is from the crystal or amorphous density.
"""
import numpy as np

from src.calc import general_calc as gen_calc
from src.calc import molecule_utils as mol_utils

def calc_all_dists(pos1, all_pos, cutoff=10):
	dist = np.linalg.norm(all_pos - pos1, axis=1)
	dist[np.abs(dist) > cutoff] = 0.0
	return dist

class Crystallinity(gen_calc.Calc_Type):
	"""
	Will calculate the Density from the data contained within a data file.

	Inputs:
		* Variable <Variable> => An instance of the Variable class

	Important Attributes:
		* required_metadata <tuple> => Any keys that are required in the metadata dictionary.
		* required_calc <tuple> => Any values that need calculating to calculate this value.
		* data <*> => The data that has been calculated.
	"""
	required_metadata = ("crystal_density", "amorphous_density", "atoms_per_molecule",
						 "number_each_atom")
	required_calc = ('density', )
	required_data_types = ('pos',)

	# Need these 3 attributes to create a new variable type
	metadata = {'file_type': 'csv'}
	name = "Crystallinity"

	def calc(self):
		self.global_density = [i[-1] for i in self.density][-1]

		self.Var['coordinate_wrapping'] = 'wrapped'
		self.get_xyz_data()
		self.get_cols(self.metadata['number_each_atom'], self.metadata['atoms_per_molecule'])

		self.C_dens = self.metadata['crystal_density']
		self.A_dens = self.metadata['amorphous_density']
		self.dens_diff = abs(self.C_dens - self.A_dens)
		
		self.global_crysallinity = self.get_degree_of_crystallinity(self.global_density)
		print(f"Global Density = {self.global_density:.2f} g cm^{-3}")
		print(f"Global Crystallinity = {self.global_crysallinity:.1f}%")

		# # Will create some reshaped arrays to more conviently loop over molecules. These are:
		# #  * self.all_mol_crds -> The molecular coords in shape (nstep, nmol, nat_per_mol, 3)
		# #  * self.mol_col      -> The atom types in shape (nmol, nat_per_mol)
		# self._get_mol_crds()

		# # Get center of masses
		# self.all_COM = mol_utils.get_COM(self.all_mol_crds, self.mol_col)

		# Loop over center of masses
		for at_crds in self.compute_data:

			self._calc_local_crystallinities(at_crds, self.cols)
			break









	def _calc_local_crystallinities(self, crds, at_types, grid_points=False):
		"""
		Will calculate the local crystalinity at grid points
		"""
		xspaces, yspaces, zspaces = 40, 100, 40
		if grid_points is False:
			self._create_grid_points(xspaces, yspaces, zspaces)
		else:
			raise SystemExit("Can't handle grid_points currently")


		dx = (self.xmax - self.xmin) / xspaces
		dy = (self.ymax - self.ymin) / yspaces
		dz = (self.zmax - self.zmin) / zspaces

		vol = dx * dy * dz * 1e-24  # convert to cm
		masses = mol_utils.get_atom_masses(at_types[0]) * 1.6602752976190478e-24  # convert to grams

		grid_points = np.zeros((xspaces+1, yspaces+1, zspaces+1))
		plt_crds, plt_mass = [], []
		for crd, mass in zip(crds, masses):
			x, y, z = crd

			ix = int( (x - self.xmin) // dx )
			iy = int( (y - self.ymin) // dy )
			iz = int( (z - self.zmin) // dz )

			grid_points[ix, iy, iz] += mass

		densities = grid_points / vol
		self.local_crystallinities = self.get_degree_of_crystallinity(densities)
		# min_crys = np.min(self.local_crystallinities)
		# self.local_crystallinities -= min_crys
		# max_crys = np.max(self.local_crystallinities)
		# self.local_crystallinities /= max_crys
		# self.local_crystallinities *= 100

		import matplotlib.pyplot as plt

		for i in range(yspaces):
			print(f"\r{i}", end="\r")
			plt.imshow(self.local_crystallinities[:, i, :],
					   extent=(self.xmin+dx/2, self.xmax-dx/2, self.zmin+dz/2, self.zmax-dz/2),
				       cmap="Oranges", interpolation="spline36", vmin=0, vmax=100)
			
			plt.colorbar()
			plt.savefig(f"frames/frame_{str(i).zfill(5)}.png")
			plt.close()

	def _create_grid_points(self, xspaces=30, yspaces=30, zspaces=30):
		"""
		Will create the grid and at every point the local crystallinity will be calculated.
		"""
		x = self.compute_data[0][:, 0]
		y = self.compute_data[0][:, 1]
		z = self.compute_data[0][:, 2]

		self.xmin, self.xmax = np.min(x), np.max(x)
		self.ymin, self.ymax = np.min(y), np.max(y)
		self.zmin, self.zmax = np.min(z), np.max(z)

		self.xlen = self.xmax - self.xmin
		self.ylen = self.ymax - self.ymin
		self.zlen = self.zmax - self.zmin

		self.x_points = np.linspace(self.xmin, self.xmax, xspaces)
		self.y_points = np.linspace(self.ymin, self.ymax, yspaces)
		self.z_points = np.linspace(self.zmin, self.zmax, zspaces)



	def get_degree_of_crystallinity(self, density):
		"""
		Will output the degree of crystallinity.

		rho_crystal = 100
		rho_amorphous = 0

		crystallinity = 100 * ( (rho - rho_amorph) / |rho_crys - rho_amorph| )

		Inputs:
			* density <float> => The density of the sample (rho in the formula)
		Outputs:
			<float> The degree of crysallinity between 0 and 100.
		"""
		return 100 * ( (density - self.A_dens) / (self.dens_diff) )

