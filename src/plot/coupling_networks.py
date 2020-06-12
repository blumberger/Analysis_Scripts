#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to plot the coupling networks from coupling data.
"""
# Plotting tools
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

import numpy as np

# Own Modules
from src.plot import general_plot as gen_plot
from src.calc import molecule_utils as mol_utils


class Coupling_Networks(gen_plot.Plot_Type):
	"""
	Will plot the coupling networks.
	"""
	metadata = {}
	_defaults = {
				 "save_coupling_connection_plot": False,
				 "plot_coupling_annotations": "auto", "plot_coupling_atoms": "auto",
				 "plot_coupling_mol_numbers": False,  "a1_elev": 90,
				 "a2_elev": 0, "a1_azim": 0, "a2_azim": 0, "elev_increment":False,
				 "azim_increment":1, "init_elev":30, "init_azim":0,
				 "do_coupling_rotate": False, 'plot_coupling_grains': False,
				 }
	required_metadata = ("reorganisation_energy",)
	required_var_attributes = ("get_xyz_data", "get_xyz_cols", "mol_centroids",
							   "data", "all_centroids", "spliced_crds",
							   "spliced_cols", "mol_nums",
							  )
	name = "Coupling Networks"

	def __init__(self, Variable):
		super().__init__(Variable)

	def _plot_(self):
		"""
		The plotting function which gets called from the generic plot class
		"""
		xyz_data = self.Var.get_xyz_data()
		cols = self.Var.get_xyz_cols()
		for istep in range(len(xyz_data)):
			self.f = plt.figure(figsize=(16, 9))
			self.a1 = self.f.add_subplot(121, projection="3d", proj_type = 'ortho');
			self.a2 = self.f.add_subplot(122, projection="3d", proj_type = 'ortho')

			self.min_Hab = self.metadata["reorganisation_energy"] / 100.


			reorg = self.metadata['reorganisation_energy']
			self.plot_params = {
						   reorg/2.: ({'color': 'r', 'lw': 3}, r"$H_{ab} \geq \frac{\lambda}{2}$"),
						   reorg/10.: ({'color': 'g', 'lw': 1.2}, r"$\frac{\lambda}{2} > H_{ab} \geq \frac{\lambda}{10}$"),
						   self.min_Hab: ({'color': 'b', 'lw': 0.3}, r"$\frac{\lambda}{10} > H_{ab} \geq \frac{\lambda}{100}$"),
						  }
			self.graph_data = self.__get_coupling_connections__(self.Var.mol_centroids[istep], self.Var.data[istep], self.plot_params)
			if not self.metadata['plot_coupling_grains']:
				self._plot_coupling_connections_(self.Var.mol_centroids[istep], self.plot_params, self.a1)
			else: 
				self._plot_coupling_connections_(self.Var.mol_centroids[istep], self.plot_params,
												 self.a1, ignore_cats=[0, 2])
				self.plot_grains(self.Var.mol_centroids[istep], self.a1)

			self._plot_mol_selection_(self.Var.all_centroids[istep], self.Var.mol_centroids[istep], self.a2)

			if self.metadata['plot_coupling_atoms']:
				self._plot_coupling_atoms_(self.Var.spliced_crds[istep], self.Var.spliced_cols[istep], self.a1)
			if self.metadata['plot_coupling_mol_numbers']:
				self._plot_mol_nums_(self.Var.mol_centroids[istep], istep, self.a1)

			self.a1.view_init(elev=self.metadata['a1_elev'], azim=self.metadata['a1_azim'])
			self.a2.view_init(elev=self.metadata['a2_elev'], azim=self.metadata['a2_azim'])

			# Rotate the plot
			if self.metadata['do_coupling_rotate']:	self.rotate_plot()

			# Decide how to handle the plotted data
			plt.tight_layout()
			if bool(self.metadata['save_coupling_connection_plot']):
				plt.savefig(self.metadata['save_coupling_connection_plot'])
				plt.close()
			else:
				plt.show()


	def __get_coupling_connections__(self, all_mol_crds, couplings, plot_params, convert=False):
		"""
		Will return dict with the coupling connection info for plotting.

		Inputs:
			* all_mol_crds <array> => The molecular coordinates
			* couplings <array> => The coupling data.
			* convert <bool> => Convert from Hartree to meV or not. (default no conversion)
		"""
		if convert is False:  conversion = 1
		else:  conversion = consts.Ha_to_meV

		# Loop over coupling mol nums (these are the integer indices in the coupling dict.)
		graph_data = []
		for mol1 in couplings:
			if not isinstance(mol1, (int, np.int64)): continue
			mol_couplings = couplings[mol1]

			max_cat = max(plot_params.keys())

			# Loop over the mols this mol is coupled with.
			for mol2 in mol_couplings:
				if mol1 != mol2:
					Hab = abs(mol_couplings[mol2] * conversion)
					if Hab < self.min_Hab: continue

					for icat, max_coup in enumerate(plot_params):
						if Hab >= max_coup:
							plot_args = plot_params[max_coup][0]
							break
					else:
						# Coupling too low so skip the mol pair
						continue

					pos1, pos2 = all_mol_crds[mol1], all_mol_crds[mol2]
					point_data = {'pos': ((pos1[0], pos2[0]),
										  (pos1[1], pos2[1]),
										  (pos1[2], pos2[2])),
								  'args': plot_args,
								  "category": icat,
								  "has_strong_coupl": abs(Hab) >= max_cat,
								  'Hab': mol_couplings[mol2],
								  "mol_nums": (self.Var.mol_nums[0][mol1],
											   self.Var.mol_nums[0][mol2],),
								  "mol_inds": (mol1, mol2),
								 }
					graph_data.append(point_data)

		return graph_data

	def rotate_plot(self):
		"""
		Will do the rotation of the coupling network axis and savefigs for stitching.
		"""
		def func(count):  plt.savefig(f"frames/rot_{str(count).zfill(5)}.png")
		self._rotate_plot_(self.a1, self.metadata['elev_increment'], self.metadata['azim_increment'],
							   self.metadata['init_elev'], self.metadata['init_azim'],
						   step_func=func)

	def _plot_coupling_atoms_(self, xyz, cols, ax=False):
		"""
		Will plot the atoms on the same graph as the coupling connections

		Inputs:
			* xyz <arr> => The xyz data of the atoms in shape (nmol, nat_per_mol, 3)
			* cols <arr> => The element symbols of the atoms (nmol, nat_per_mol)
			* ax <plt.ax> OPTIONAL => The plot axis

		Outputs:
			<plt.axis> The axis the data is plotted on.
		"""
		plot_atoms = True
		if self.metadata['plot_coupling_atoms'] == "auto":
			plot_atoms = len(xyz) < 7
		elif type(self.metadata['plot_coupling_atoms']) == bool:
			plot_atoms = self.metadata['plot_coupling_atoms']
		
		if plot_atoms is False:
			return

		if ax is False:
			fig = plt.figure()
			ax = fig.add_subplot(111, projection="3d")

		xyz = np.reshape(xyz, (xyz.shape[0] * xyz.shape[1], 3))
		cols = np.reshape(cols, (cols.shape[0] * cols.shape[1]))

		# colors, sizes = cols.copy(), cols.copy()
		for elm in set(cols):
			mask = cols == elm
			color = mol_utils.PT_abbrv[elm]['plot_color']
			size = 6 * (mol_utils.PT_abbrv[elm]['atomic_weight']) ** 0.2
			crds = xyz[mask]
			ax.plot(crds[:,0], crds[:,1], crds[:,2], '.',
					ls="none", ms=size, color=color, alpha=0.5)

		# Reshape axes
		xlim = ax.get_xlim(); ylim = ax.get_ylim(); zlim = ax.get_zlim()
		xdiff = np.diff(xlim); ydiff = np.diff(ylim); zdiff = np.diff(zlim)
		max_diff = max([xdiff, ydiff, zdiff])
		x_ext = abs(max_diff - xdiff) / 2.
		y_ext = abs(max_diff - ydiff) / 2.
		z_ext = abs(max_diff - zdiff) / 2.
		ax.set_xlim([xlim[0] - x_ext, xlim[1] + x_ext])
		ax.set_ylim([ylim[0] - y_ext, ylim[1] + y_ext])
		ax.set_zlim([zlim[0] - z_ext, zlim[1] + z_ext])

		return ax

	def _plot_mol_nums_(self, all_mol_centre, istep=0, ax=False):
		"""
		Will plot the molecule  the graph

		Inputs:
			* all_mol_centre <array> => The coordinates of the molecular centers in shape (nmol, 3)
			* ax <plt.axis> OPTIONAL => The axis to plot on
		"""
		if ax is False:
			fig = plt.figure()
			ax = fig.add_subplot(111, projection="3d")

		for i, xyz in enumerate(all_mol_centre):
			ax.text(xyz[0], xyz[1], xyz[2], r"$\mathbf{%i}$" % self.Var.mol_nums[istep][i])

	def _plot_coupling_connections_(self, all_mol_crds, plot_params, ax=False,
									ignore_cats=()):
		"""
		Will plot the connections via coupling for the given molecular system.

		Inputs:
			* all_mol_crds <array> => The molecular coordinates
			* plot_params <dict> => What categories are used in the plotting.
			* ax <plt.axis> OPTIONAL => The axis to plot on
			* plot_large <bool> => Whether to plot the largest couplings or not.
		"""
		if ax is False:
			fig = plt.figure()
			ax = fig.add_subplot(111, projection="3d")

		ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z");
		ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([]);
		reorg = self.metadata['reorganisation_energy']

		self._plot_xyz_data(all_mol_crds, ax,
							args={'ls': 'none', 'marker': '.', 'color': 'k', 'ms': 3})
		for plot_data in self.graph_data:
			if plot_data['category'] in ignore_cats:
				continue

			ax.plot(*plot_data['pos'], **plot_data['args'])


		annotate = True
		if self.metadata['plot_coupling_annotations'] == 'auto':
			annotate =  len(all_mol_crds) < 42
		elif type(self.metadata['plot_coupling_annotations']) == bool:
			annotate = self.metadata['plot_coupling_annotations']

		if annotate:
			for plot_data in self.graph_data:
				if abs(plot_data['Hab']) >= reorg/10:
					pos = np.mean(plot_data['pos'], axis=1)
					ax.text(pos[0], pos[1], pos[2],	f"{plot_data['Hab']:.1f}",
							horizontalalignment='center', verticalalignment='center')

		ax.set_title(f"$\lambda = $ %.1f meV" % reorg)

		# if self.metadata['CC_plot_title']: ax.set_title(self.metadata['CC_plot_title'].replace("Layer", "").replace("_", " ").strip())

		legend_elements = [Line2D([0], [0], label=plot_params[i][1], **plot_params[i][0]) for i in plot_params]
		ax.legend(handles=legend_elements, loc="best")

	def _plot_mol_selection_(self, all_mols, mol_selection, ax=False):
		"""
		Will plot the molecules selected within the full system.

		Inputs:
			* all_mols <array> => The full system
			* mol_selection <array> => The selected molecules
			* ax <plt.axis> OPTIONAL => The axis to plot on
		"""
		if ax is False:
			fig = plt.figure()
			ax = fig.add_subplot(111, projection="3d")

		ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z");
		ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([]);

		# Plot the full system of molecules as black dots
		self._plot_xyz_data(all_mols, ax,
							args={'ls': 'none', 'marker': '.', 'color': 'k'})

		# Plot the selected molecules as red blobs
		self._plot_xyz_data(mol_selection, ax,
							args={'ls': 'none', 'marker': 'o', 'color': 'r'})

		ax.set_title("Mols selected in red")
		return ax

	def _get_grains_(self):
		"""
		Will get the groups of mols that have the highest category of coupling
		"""
		# Define some vars
		self.grains = {}

		# First group up the coupled mols in a more easy to use way
		strongly_coupled_mols = list(set([i['mol_inds'] for i in self.graph_data if i['has_strong_coupl']]))
		couplings_by_mol = {}
		for mol1, mol2 in strongly_coupled_mols:
			# Check if mol1 has been categorised
			if mol1 not in couplings_by_mol:
				couplings_by_mol.setdefault(mol1, []).append(mol2)
			else:
				if mol2 not in couplings_by_mol[mol1]:
					couplings_by_mol[mol1].append(mol2)

			# Check if mol2 has been categorised
			if mol2 not in couplings_by_mol:
				couplings_by_mol.setdefault(mol2, []).append(mol1)
			else:
				if mol1 not in couplings_by_mol[mol2]:
					couplings_by_mol[mol2].append(mol1)		

		# Now group them by connections
		mols_to_do = list(couplings_by_mol.keys())
		grain_count = 0
		mol = mols_to_do[0]

		# Loop over the maximum number of grains possible
		for i in range(len(mols_to_do) // 2):
			self.grains[grain_count] = self._follow_dict_chain_(couplings_by_mol, mol)
			for imol in self.grains[grain_count]:
				mols_to_do.remove(imol)

			if len(mols_to_do) == 0:
				break

			mol = mols_to_do[0]
			grain_count += 1

	def plot_grains(self, all_mol_crds, ax=False):
		"""
		Will plot the groups of mols with the highest category of coupling.
		"""
		if ax is False:
			fig = plt.figure()
			ax = fig.add_subplot(111, projection="3d")

		self._get_grains_()

		print("\nNum Grains = %s" % len(self.grains))

		for grain_id in self.grains:
			xyz = all_mol_crds[self.grains[grain_id]]
			ms = np.sqrt(len(xyz)) * 3
			mean = np.mean(xyz, axis=0)
			ax.plot([mean[0]], [mean[1]], [mean[2]], 'ro', ms=ms)



	def _follow_dict_chain_(self, mol_couplings, mol, mols_clustered=False):
		"""Will cluster the high coupling mols from a dict."""
		if mols_clustered is False: mols_clustered = []

		if mol not in mols_clustered: mols_clustered.append(mol)

		for new_mol in mol_couplings[mol]:
			if new_mol not in mols_clustered:
				mols_clustered.append(new_mol)
				mols_clustered = self._follow_dict_chain_(mol_couplings, new_mol, mols_clustered)

		return mols_clustered