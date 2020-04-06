#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate couplings from hamiltonian data
"""
import numpy as np
import pandas as pd	
import re

# Plotting tools
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from src.calc import general_calc as gen_calc
from src.calc import molecule_utils as mol_utils
from src.calc import geometry as geom

from src.io_utils import CP2K_inp_files as cp2k_inp

from src.data import consts


class Coupling_Connections(gen_calc.Calc_Type):
	"""
	Will get the distribution of couplings from hamiltonian data for each predefined layer. 
	"""
	required_data_types = ('pos', 'pseudo_ham',)
	required_calc = ("mol_layers", )
	_defaults = {'zmin': 'all', 'zmax': 'all', 'ymin': 'all', 'ymax': 'all', 'xmin': 'all', 'xmax': 'all',
				 "plot_coupling_connections": True, "CC_plot_title": False, "CC_savefig": False}
	required_metadata = ('xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax', "atoms_per_molecule", "reorganisation_energy")
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

		# Get any molecular mapping (if molecules from the decomp section are not spaced by
		#                            a single integer then this should be provided as an index
		#							 mapping. This can be done in the input file as:
		#							 'read "DECOMP.inp" cp2k_inp into var')
		mol_nums = np.arange(len([i for i in couplings if type(i) == int]))
		mol_map = {i: i for i in mol_nums}
		if 'cp2k_inp' in self.Var.data:
			cp2k_params = self.Var.data['cp2k_inp'].all_data['params']
			if 'ENERGY_DECOMP' in cp2k_params.keys():
				decomp_sect = cp2k_params['ENERGY_DECOMP']['INDEX_MOL_DECOMP'].split()
				mol_map = {i: int(decomp_sect[i]) - 1 for i in mol_nums}
				mol_nums = np.array([mol_map[i] for i in mol_nums])

		# Get coupling limit to consider.
		reorg = self.metadata['reorganisation_energy']
		min_Hab = reorg / 100.
		plot_params = {reorg/4.: ({'color': 'r', 'lw': 3}, r"$H_{ab} \geq \frac{\lambda}{4}$"),
					    reorg/10.: ({'color': 'g', 'lw': 1.5}, r"$\frac{\lambda}{4} > H_{ab} \geq \frac{\lambda}{10}$"),
					    min_Hab: ({'color': 'b', 'lw': 0.3}, r"$\frac{\lambda}{10} > H_{ab} \geq \frac{\lambda}{100}$"),
					  }

		# Loop over all the files that contain xyz data.
		for ifile in range(len(xyz_data)):
			# Loop over all steps in the xyz data.
			for istep in range(len(xyz_data[ifile])):

				# Do some data reshaping
				step_data = xyz_data[ifile][istep]
				cols = cols[ifile][istep]
				
				mol_crds = mol_utils.atoms_to_mols(step_data, self.metadata['atoms_per_molecule'])
				mol_col = mol_utils.cols_to_mols(cols, self.metadata['atoms_per_molecule'])

				COM = mol_utils.get_COM_split_mols(mol_crds, mol_col)

				# Loop over coupling mol nums (these are the integer indices in the coupling dict.)
				graph_data = []
				for mol1 in couplings:
					if type(mol1) != int: continue
					mol_couplings = couplings[mol1]

					# These only contain site-energies -this is just a slight optimisation.
					if len(mol_couplings) == 1: continue

					# Loop over the mols this mol is coupled with.
					for mol2 in mol_couplings:
						if mol1 != mol2:
							Hab = mol_couplings[mol2] * consts.Ha_to_meV
							if Hab < min_Hab: continue

							for max_coup in plot_params:
								if Hab >= max_coup:
									plot_args = plot_params[max_coup][0]
									break
							else:
								raise SystemExit("Something went wrong categorising the coupling colors")

							pos1, pos2 = self.mol_layers.rotated_COM[mol_map[mol1]], self.mol_layers.rotated_COM[mol_map[mol2]]
							point_data = {'pos': ((pos1[0], pos2[0]),
												  (pos1[1], pos2[1]),
												  (pos1[2], pos2[2])),
										  'args': plot_args,
										 }
							graph_data.append(point_data)



				do_plot = self.metadata['plot_coupling_connections'] + bool(self.metadata['CC_savefig'])
				if do_plot:
					curr_mol = self.mol_layers.rotated_COM[mol_nums]
					a = self._plot_xyz_data(curr_mol, args={'ls': 'none', 'marker': '.', 'color': 'k'})
					for plot_data in graph_data:
						a.plot(*plot_data['pos'], **plot_data['args'])

					a.view_init(azim=0, elev=90)

					if self.metadata['CC_plot_title']: a.set_title(self.metadata['CC_plot_title'])


					legend_elements = [Line2D([0], [0], label=plot_params[i][1], **plot_params[i][0]) for i in plot_params]

					a.legend(handles=legend_elements, loc="best")
					plt.tight_layout()

					if bool(self.metadata['CC_savefig']):
						plt.savefig(self.metadata['CC_savefig'])
						plt.close()
					else:
						plt.show()



class Couplings(gen_calc.Calc_Type):
	"""
	Will get the distribution of couplings from hamiltonian data.
	"""
	required_data_types = ("pseudo_ham", )
	_defaults = {'plot_global_couplings': False}

	def _calc_(self):
		ham_data = self.Var.data['pseudo_ham'].data

		couplings = [j for i in ham_data for j in i['couplings']]
		self.couplings = np.array(couplings) * consts.Ha_to_meV

		bins, bin_edges = np.histogram(self.couplings, bins='auto')

		print(f"Nbins = {len(bins)}")
		self.distribution = pd.DataFrame({'bins': bins, 'bin_edges': bin_edges[:-1]})
		self.smoothed_distribution = self.distribution.rolling(10, center=True).mean()

		print(self.Var.data['pseudo_ham'].filepath)

		if self.metadata['plot_global_couplings']:
			plt.show()
			# f.savefig(f'Coupling_Layer_{re.sub("[a-zA-Z ]", "", self.Var["title"])}.png')
			# plt.close()

	def plot_couplings(self, a=False):
		if a is False:
			f, a = plt.subplots()

		a.plot(self.smoothed_distribution['bin_edges'], self.smoothed_distribution['bins'], 'k-')
		# a.ylim([0, 30])
		a.set_ylabel(r"Num Counts")
		a.set_xlabel(r"H$_{ab}$ [meV]")
		a.set_title(self.Var['title'])
		ylim = a.get_ylim()
		a.set_ylim([ylim[0], ylim[1]/6])
		# a.set_xlim([-500, 300])

		yticks, xticks = a.get_yticks(), a.get_xticks()
		max_y = np.max(self.smoothed_distribution['bins'])
		x_peak  = self.smoothed_distribution['bin_edges'][self.smoothed_distribution['bins'] == max_y]
		# a.annotate(f"Max Peak at ({float(x_peak):.1f}, {int(round(max_y, 1)):.1f})",
		# 			(xticks[1], yticks[-2]), fontsize=20)
		plt.tight_layout()

		return a
