#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate couplings from hamiltonian data using the AOM method from Orestis' code.
"""
import numpy as np
import pandas as pd
import re
import os
import subprocess
import time
import multiprocessing as mp

from src.calc import general_calc as gen_calc
from src.calc import molecule_utils as mol_utils
from src.calc import geometry as geom

from src.io_utils import CP2K_inp_files as cp2k_inp
from src.io_utils import general_io as gen_io

from src.data import consts

def calc_all_dists(pos_ind, all_pos, cutoff):
	dist = np.linalg.norm(all_pos - all_pos[pos_ind], axis=1)
	mol_nums = np.arange(len(all_pos))[dist < cutoff]
	return mol_nums[mol_nums > pos_ind]


class Calc_All_AOM_Couplings(gen_calc.Calc_Type):
	"""
	Will calculate the couplings between all mols and their nearest neighbours.

	This will take position data for a molecular system and loop over all molecules
	finding their nearest neighbours (within a cutoff). Once these have been found
	a configuration file will be created for each pair for Orestis's AOM coupling
	calculator.


	N.B.

	At the moment the interface with Orestis's AOM coupling calculator is via files
	written to disk. Obviously this isn't a very optimised interface and will be painfully
	slow on computers with slow disk writes so at some point I'll write a C wrapper for
	the C code like the RDF calculator.
	"""
	required_data_types = ('pos',)
	_defaults = {
				 'NN_cutoff': 10, "coupling_calc_exe": "./bin/STO_proj_AOM_overlap",
				 "show_timing": False, "number_processors": False, 'molecule_numbers': 'all',
				 "plot_coupling_connections": False, "save_coupling_connection_plot": False,
				 "plot_coupling_annotations": "auto", "plot_coupling_atoms": "auto",
				 "plot_coupling_mol_numbers": False, 'xmax': False, 'xmin': False,
				 'ymax': False, 'ymin': False, 'zmax': False, 'zmin': False,
				 "coordinate_wrapping": "wrapped", "delete_files": True, "a1_elev": 90,
				 "a2_elev": 0, "a1_azim": 0, "a2_azim": 0, "xy_rotation_matrix": False,
				 "elev_increment":False, "azim_increment":1, "init_elev":30, "init_azim":0,
				 "coupling_do_rotate": False,
				 }
	required_metadata = ("atoms_per_molecule", "AOM_COEFF_file",)

	AOM_mu = {
			  'H': '1.0000',
			  'C': '1.6083  1.0000',
			  'N': '1.9237  1.5000',
			  'O': '2.2458  2.2266',
			  'S': '2.1223  1.8273',
			 }

	metadata = {'file_type': 'json'}
	coup_calc_count = 0

	name = "Create AOM Couplings Config"

	l_col_len = 13
	nproc = mp.cpu_count() - 4
	config_filepaths = tuple(f"tmp_config_{i}.txt" for i in range(nproc))
	dimer_filepaths = tuple(f"tmp_dimer_{i}.xyz" for i in range(nproc))

	times = {}

	def _calc_(self):

		t1 = time.time()


		self.ats_per_mol = self.metadata['atoms_per_molecule']

		# Only act on the first pos file.
		xyz_data = self.Var.data.get_xyz_data()[0]
		if self.metadata['xy_rotation_matrix'] is not False:
			xyz_data = geom.rotate_crds(xyz_data, self.metadata['xy_rotation_matrix'])
		cols = self.Var.data.get_xyz_cols()[0]

		# Get the molecular coordinate
		all_mol_crds_orig = mol_utils.atoms_to_mols(xyz_data, self.ats_per_mol,
													nstep=len(xyz_data))
		self.nmol = np.shape(all_mol_crds_orig)[1]
		self.nstep = len(xyz_data)
		all_mol_cols_orig = mol_utils.cols_to_mols(cols, self.ats_per_mol)


		# Decide which molecules to calculate the coupling for.
		self.all_centroids = mol_utils.get_COM_split_mols(all_mol_crds_orig, all_mol_cols_orig[0])
		all_mol_crds, all_mol_cols = self._get_mol_nums_(all_mol_crds_orig, all_mol_cols_orig)
		self.mol_centroids = np.array([cent[mask] for cent, mask in zip(self.all_centroids, self.mol_nums)])

		self.mol_centroids, mol_mask = self._apply_boundaries_(self.mol_centroids, do_mol_nums=True)
		mol_cols = np.array([col[mask] for col, mask in zip(all_mol_cols, mol_mask)])
		mol_crds = np.array([crds[mask] for crds, mask in zip(all_mol_crds, mol_mask)])
		

		# Decide how many procs to use
		self._set_nproc_(mol_crds)
		self.config_filepaths = tuple(f"tmp_config_{i}.txt"
									  for i in range(self.nproc))
		self.dimer_filepaths = tuple(f"tmp_dimer_{i}.xyz"
									 for i in range(self.nproc))

		# Get the nearest neighbours between molecular centroids (just for writing)
		self.spliced_crds = np.array([data[inds] for data, inds in zip(all_mol_crds_orig, self.mol_nums)])
		self.spliced_cols = np.array([data[inds] for data, inds in zip(all_mol_cols_orig, self.mol_nums)])

		# Get the Nearest neighbour dict
		mol_NN = np.array([[calc_all_dists(pos_ind, step_data, self.metadata['NN_cutoff'])
							for pos_ind in range(len(step_data))]
						   for step_data in self.mol_centroids])

		# Create all the configuration files that will be needed.
		unique_elm = self.__get_unique_elements__(all_mol_cols)
		for dimer_file, config_file in zip(self.dimer_filepaths, self.config_filepaths):
			self._create_config_file_(dimer_file, config_file, unique_elm)

		# Loop over all steps and molecules and calculate couplings
		self.data = [{} for istep in range(len(xyz_data))]
		for istep in range(len(xyz_data)):
			self.step_NN = mol_NN[istep]
			self.num_calcs = sum([len(self.step_NN[i]) for i in self.step_NN])
			self.mol_cols = mol_cols[istep]
			self.mol_crds = mol_crds[istep]

			len_data_div = len(self.step_NN) // self.nproc
			self.tot_mols = len(self.step_NN)
			mol_inds_for_each_proc = [np.arange(len_data_div*(i), len_data_div*(i+1))
									  for i in range(self.nproc)]

			args = [(mol_inds_for_each_proc[i], i) for i in range(self.nproc)]
			if self.nproc > 1:
				pool = mp.Pool(self.nproc)
				all_data = pool.map(self.__calc_mol_couplings__, args)
				for i in all_data:	self.data[istep].update(i)
			else: self.data[istep] = self.__calc_mol_couplings__(args[0])

		self.coup_calc_count = 0
		# For some reason multiprocessing isn't terminating the processes before here.
		if self.nproc > 1: pool.terminate()

		# Show some timing data if necessary
		if self.metadata['show_timing']:
			mean_times = {i: np.mean(self.times[i]) for i in self.times}
			time_per_step = sum(mean_times.values())
			for i in self.times:
				print(f"Timing {i}: {mean_times[i]} s ({100. * mean_times[i] / time_per_step:.1f}%)")

			tot_time = time.time() - t1
			print(f"Total Time Taken: {tot_time:.1f}s")

		if self.metadata['delete_files']:
			os.remove("AOM_dimer.include")
		return self.data

	def _apply_boundaries_(self, mol_centres, do_mol_nums=False):
		"""
		Will splice the molecular centres with the x, y, z min/max.

		This is careful to keep track of the molecular numbers too.

		Inputs:
			* mol_centres <arr> => The centers of each molecule to be spliced of shape (nstep, nmol, 3)

		Outputs:
			* mol_centres <arr> => The spliced mol centres
		"""
		xmin, xmax = self.metadata['xmin'], self.metadata['xmax']
		ymin, ymax = self.metadata['ymin'], self.metadata['ymax']
		zmin, zmax = self.metadata['zmin'], self.metadata['zmax']

		# Add steps to the mol_nums if they aren't there
		self.mol_nums = np.array(self.mol_nums)
		if len(self.mol_nums.shape) == 1:
			self.mol_nums = [self.mol_nums for i in range(len(mol_centres))]

		# Remove the coordinate axis in the mask's shape
		crd_ax = geom._get_crd_ax_(mol_centres, 3)
		shape = list(np.shape(mol_centres))
		del shape[crd_ax]
		shape = tuple(shape)

		# Create the mask
		mask = np.ones(shape, dtype=bool)
		if xmin: mask = mask & (np.take(mol_centres, 0, crd_ax) > xmin)
		if xmax: mask = mask & (np.take(mol_centres, 0, crd_ax) < xmax)
		if ymin: mask = mask & (np.take(mol_centres, 1, crd_ax) > ymin)
		if ymax: mask = mask & (np.take(mol_centres, 1, crd_ax) < ymax)
		if zmin: mask = mask & (np.take(mol_centres, 2, crd_ax) > zmin)
		if zmax: mask = mask & (np.take(mol_centres, 2, crd_ax) < zmax)

		# Apply the mask
		if do_mol_nums:
			self.mol_nums = np.array([nums[mask[i]] for i, nums in enumerate(self.mol_nums)])

		return np.array([cent[mask[i]] for i, cent in enumerate(mol_centres)]), mask

	def _set_nproc_(self, xyz_data):
		""" Will refine the number of processors to use."""
		min_mol_per_proc = 15 # absolute min is actually 16/2 per proc

		if type(self.metadata['number_processors']) == int:
			self.nproc = self.metadata['number_processors']

		# Decide how many procs to use.
		if self.nproc > mp.cpu_count(): self.nproc = mp.cpu_count()
		elif self.nproc < 1: self.nproc = 1
		
		max_proc = np.ceil(xyz_data.shape[1] / min_mol_per_proc)
		if self.nproc > max_proc: self.nproc = max_proc

		if self.nproc < 1: self.nproc = 1

		self.nproc = int(self.nproc)

	def _get_mol_nums_(self, all_mol_crds, all_mol_cols):
		"""
		Will create the molecule numbers array.

		Inputs:
			* all_mol_crds <array> => This is the molecular coordinates array
			* all_mol_crds <array> => This is the element type of each atom in the shape (nmol, nat_per_mol, 3)
		"""
		if type(self.metadata['molecule_numbers']) == str:
			if self.metadata['molecule_numbers'] == 'all':
				# Create the molecule numbers/indicess
				self.mol_nums = np.array([np.arange(len(i)) for i in all_mol_crds])
			else:
				raise SystemExit(f"String argument '{self.metadata['molecule_numbers']}' not understood"
								 + " for parameter 'molecule_numbers'.")

		# If the molecular numbers type is a list then use them as indices
		elif isinstance(self.metadata['molecule_numbers'], (list, type(np.array(1)))):
			all_mol_crds = all_mol_crds[:, self.metadata['molecule_numbers']]
			all_mol_cols = all_mol_cols[:, self.metadata['molecule_numbers']]
			self.mol_nums = np.array([self.metadata['molecule_numbers'] for i in all_mol_crds])

		else:
			raise SystemExit(f"Argument '{self.metadata['molecule_numbers']}' of type"
							 + f"'{type(self.metadata['molecule_numbers'])}' not understood "
							 + " for parameter 'molecule_numbers'.")

		return all_mol_crds, all_mol_cols

	def __calc_mol_couplings__(self, args):#mol_nums, proc_num=0):
		"""
		Will calculate the couplings for N molecules and their nearest neighbours.

		This is a function designed to be used with multiprocessing. This is why the
		input arguments are a little obscure. The inputs should be:
			(mols, proc_num)
		where:
			* mols <list> => list of mol indices to calculate their couplings with their nearest neighbours
			* proc_num <int> => integer giving the process index
		"""
		mol_nums, proc_num = args
		data = {}
		for imol1 in mol_nums:
			nearest_mols = self.step_NN[imol1]

			for imol2 in nearest_mols:
				# The config filepaths is to make life easier later if this needs
				#  parallelising.
				coupling = self.get_coupling(imol1, imol2, self.mol_crds, self.mol_cols,
											 self.config_filepaths[proc_num],
											 self.dimer_filepaths[proc_num])

				data.setdefault(imol1, {})[imol2] = coupling * self.metadata['AOM_scaling_factor']

				self.coup_calc_count += 1
				
				num_done_blocks = int((self.coup_calc_count * self.nproc)/self.num_calcs * 19) + 1
				# num_to_do_blocks = 20 - num_done_blocks
				print("\r" + f"{self.num_calcs} vals to calculate ({self.tot_mols} mols) :   " 
					  + "*"*num_done_blocks, end="\r")

		print("\rFinished                                                               ", end="\r")

		if not os.path.isfile(self.dimer_filepaths[proc_num]):
			print("\n\nWarning no coupling values found within given cutoff"
				  + f"of {self.metadata['NN_cutoff']} for given mols.")
		else:
			if self.metadata['delete_files']:
				os.remove(self.dimer_filepaths[proc_num])

		if self.metadata['delete_files']:
			os.remove(self.config_filepaths[proc_num])

		return data

	def __get_unique_elements__(self, all_mol_cols):
		"""
		Will get the unique elements in the element symbols

		Inputs:
			* all_mol_cols <arr> => The element symbols of shape (nmol, nat_per_mol)
		"""
		# Error check
		for elm in set(all_mol_cols[0, 0]):
			if elm not in self.AOM_mu: raise SystemExit(f"Need coefficients for AOM_mu {elm}")

		# Create all the configuration files
		return tuple(i for i in self.AOM_mu if i in all_mol_cols[0, 0])


	def get_coupling(self, imol1, imol2, mol_crds, mol_cols,
					 config_filepath=False, xyz_filepath=False):
		"""
		Will call Orestis's coupling calculator and return the coupling between mol1 and mol2.

		Steps:
			1) Get mol coords for mol1 and mol2
			2) Write dimer crds as xyz file
			IF config_filepath hasn't been created:
				3) Create config file for coupling calculator
			4) Run coupling calculator and return coupling value.
			5) Tidy up by removing any temporary files created.

		Inputs:
			* imol1 <int> => Molecular index for mol 1
			* imol2 <int> => Molecular index for mol 2
			* mol_crds <array> => Molecular coordinates. Should have shape (nmol, natom_per_mol, 3)
			* mol_cols <array> => Atom element symbols. Should have shape (nmol, natom_per_mol, 3)
			* config_filepath <str> => Where to look for the configuration file
			* xyz_filepath <str> => Where to save the xyz file

		Ouputs:
			<float> The coupling between mol1 and mol2
		"""
		t1 = time.time()
		# Get mol crds
		dimer_crds = mol_crds[[imol1, imol2]]
		dimer_cols = mol_cols[[imol1, imol2]]

		natom = np.shape(dimer_crds)[1] * 2
		t2 = time.time()

		# Write dimer as xyz file
		if xyz_filepath is False: xyz_filepath = "tmp_dimer.xyz"
		self._write_xyz_dimer_(dimer_crds, dimer_cols, xyz_filepath)
		t3 = time.time()

		# Create the config file
		if config_filepath is False:
			config_filepath = "tmp_config.txt"
			unique_elm = self.__get_unique_elements__(mol_cols)
			self._create_config_file_(xyz_filepath, config_filepath, unique_elm)

		exe_filepath = self.metadata['coupling_calc_exe']

		exe_cmd = [exe_filepath, config_filepath]
		res = subprocess.run(exe_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		if res.stderr:
			raise SystemExit(f"Bad mol {imol1} {imol2}")
		t4 = time.time()

		self.times.setdefault('getCrds', []).append(t2 - t1)
		self.times.setdefault('writeCrds', []).append(t3 - t2)
		self.times.setdefault('calcCoupl', []).append(t4 - t3)

		return float(res.stdout.decode("utf-8"))

	def _create_config_file_(self, dimer_file, filepath, unique_elm):
		"""
		Will create the configuration file for the coupling calculator.

		This will also create the AOM_dimer file for the coupling calculator too.
		"""
		# unique_elm = ('H', 'C', 'N', 'O', 'S')
		cols = [
				 ('', ''),
				 ('mode', 'dimer'),
				 ('verb', 'none'),
				 ('', ''),
				 ('dimer', dimer_file),
				 ('AOM_include', "AOM_dimer.include"),
				 ('atoms_frag1', self.ats_per_mol),
				 ('', ''),
				]

		if not os.path.isfile(self.metadata['AOM_COEFF_file']):
			raise SystemError(f"Can't find single mol AOM file. The path {self.metadata['AOM_COEFF_file']} doesn't exist.")

		with open(self.metadata['AOM_COEFF_file'], "r") as f_read:
			AOM_txt = f_read.read() * 2
			with open("AOM_dimer.include", "w") as f_write:
				f_write.write(AOM_txt)

		for elm in unique_elm:
			cols.append((f"AOM_mu  {elm}", self.AOM_mu[elm]))

		s = [(i[0].ljust(self.l_col_len) + str(i[1])).strip() for i in cols]
		with open(filepath, "w") as f:
			f.write('\n'.join(s))


	def _write_xyz_dimer_(self, dimer_crds, dimer_cols, filepath="tmp_dimer.xyz"):
		"""
		Will write a dimer's coords as an xyz file.

		Inputs:
			* dimer_crds <array> => The dimer crds of shape (2, natom_per_mol, 3)
			* dimer_cols <array> => The element symbols of shape (2, natom_per_mol, 3)
		"""
		natom = np.shape(dimer_crds)[1] * 2
		# Reshape the arrays to put them in a more convienent format
		dimer_write_crds = np.reshape(dimer_crds, (natom, 3))
		dimer_write_crds = [''.join([("%.3f" % j).rjust(9) for j in i])
							for i in dimer_write_crds]

		dimer_write_cols = np.reshape(dimer_cols, (natom,))

		# Actually create and write the filetxt
		xyz_body = np.char.add(dimer_write_cols, dimer_write_crds)
		with open(filepath, "w") as f:
			f.write(f"{natom}" + "\n\n" + '\n'.join(xyz_body))

	def json_data(self):
		"""
		Will return the json data.
		"""

		write_data = {}
		for istep, step_data in enumerate(self.data):
			for mol1 in step_data:
				for mol2 in step_data[mol1]:
					step_str = f"step_{istep}"
					if abs(step_data[mol1][mol2] - 0) > 1e-12:
						write_data.setdefault(istep, {}).setdefault(str(mol1), {})[str(mol2)] = step_data[mol1][mol2]

		return write_data

	def get_xyz_data(self):
		"""Will return the xyz data from the class."""
		return np.array([np.reshape(self.spliced_crds,
									(self.spliced_crds.shape[1] * self.spliced_crds.shape[2], 3))])

	def get_xyz_cols(self):
		"""Will return the columns of the xyz data"""
		return np.array([np.reshape(self.spliced_cols,
									(self.spliced_cols.shape[1] * self.spliced_cols.shape[2]))])

	def get_xyz_timesteps(self):
		"""Will return the timesteps or step number if the timesteps aren't available."""
		if hasattr(self.Var.data, 'get_xyz_timesteps'):
			return self.Var.data.get_xyz_timesteps()
		else:
			return [i for i in range(len(self.nstep))]

class Coupling_Connections(gen_calc.Calc_Type):
	"""
	Will get the distribution of couplings from hamiltonian data for each predefined layer.
	"""
	required_data_types = ('pos', 'pseudo_ham',)
	required_calc = ("long_ax_rotation", )
	_defaults = {'zmin': 'all', 'zmax': 'all', 'ymin': 'all', 'ymax': 'all', 'xmin': 'all', 'xmax': 'all',
				 "plot_coupling_connections": True, "CC_plot_title": False, "CC_savefig": False}
	required_metadata = ('xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax', "atoms_per_molecule", "reorganisation_energy")
	name = "Layer Couplings"

	def _calc_(self):
		"""
		Will plot the coupling connections in each layer.
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
			cp2k_params = self.Var.data['cp2k_inp'].data['params']
			if 'ENERGY_DECOMP' in cp2k_params.keys():
				decomp_sect = cp2k_params['ENERGY_DECOMP']['INDEX_MOL_DECOMP'].split()
				mol_map = {i: int(decomp_sect[i]) - 1 for i in mol_nums}
				mol_nums = np.array([mol_map[i] for i in mol_nums])

		# Get coupling limit to consider.
		reorg = self.metadata['reorganisation_energy']
		min_Hab = reorg / 100.
		plot_params = {reorg/2.: ({'color': 'r', 'lw': 3}, r"$H_{ab} \geq \frac{\lambda}{2}$"),
						reorg/10.: ({'color': 'g', 'lw': 1.5}, r"$\frac{\lambda}{2} > H_{ab} \geq \frac{\lambda}{10}$"),
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
				rotated_COM = geom.rotate_crds(COM, self.long_ax_rotation.xy_rotation_matrix)

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

							pos1, pos2 = rotated_COM[mol_map[mol1]], rotated_COM[mol_map[mol2]]
							point_data = {'pos': ((pos1[0], pos2[0]),
												  (pos1[1], pos2[1]),
												  (pos1[2], pos2[2])),
										  'args': plot_args,
										 }
							graph_data.append(point_data)



				do_plot = self.metadata['plot_coupling_connections'] + bool(self.metadata['CC_savefig'])
				if do_plot:
					curr_mol = rotated_COM[mol_nums]
					f = plt.figure()
					a1 = f.add_subplot(121, projection="3d")
					a2 = f.add_subplot(122, projection="3d")
					a1.set_xlabel("X"); a1.set_ylabel("Y"); a1.set_zlabel("Z");
					a1.set_xticks([]); a1.set_yticks([]); a1.set_zticks([]);
					a2.set_xlabel("X"); a2.set_ylabel("Y"); a2.set_zlabel("Z");
					a2.set_xticks([]); a2.set_yticks([]); a2.set_zticks([]);

					self._plot_xyz_data(rotated_COM, a2, args={'color': "k", 'ls': "none", "marker": '.', 'alpha': 0.5})
					self._plot_xyz_data(rotated_COM[mol_nums], a2, args={'color': "r", 'ls': "none", "marker": 'o', 'alpha': 1})

					self._plot_xyz_data(curr_mol, a1, args={'ls': 'none', 'marker': '.', 'color': 'k'})
					for plot_data in graph_data:
						a1.plot(*plot_data['pos'], **plot_data['args'])


					a1.view_init(azim=self.metadata['a1_azim'], elev=self.metadata['a1_elev'])
					a2.view_init(azim=self.metadata['a2_azim'], elev=self.metadata['a2_elev'])

					if self.metadata['CC_plot_title']: a1.set_title(self.metadata['CC_plot_title'].replace("Layer", "").replace("_", " ").strip())


					legend_elements = [Line2D([0], [0], label=plot_params[i][1], **plot_params[i][0]) for i in plot_params]

					a1.legend(handles=legend_elements, loc="best")
					plt.tight_layout()

					if bool(self.metadata['CC_savefig']):
						plt.savefig(self.metadata['CC_savefig'])
						plt.close()
					else:
						plt.show()


					break



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
