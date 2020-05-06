#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate couplings from hamiltonian data
"""
import numpy as np
import pandas as pd	
import re, os

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

from statistics import NormalDist



coupling_pairs = {
"P_dir" : [(i, i+40) for i in range(1, 761)],
"T1_dir" : [(i, i+3) for i in range(1, 798, 2)] + [(i, i+41) for i in range(1, 760, 2)],
"T2_dir" : [(i, i+1) for i in range(1, 800, 2)] + [(i, i+43) for i in range(1, 760, 2)],
}
class TMP(gen_calc.Calc_Type):
	"""
	A class to be alterred however anyone wants to.

	This is just a convienent way to use your own code that is
	too complex for a little inline-script etc... but not useful
	as a repeatable function.

	To use just edit the _calc_ function below.
	"""
	required_calc = ('mol_layers')

	def _calc_coupling_histograms_(self):
		
		print("Sorting the couplings into directions")

		data = {i: [] for i in coupling_pairs}
		for coup_dir in coupling_pairs:
			for (imol1, imol2) in coupling_pairs[coup_dir]:
				i = imol1 - 1
				j = imol2 - 1
				if i not in self.all_AOM_couplings.data[0]:
					if j not in self.all_AOM_couplings.data[0]:
						print(f"Can't find {i}, {j}")
						continue

					couplings = self.all_AOM_couplings.data[0][j]
					next_mol = i

				else:
					couplings = self.all_AOM_couplings.data[0][i]
					next_mol = j

				if next_mol not in couplings:
					print(f"Can't find {i}, {j}")
					continue

				data[coup_dir].append(couplings[next_mol])


		for i in data:
			plt.hist(data[i], bins="auto")

			# data = [0.7237248252340628, 0.6402731706462489, -1.0616113628912391, -1.7796451823371144, -0.1475852030122049, 0.5617952240065559, -0.6371760932160501, -0.7257277223562687, 1.699633029946764, 0.2155375969350495, -0.33371076371293323, 0.1905125348631894, -0.8175477853425216, -1.7549449090704003, -0.512427115804309, 0.9720486316086447, 0.6248742504909869, 0.7450655841312533, -0.1451632129830228, -1.0252663611514108]
			norm = NormalDist.from_samples(data[i])
			# NormalDist(mu=-0.12836704320073597, sigma=0.9240861018557649)
			print(norm.mean)
			# -0.12836704320073597
			print(norm.stdev)
			print(norm)

		plt.show()









	def _calc__crystal_layers_(self):
		"""
		Will divide the system up in slices determined by the layering algorithm.

		This was made as a quick hack to make the crystal system for FSSH work.
		"""
		crds = self.Var.data.get_xyz_data()[0][0]
		cols = self.Var.data.get_xyz_cols()[0][0]

		natom = len(crds)
		run_inp = self.Var.data['cp2k_inp'].file_data['run.inp']
		top_charge = self.Var.data['cp2k_inp'].file_data['TOPOLOGY-CHARGE-ONLY.inp']
		top_neutral = self.Var.data['cp2k_inp'].file_data['TOPOLOGY-NEUTRAL-ONLY.inp']

		mol_crds = mol_utils.atoms_to_mols(crds, 36)
		mol_cols = mol_utils.cols_to_mols(cols, 36)
		nmol = len(mol_crds)

		COM = mol_utils.get_COM_split_mols(mol_crds, mol_cols)
		
		system_info = geom.get_system_size_info(COM)

		dim = 'z'
		idim = {'x': 0, 'y': 1, 'z': 2}
		windows = self.mol_layers.layer_starts

		for i in range(len(windows)-1):
			start, end = windows[i], windows[i+1]

			mask = (COM[:, idim[dim]] > start) & (COM[:, idim[dim]] < end)

			masked_crds = COM[mask]
			mol_nums = np.arange(len(COM))[mask]
			
		 	# Create the directory.
			dir_ = f"crystal_layers/Layer_{start:.2f}<{dim}<{end:.2f}"
			if not os.path.isdir(dir_):
				os.makedirs(dir_)

			# Create the DECOMP file.
			cp2k_inp.create_DECOMP_inp(mol_nums, 36, f"{dir_}/DECOMP.inp")

			# Create the run.inp file.
			run_inp.change_param(['MOTION', 'MD', 'STEPS'], 3000)
			run_inp.change_param(['FORCE_EVAL', 'MIXED', 'ADIABATIC', 'AOM', 'NUMBER_DIABATIC_STATES'],
								  len(mol_nums))

			run_inp.write(f"{dir_}/run.inp")

			# Create the AOM_COEFF.include file
			with open("./src/data/AOM_COEFFs/Pentancene_Lammps_Single_Mol.include", "r") as f:
				cp2k_inp.create_AOM_include(nmol, mol_nums, 36, f.read(), f"{dir_}/AOM_COEFF.include")


			# Create the TOPOLOGY*.inp
			for fpath, inp in (("crystal_layers/Shared_Files/TOPOLOGY-CHARGE-ONLY.inp", top_charge),
						("crystal_layers/Shared_Files/TOPOLOGY-NEUTRAL-ONLY.inp", top_neutral)):
				
				inp.remove(['CELL', 'ABC'])
				inp.remove(['CELL', 'ALPHA_BETA_GAMMA'])

				inp.change_param(['CELL', 'PERIODIC'], "NONE")
				inp.change_param(['TOPOLOGY', 'NUMBER_OF_ATOMS'], natom)
				inp.change_param(['TOPOLOGY', 'MOL_SET', 'MOLECULE', 'NMOL'], nmol)
				inp.change_param(['TOPOLOGY', 'COORD_FILE_NAME'], "crystal.xyz")

				# inp.add(['CELL', 'C'], '  '.join(self.Var.data['lammps_dump'].metadata['c'].astype(str)))
				# inp.add(['CELL', 'B'], '  '.join(self.Var.data['lammps_dump'].metadata['b'].astype(str)))
				# inp.add(['CELL', 'A'], '  '.join(self.Var.data['lammps_dump'].metadata['a'].astype(str)))

				inp.write(fpath)
				
			# self._plot_xyz_data(COM[mask])
			# plt.show()

	def _calc_no_layers_(self):
		"""
		Will just divide up the system in 6 even slices in the Z axis.
		"""
		crds = self.Var.data['xyz'].xyz_data[0]
		cols = self.Var.data['xyz'].cols[0]
		run_inp = self.Var.data['cp2k_inp'].file_data['run_inp']

		mol_crds = mol_utils.atoms_to_mols(crds, 36)
		mol_cols = mol_utils.cols_to_mols(cols, 36)
		nmol = len(mol_crds)
		print(nmol)

		COM = mol_utils.get_COM_split_mols(mol_crds, mol_cols)
		
		system_info = geom.get_system_size_info(COM)

		dim = 'x'
		idim = {'x': 0, 'y': 1, 'z': 2}
		windows = np.linspace(system_info[f'{dim}min'], system_info[f'{dim}max'], 8)
		for i in range(1, len(windows)-2):
			start, end = windows[i], windows[i+1]

			mask = (COM[:, idim[dim]] > start) & (COM[:, idim[dim]] < end)

			masked_crds = COM[mask]
			mol_nums = np.arange(len(mol_cols))[mask]
			
			# Create the directory.
			dir_ = f"yz_planes/Layer_{start:.2f}<{dim}<{end:.2f}"
			if not os.path.isdir(dir_):
				os.makedirs(dir_)

			# Create the DECOMP file.
			cp2k_inp.create_DECOMP_inp(mol_nums, 36, f"{dir_}/DECOMP.inp")

			# Create the run.inp file.
			cp2k_inp.change_param(run_inp, ['MOTION', 'MD', 'STEPS'], 2000)
			cp2k_inp.change_param(run_inp,
								  ['FORCE_EVAL', 'MIXED', 'ADIABATIC', 'AOM', 'NUMBER_DIABATIC_STATES'],
								  len(mol_nums))
			cp2k_inp.write_inp(run_inp, f"{dir_}/run.inp")

			# Create the AOM_COEFF.include file
			with open("./src/data/AOM_COEFFs/Pentancene_Lammps_Single_Mol.include", "r") as f:
				cp2k_inp.create_AOM_include(len(mol_crds), mol_nums, 36, f.read(), f"{dir_}/AOM_COEFF.include")


