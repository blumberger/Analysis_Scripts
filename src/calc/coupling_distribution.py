#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate couplings from hamiltonian data
"""
import numpy as np
import pandas as pd	
import re

from src.calc import general_calc as gen_calc
from src.data import consts


class Couplings(gen_calc.Calc_Type):
	"""
	Will get the distribution of couplings from hamiltonian data.
	"""
	required_data_types = ("pseudo_ham", )

	def _calc_(self):
		ham_data = self.Var.data['pseudo_ham'].data

		couplings = [j for i in ham_data for j in i['couplings']]
		self.couplings = np.array(couplings) * consts.Ha_to_meV

		bins, bin_edges = np.histogram(self.couplings, bins='auto')

		print(f"Nbins = {len(bins)}")
		self.distribution = pd.DataFrame({'bins': bins, 'bin_edges': bin_edges[:-1]})
		self.smoothed_distribution = self.distribution.rolling(10, center=True).mean()

		print(self.Var.data['pseudo_ham'].filepath)

		import matplotlib.pyplot as plt

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
		# plt.show()
		f.savefig(f'Coupling_Layer_{re.sub("[a-zA-Z ]", "", self.Var["title"])}.png')
		plt.close()
