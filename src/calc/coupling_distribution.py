#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate couplings from hamiltonian data
"""
import numpy as np
import pandas as pd	

from src.calc import general_types as gen_type


class Couplings(gen_type.Calc_Type):
	"""
	Will get the distribution of couplings from hamiltonian data.
	"""
	required_data_types = ("pseudo_ham", )

	def calc(self):
		ham_data = self.Var.required_data.get_data()

		couplings = [j for i in ham_data for j in i['couplings']]
		self.couplings = np.array(couplings)

		bins, bin_edges = np.histogram(self.couplings, bins="auto", density=True)

		self.distribution = pd.DataFrame({'bins': bins, 'bin_edges': bin_edges[:-1]})
		self.smoothed_distribution = self.distribution.rolling(100, center=True).mean()

		# import matplotlib.pyplot as plt

		# plt.plot(rolling['bin_edges'], rolling['bins'], 'k-')
		# plt.ylim([0, 30])
		# plt.ylabel(r"Count Density")
		# plt.xlabel(r"H$_{ab}$ [Ha]")
		# plt.title("Global Coupling Distribution")
		# plt.show()
		# plt.hist(couplings, bins=1000, density=True)
