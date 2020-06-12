#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module that contains general types that can be subclassed to create new plot types.
"""
import copy
import numpy as np

from src.input_file import input_file_types as inp_types

class Plot_Type(object):
	"""
	A type containing some useful functions that can be inheritted when creating a new calc type.

	Inputs:
		* Variable <Variable> => An instance of the Variable class.

	Public Methods:
		* calc => To be overridden to calculate the property in question.
	"""
	required_metadata = ()
	required_var_attributes = ()

	name = "General Plot Type"

	def __init__(self, Variable):
		"""
		Just check we have the required properties for calculating the quantitiy.
		"""
		self.Var = Variable

		# Copy all data to make it unique for each instance
		all_vars = [i for i in dir(self) if i[0] != '_']
		for i in all_vars:
			if i[0] != '_':
				var = getattr(self, i)
				if not callable(var) and isinstance(var, (dict, list, tuple)):
					setattr(self, i, copy.deepcopy(var))

		# Set the default parameters
		for key in self._defaults:
			if key not in self.Var.metadata:
				self.Var.metadata[key] = self._defaults[key]

			self.metadata[key] = self.Var.metadata[key]

		# Check we have all the data we need to calculate the property
		for key in self.required_metadata:
			if key not in self.Var.metadata:
			   raise KeyError(f"Please load the data '{key}' into the variable '{self.Var.name}'")
			else:
			   self.metadata[key] = self.Var[key]

		# Check the input variable has the required attributes
		for attr in self.required_var_attributes:
			if not hasattr(self.Var, attr):
				if hasattr(self.Var, "data"):
					if type(self.Var.data) == inp_types.Vars:
						found = False
						for key in self.Var.data:
							if hasattr(self.Var.data[key], attr):
								setattr(self.Var, attr, getattr(self.Var.data[key], attr))
								found = True
								break
						if found: continue


					elif hasattr(self.Var.data, attr):
						setattr(self.Var, attr, getattr(self.Var.data, attr))
						continue

					all_attrs = [i for i in dir(self.Var) if i[:2] != "__"]
					msg = f"Can't plot '{self.name}' from the variable '{self.Var.name}'\n\n"
					msg += f"{Variable.name} doesn't have the attribute '{attr}'" + "\n"
					msg += "\nAll Attributes\n\t* " + "\n\t* ".join(all_attrs)
					raise SystemError("\n\n" + msg)

		self._plot_()

	def _plot_(self):
		"""A function to be overriden"""
		print(f"Please override the plot func in {self.name}")

	def _plot_xyz_data(self, ats, a=False, args={'ls': 'None', 'marker': '.'}):
		if a is False:
			f = plt.figure()
			a = f.add_subplot(111, projection="3d")

			a.set_xlabel("X")
			a.set_ylabel("Y")
			a.set_zlabel("Z")

			a.view_init(elev=-3, azim=64)

			a.set_xticks([])
			a.set_yticks([])
			a.set_zticks([])

		a.plot(ats[:, 0], ats[:, 1], ats[:, 2], **args)

		return a

	def _rotate_plot_(self, ax, elev_increment=False, azim_increment=1, init_elev=30, init_azim=0, step_func=False, args=()):
		"""
		Will rotate a 3D plot.

		Inputs:
			* ax <plt.axis> => The axis to rotate
			* elev_increment <float> => The value to increment the elevation axis by
			* azim_increment <float> => The value to increment the azimuthal axis by
			* init_elev <float> => The initial elevation of the axis
			* init_azim <float> => The initial azimuth angle of the axis
			* step_func <function> => A function to carry out at each step
		"""
		count = 0
		if elev_increment and azim_increment:
			for elev_ang in np.arange(init_elev, 360 + init_elev, elev_increment):
				for azim_ang in np.arange(init_azim, 360 + init_azim, azim_increment):
					count += 1
					ax.view_init(azim=azim_ang, elev=elev_ang)
					if step_func:
						step_func(count, *args)

		elif elev_increment and not azim_increment:
			for elev_ang in np.arange(init_elev, 360 + init_elev, elev_increment):
				count += 1
				ax.view_init(azim=init_azim, elev=elev_ang)
				if step_func:
					step_func(count, *args)

		elif azim_increment and not elev_increment:
			for azim_ang in np.arange(init_azim, 360 + init_azim, azim_increment):
				count += 1
				ax.view_init(azim=azim_ang, elev=init_elev)
				if step_func:
					step_func(count, *args)