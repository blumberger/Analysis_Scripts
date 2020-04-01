#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate the average rotation of the long axis of a molecular system
"""
import numpy as np
import pandas as pd

from src.calc import general_calc as gen_calc
from src.calc import geometry as geom

class Long_Ax_Rot(gen_calc.Calc_Type):
    """
    Will calculate the average long axis rotation of a molecular system. The result will be a 3-vector.

    This calculates the rotation angle by calculating the histogram each component of the long axis
    rotation vectors then finding the maximum.

    Inputs:
        * Variable <Variable> => An instance of the Variable class

    Important Attributes:
        * required_metadata <tuple> => Any keys that are required in the metadata dictionary.
        * required_calc <tuple> => Any values that need calculating to calculate this value.
        * data <*> => The data that has been calculated.
    """
    required_metadata = ("atoms_per_molecule", "number_atoms")
    required_data_types = ('pos',)
    required_calc = ("angular_dist",)

    # Need these 3 attributes to create a new variable type
    metadata = {}
    name = "Long Axis Rotation"

    def _calc_(self):
        """
        Will do the calculating and return the standard long axis vector.
        """
        counts_x, bin_x = np.histogram(np.abs(self.angular_dist.long_vecs[:,0]), bins='auto')
        counts_y, bin_y = np.histogram(np.abs(self.angular_dist.long_vecs[:,1]), bins='auto')
        counts_z, bin_z = np.histogram(np.abs(self.angular_dist.long_vecs[:,2]), bins='auto')

        max_bin_x = bin_x[np.argmax(counts_x)]
        max_bin_y = bin_y[np.argmax(counts_y)]
        max_bin_z = bin_z[np.argmax(counts_z)]

        self.long_ax_vector = np.array([max_bin_x, max_bin_y, max_bin_z])

        xy_plane = np.array([0, 0, 1])
        dot_prod = np.dot(self.long_ax_vector, xy_plane)
        mags = np.linalg.norm(self.long_ax_vector)

        self.yz_plane_angle = geom.angle_between_vecs(self.long_ax_vector, [1, 0, 0])
        self.xz_plane_angle = geom.angle_between_vecs(self.long_ax_vector, [0, 1, 0])
        self.xy_plane_angle = geom.angle_between_vecs(self.long_ax_vector, [0, 0, 1])

        self.yz_rotation_matrix = geom.map_vec1_to_unit(self.long_ax_vector, [1, 0, 0])
        self.xz_rotation_matrix = geom.map_vec1_to_unit(self.long_ax_vector, [0, 1, 0])
        self.xy_rotation_matrix = geom.map_vec1_to_unit(self.long_ax_vector, [0, 0, 1])

        return self.long_ax_vector



