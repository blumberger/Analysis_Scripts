#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate nearest neighbour lists
"""
import numpy as np

from src.calc import general_types as gen_type

class NN(gen_type.Calc_Type):
    """
    Will calculate the Nearest neighbour list from the data contained within a data file.

    Inputs:
        * Variable <Variable> => A ass or a class that has been derived from Data_File

    Important Attributes:
        * required_metadata <tuple> => Any keys that are required in the metadata dictionary.
        * required_calc <tuple> => Any values that need calculating to calculate this value.
        * data <*> => The data that has been calculated.
    """
    required_metadata = ()

    # Need these 3 attributes to create a new variable type
    data = {}
    metadata = {'file_type': 'json'}
    name = "Nearest Neighbour Calculator"

    def calc(self):
        """
        Will calculate the pvecs from self.Data_File.numeric_data
        """
        XYZFile = self.Var.data
        cols = XYZFile.cols
        at_crds = np.array([i[cols[0] != 'Ne'] for i in XYZFile.numeric_data])
        natom = len(at_crds[0])
        nstep = len(at_crds)
       
        # Calculate the nearest neighbour lists for each step
        for step in range(nstep):
            crds = at_crds[step]
            NN_step = {'distances': {},
                       'atom_indices': {}}

            # Loop over upper triangle of atom pairs
            for iat in range(natom-1):
                # Get the atom indices
                at_inds = np.arange(len(crds))

                # Get distances between atoms (only upper triangle though)
                at_msk = at_inds > iat
                all_dist = crds[at_msk] - crds[iat]
                all_dist = np.linalg.norm(all_dist, axis=1)

                # Sort the data by distance
                at_inds  = at_inds[at_msk]
                sorted_data = sorted(zip(all_dist.tolist(),   # Use tolist for writing json later
                                         at_inds.tolist()))
                at_inds = [i[1] for i in sorted_data]
                all_dist = [i[0] for i in sorted_data]

                # Save the data
                NN_step['distances'][iat] = all_dist
                NN_step['atom_indices'][iat] = at_inds

            self.data[step] = NN_step
        return self.data

