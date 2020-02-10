#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate pvecs from some xyz coordinates.
"""
import numpy as np

from src.calc import general_types as gen_type

class PVecs(gen_type.Calc_Type):
    """
    Will calculate the pvecs from the data contained within a data file.

    Inputs:
        * Variable <Variable> => A ass or a class that has been derived from Data_File

    Important Attributes:
        * required_metadata <tuple> => Any keys that are required in the metadata dictionary.
    """
    required_metadata = ('num_at_per_mol',)

    def calc(self):
        """
        Will calculate the pvecs from self.Data_File.numeric_data
        """
        XYZFile = self.Var.data
        at_crds = XYZFile.numeric_data
        natom_per_mol = self.Var.metadata['num_at_per_mol']
        nstep = XYZFile.nstep
        
        # First remove 'Ne' atoms
        cols = XYZFile.cols
        at_crds = np.array([i[cols[0] != 'Ne'] for i in XYZFile.numeric_data])
        print(at_crds)



