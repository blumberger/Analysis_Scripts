#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate nearest neighbour lists
"""
import numpy as np
import json

# Import calculating functions
from src.calc import general_types as gen_type
from src.calc import geometry as geom
from src.calc import molecule_utils as mol_utils

from src.data import consts

# import type checking functions
from src.system import type_checking as type_check

class Angular_Dist(gen_type.Calc_Type):
    """
    Will calculate the angular distributions of the molecules in a system.

    Inputs:
        * Variable <Variable> => An instance of the Variable class

    Important Attributes:
        * required_metadata <tuple> => Any keys that are required in the metadata dictionary.
        * required_calc <tuple> => Any values that need calculating to calculate this value.
        * data <*> => The data that has been calculated.
    """
    required_metadata = ('long_axis_atoms', 'short_axis_atoms',
                         'atoms_per_molecule')

    # Need these 3 attributes to create a new variable type
    data = {}
    metadata = {'file_type': 'json'}
    name = "Angular Distributions"
    with open(consts.PT_FILEPATH) as f: PT = json.load(f)

    def calc(self):
        """
        Will calculate the angular distribution of the molecular system.

        This will loop over each molecule, find the rotation (wrt long ax, short
        ax of central molecule) of the long and short axes of the molecule then
        create a histogram of this data.
        """
        ats_per_mol = self.Var.metadata['atoms_per_molecule']
        all_at_crds = self.Var.data.xyz_data

        for at_crds in all_at_crds:
            # Divide atomic coordinates into molecular coordinates
            mol_crds = mol_utils.atoms_to_mols(at_crds, ats_per_mol,
                                                    nstep=len(at_crds))

            # Get the center mol (the one to compare to)
            avg_mol_crds = np.mean(mol_crds, axis=1)
            cent_ind, center = geom.find_center_atom(avg_mol_crds)

            # Get the long axis atoms
            long_ax_ats = self.Var.data.metadata['long_axis_atoms']
            short_ax_ats = self.Var.data.metadata['short_axis_atoms']

            long_axes = self.get_displacement_vecs(mol_crds, long_ax_ats)
            break

    def get_displacement_vecs(self, mol_crds, at_inds):
        """
        Will get the vectors describing the displacement between 2 atoms.

        Inputs:
            * mol_crds <np.NDArray> => (nmol, nat_per_mol, 3) The molecular coordinates
            * at_inds <list<int>> => The atoms to get the displacement vec for.
        Outputs:
            <np.NDArray> (nmol, 3) The displace from at1 to at2 for each mol.
        """
        disp_vec = mol_crds[:, at_inds[0], :] - mol_crds[:, at_inds[1], :]
        return disp_vec
