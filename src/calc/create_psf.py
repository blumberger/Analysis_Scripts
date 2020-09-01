#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to create a psf file from an xyz file.
"""
import numpy as np

from src.calc import general_calc as gen_calc
from src.calc import molecule_utils as mol_utils

from src.system import type_checking as type_check


def calc_all_dists(pos1, all_pos, types=False):
    dist = np.linalg.norm(all_pos - pos1, axis=1)
    if types is not False:
        sorting = sorted(zip(dist, np.arange(1, len(dist)+1), types))
    else:
        sorting = sorted(zip(dist, np.arange(1, len(dist)+1)))
    return sorting


class Create_PSF(gen_calc.Calc_Type):
    """
    Will create a psf file from
    """
    required_metadata = ('atoms_per_molecule', 'number_each_atom', 'bonds',
                         'dihedrals', 'angles', 'atom_types',)
    required_data_types = ('pos', )
    _write_types = ('psf', )

    # Need these 3 attributes to create a new variable type
    metadata = {'file_type': 'psf'}
    name = "Create PSF File"

    def _calc_(self):
        """
        Will call the relevant functions to calculate the psf file info.

        This will calculate:
            * Bonds => Which atom is bonded to which. We assume elements num of
                        bonds are given by the periodic table.
            * Angles => Which elements form the angles are given by the metadata
            * Dihedrals => Which elements form the dihedrals are given by the metadata
        """
        self.ats_per_mol = self.metadata['atoms_per_molecule']

        self.xyz_data = self.Var.data.get_xyz_data()
        self.cols = self.Var.data.get_xyz_cols(self.metadata['number_each_atom'])

        self.each_file_data = []
        for ifile in range(len(self.xyz_data)):
            file_xyz = self.xyz_data[ifile]
            file_cols = self.cols[ifile]

            self.natom = len(self.xyz_data[0][0])
            print(np.shape(self.xyz_data))
            self.all_mol_crds = mol_utils.atoms_to_mols(file_xyz, self.ats_per_mol)
            nmol = self.all_mol_crds.shape[1]
            self.mol_col = np.reshape(file_cols[0], (nmol, self.ats_per_mol))

            for mol_crd in self.all_mol_crds:
                # Calculate the NN list for 1 molecule
                NN = np.apply_along_axis(calc_all_dists, 1, mol_crd[0], mol_crd[0],
                                         self.mol_col[0])

                # Get the bonding info for the first molecule only. This can be replicated
                #  later
                bond_types = self.get_bonding_types(self.metadata['bonds'])
                all_bonds = mol_utils.get_bonding_info(mol_crd[:1], bond_types,
                                                       self.mol_col, self.metadata['atom_types'],
                                                       NN=NN, cutoff=4)

                bonds = self.remove_bond_info_dupes(all_bonds[0])
                nbonds = len(bonds)

                # Get the angle info
                all_angles = mol_utils.get_topo_info(mol_crd[:1], self.metadata['angles'],
                                                          all_bonds, self.mol_col,
                                                          self.metadata['atom_types'],
                                                          NN=NN)
                angles = self.remove_top_dict_dupes(all_angles[0])
                nangles = len(angles)

                # Get the dihedral info
                all_dihedrals = mol_utils.get_topo_info(mol_crd[:1], self.metadata['dihedrals'],
                                                             all_bonds, self.mol_col,
                                                             self.metadata['atom_types'],
                                                             NN=NN)
                dihedrals = self.remove_top_dict_dupes(all_dihedrals[0])
                ndihedrals = len(dihedrals)

            self.each_file_data.append(
                                    {
                                     'all_bonds': all_bonds, 'nbonds': nbonds, 'bonds': bonds,
                                     'all_angles': all_angles, 'nangles': nangles, 'angles': angles,
                                     'all_dihedrals': all_dihedrals, 'ndihedrals': ndihedrals, 'dihedrals': dihedrals,
                                    }
                                      )
            self.plot_bond_structure(bonds)
            print("\n\nINFO:\nTo see if the bonding looks correct open the bond_struct.png picture!\n\n")


    def remove_top_dict_dupes(self, top_dict):
        """
        Will remove the duplicates from the all angles and dihedrals dictionary.

        This will remove any chains of atom indices that are equivalent. E.g:
            if [2, 1, 11] is already in the list of atoms then we can't allow
            [11, 1, 2] (its reverse).

        Inputs:
            * top_dict <dict> => The dictionary of the topological quantities
        Ouputs:
            <array> Every value in an array where there aren't any items that
                    are the reverse of any other item.
        """
        all_inds = []
        for iat in top_dict:
            for vals in top_dict[iat]:
                if vals not in all_inds and list(reversed(vals)) not in all_inds:
                    all_inds.append(vals)

        return np.array(all_inds)


    def remove_bond_info_dupes(self, all_bonds):
        """
        Will remove duplicate bonds from the all_bonds dict

        Inputs:
            * all_bonds <dict> => The dictionary holding info on all bonds
        """
        all_pairs = []
        for i in all_bonds:
            for j in all_bonds[i]:
                pair = min([i, j]), max([i, j])
                if pair not in all_pairs:
                    all_pairs.append(pair)
        return np.array(all_pairs)



    def get_bonding_types(self, bond_info):
        """
        Will return a dict with every allowed bonding configuration.

        Inputs:
            * bond_info <dict> => The dict declared in the input parameters
        Outputs:
            <dict> bonding types
        """
        bonding_types = {}
        for i in bond_info:
            if len(bond_info[i]) == 2:
                if bond_info[i][0] == bond_info[i][1]:
                    bonding_types.setdefault(bond_info[i][0],
                                             []).append(bond_info[i][1])
                else:
                    bonding_types.setdefault(bond_info[i][0],
                                             []).append(bond_info[i][1])
                    bonding_types.setdefault(bond_info[i][1],
                                             []).append(bond_info[i][0])
            else:
                raise SystemExit("\n\n\n\nThe 'BONDS' section should declare how each "
                                  + "atom bonds with other atoms.\n\n"
                                  + "To do this you should use the format:\n\n\t"
                                  + "BONDS:\n\t1 = <at type> <at type>\n\t"
                                  + "2 = <at type> <at type>\n\t.\n\t.\n\t.\n\n"
                                  + "You haven't declared the correct number of "
                                  + "atom types in the params file, you declared"
                                  + f" {len(bond_info[i])} for var {i}.")
        return bonding_types


    def plot_bond_structure(self, bonds, show=False):
        """
        Will plot the atoms in a 3D scatter plot with bonds that have been calculated plotted too
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure(figsize=(16,9))
        ax = fig.add_subplot(111, projection="3d")

        crds = self.all_mol_crds[0][0]
        cols = [self.metadata['atom_types'][i] for i in self.mol_col[0]]
        masses = np.array([mol_utils.PT_abbrv[i]['atomic_weight'] for i in cols])
        colors = np.array([mol_utils.PT_abbrv[i]['plot_color'] for i in cols])

        # Plot coords in a cube space
        max_len = max(np.abs([max(crds[:, 0]) - min(crds[:, 0]),
                       max(crds[:, 1]) - min(crds[:, 1]),
                       max(crds[:, 2]) - min(crds[:, 2])]))
        min_vals = [min(crds[:, 0]), min(crds[:, 1]), min(crds[:, 2])]

        u_mass = np.unique(masses)
        for m in u_mass:
            mask = masses == m
            c = colors[mask][0]
            new_crds = crds[mask]
            ax.plot(new_crds[:, 0], new_crds[:, 1], new_crds[:, 2], 'o',
                     color=c, ms=(m**0.3)*5, ls='none')

        ax.set_xlim([min_vals[0]-1, min_vals[0]+max_len+1])
        ax.set_ylim([min_vals[1]-1, min_vals[1]+max_len+1])
        ax.set_zlim([min_vals[2]-1, min_vals[2]+max_len+1])

        # Now add the bonds
        for bnd in bonds:
            at1, at2 = crds[bnd[0]-1], crds[bnd[1]-1]
            ax.plot([at1[0], at2[0]], [at1[1], at2[1]],
                    [at1[2], at2[2]], 'k--')

        plt.savefig("bond_struct.png")

    def create_atoms_section(self, spaces=10):
        """
        Will create the atoms section of the PSF file.

        Inputs:
            * space <str> OPTIONAL => The spacing from from the start line to the end of the first int
        Outputs:
            <str> The atoms section as a string
        """
        at_counts = {}
        print(self.natom)
        s = f"{self.natom}".rjust(spaces) + " !NATOM\n"
        for istep in range(len(self.all_mol_crds)):
            for imol, cols in enumerate(self.mol_col):
                for iat, at in enumerate(cols):
                    at_counts[at] = at_counts.setdefault(at, 0) + 1
                    at_num = at_counts[at]
                    name = f"{at}{at_num}".ljust(9)
                    at_type = f"{at}".ljust(6)
                    num = f"{iat + 1}".rjust(spaces)
                    mass = f"{mol_utils.PT_abbrv[self.metadata['atom_types'][at]]['atomic_weight']}".rjust(14)
                    s += f"{num} SYS      {imol+1}        UNK      {name}{at_type}0.000000{mass}"
                    s += "\n"

        return s

    def create_topo_section(self, inds, inds_in_row, spaces=10):
        """
        Will create the bonds, angles, dihedrals etc sections.

        Inputs:
            * inds <array> => The atom indices that should be saved in shape (nat_per_mol, N)
            * inds_in_row <int> => The number of entries in 1 row
            * space <str> OPTIONAL => The spacing from from the start line to the end of the first int
        
        Outputs:
            <str> The section text
        """
        s = ""
        inds = inds.astype(str)
        for i in range(len(inds)):
            if i % inds_in_row == 0: s += "\n"
            for ind in inds[i]:
                s += ind.rjust(spaces)
        return s


    def create_file_str(self):
        """
        Will create the string that can be written as a file.

        Outputs:
            <str> The file text
        """
        len_start = 10

        all_file_txt = []
        for file_data in self.each_file_data:
            bonds, nbonds = file_data['bonds'], file_data['nbonds']
            angles, nangles = file_data['angles'], file_data['nangles']
            dihedrals, ndihedrals = file_data['dihedrals'], file_data['ndihedrals']

            # Title
            s = "PSF\n\n"
            s += "%s !NTITLE\n   PYTHON UTILITY CREATED PSF FILE\n\n" % ("1".rjust(len_start))
            
            # Atoms section
            s += self.create_atoms_section(len_start)

            # Bonds section
            s += "\n" + f"{nbonds}".rjust(len_start) +" !NBOND: bonds\n"
            s += self.create_topo_section(bonds, 4, len_start).lstrip("\n")

            # Angles section
            s += "\n\n" + f"{nangles}".rjust(len_start) +" !NTHETA: angles\n"
            s += self.create_topo_section(angles, 3, len_start).lstrip("\n")

            # Dihedrals section
            s += "\n\n" + f"{ndihedrals}".rjust(len_start) +" !NPHI: dihedrals\n"
            s += self.create_topo_section(dihedrals, 2, len_start).lstrip("\n")

            # The other parts
            for i in ("!NIMPHI: impropers", "!NDON: donors", "!NACC: acceptors", "!NNB", "!NCRTERM: cross-terms"):
                s += "\n\n" + "0".rjust(len_start) + f" {i}"
                
            all_file_txt.append(s)

        return all_file_txt


    def __str__(self):
        """
        Return the file txt
        """
        return "\n\n".join([f"File {i}:"+"\n\n" + f"{txt}" for i, txt in enumerate(self.create_file_str())])
