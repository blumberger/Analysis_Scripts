import load_xyz
import topology
import geom

import numpy as np
import matplotlib.pyplot as plt
import os


def get_all_long_ax_angles(ats, crds, long_ax_ats=[5, 23], ats_per_mol=36):
    """
    Will get the angular distribution of the long axes of all molecules with respect to x, y and z.

    Inputs:
        * atom_types => The atom types (array of strings giving element type)
        * crds => The atomic coordinates of shape (num_atoms, 3)
        * long_ax_ats => Atom indices on the long axis of the molecule (list len 2)
        * ats_per_mol => The number of atoms in a single molecule

    Output:
        A 2D array of shape (nmol, 3) with angles of rotation wrt the x, y and z axes for each molecule.
    """
    if len(long_ax_ats) != 2:
        raise SystemExit("Please only specify 2 atoms sat on the long axis for the argument 'long_ax_ats' in func 'get_angle_distribution'")

    allMolCrds = topology.atoms_to_mols(crds, ats_per_mol)
    mol_COMs = topology.get_mol_COMs(crds, ats, ats_per_mol)
    mol_ind, _ = geom.find_center_atom(mol_COMs)
    at1, at2 = long_ax_ats[0], long_ax_ats[1]

    long_ax_vecs = allMolCrds[:, at1] - allMolCrds[:, at2] 
    center_vec = long_ax_vecs[mol_ind]

    all_mags = np.linalg.norm(long_ax_vecs, axis=1) * np.linalg.norm(center_vec)
    all_dots = np.sum(long_ax_vecs * center_vec, axis=1)
    return np.arccos(all_dots / all_mags), long_ax_vecs
    

def get_short_axes_angles(ats, crds, short_ax_ats=[[18, 19], [34, 10], [3, 6], [24, 23], [25, 22]
                                                   [20, 27], [9, 1]], ats_per_mol=36):
    """
    Will get the angles that the short axes vector of each molecule makes with the short axes of the
    first molecule selected.
    
    Inputs:
        * atom_types => The atom types (array of strings giving element type)
        * crds => The atomic coordinates of shape (num_atoms, 3)
        * long_ax_ats => Atom indices on the long axis of the molecule (list of len 2 lists.
        * ats_per_mol => The number of atoms in a single molecule

    Outputs:
        * An array of shape num_ats
    """
    allMolCrds = topology.atoms_to_mols(crds, ats_per_mol)
    mol_COMs = topology.get_mol_COMs(crds, ats, ats_per_mol)
    mol_ind, _ = geom.find_center_atom(mol_COMs)
    short_ax_ats = np.array(short_ax_ats)
    
    all_disp_vecs = allMolCrds[:, short_ax_ats[:, 0]] - allMolCrds[:, short_ax_ats[:, 1]]
    avg_short_vec = np.mean(all_disp_vecs, axis=1)
    center_vec = avg_short_vec[mol_ind]

    all_dots = np.sum(avg_short_vec * center_vec, axis=1)
    all_mags = np.linalg.norm(avg_short_vec, axis=1) * np.linalg.norm(center_vec)

    return np.arccos(all_dots / all_mags), avg_short_vec

        
def get_all_perp_angles(ats, crds, short_ax_vecs, long_ax_vecs):
    """
    Will get all the out-of-plane perpendicular angles to the short axis and long axis vectors.

    Inputs:
        * short_ax_vecs => A 2D array of shape (num_mol, 3) with the short axis vector for each molecule
        * long_ax_vecs => A 2D array of shape (num_mol, 3) with the long axis vector for each molecule
    Outputs:
        A 2D array of shape (num_mol, 3) with the perpendicular vector to the short and long axes.
    """
    perp_vecs = np.cross(short_ax_vecs, long_ax_vecs, axis=0)

    allMolCrds = topology.atoms_to_mols(crds, ats_per_mol)
    mol_COMs = topology.get_mol_COMs(crds, ats, ats_per_mol)
    mol_ind, _ = geom.find_center_atom(mol_COMs)
    center_vec = perp_vecs[mol_ind]

    all_dots = np.sum(avg_short_vec * center_vec, axis=1)
    all_mags = np.linalg.norm(perp_vecs, axis=1) * np.linalg.norm(center_vec)

    return np.arccos(all_dots/all_mags), perp_vecs



xyz_folder = "/home/kev-boi-2k18/Documents/PhD/AMORPHOUS_PENTACENE/100ps_quenched_NVT_equilibriation"
xyz_filename = "100ps_quench.xyz"
xyz_filepath = os.path.join(xyz_folder, xyz_filename)
if not os.path.isfile(xyz_filepath):
    raise SystemExit("Can't find the file: '%s'!" % xyz_filepath)


