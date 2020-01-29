import load_xyz
import topology
import geom

import numpy as np
import matplotlib.pyplot as plt
import os


def get_long_ax_angles(ats, crds, long_ax_ats=[5, 23], ats_per_mol=36):
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
    all_mags += 1e-12  # Just to remove silly numerical errors in arccos
    all_dots = np.sum(long_ax_vecs * center_vec, axis=1)

    return np.arccos(all_dots / all_mags), long_ax_vecs
    

def get_short_axes_angles(ats, crds, short_ax_ats=[[18, 19], [34, 10], [3, 6], [24, 23], [25, 22],
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
    all_mags += 1e-12  # Just to remove silly numerical errors in arccos

    return np.arccos(all_dots / all_mags), avg_short_vec

        
def get_perp_angles(ats, crds, short_ax_vecs, long_ax_vecs, ats_per_mol=36):
    """
    Will get all the out-of-plane perpendicular angles to the short axis and long axis vectors.

    Inputs:
        * short_ax_vecs => A 2D array of shape (num_mol, 3) with the short axis vector for each molecule
        * long_ax_vecs => A 2D array of shape (num_mol, 3) with the long axis vector for each molecule
    Outputs:
        A 2D array of shape (num_mol, 3) with the perpendicular vector to the short and long axes.
    """
    perp_vecs = np.cross(short_ax_vecs, long_ax_vecs, axis=1)

    mol_COMs = topology.get_mol_COMs(crds, ats, ats_per_mol)
    mol_ind, _ = geom.find_center_atom(mol_COMs)
    center_vec = perp_vecs[mol_ind]

    all_dots = np.sum(perp_vecs * center_vec, axis=1)
    all_mags = np.linalg.norm(perp_vecs, axis=1) * np.linalg.norm(center_vec)
    all_mags += 1e-12  # Just to remove silly numerical errors in arccos

    return np.arccos(all_dots/all_mags), perp_vecs

   
def plot_hist_percent(data, a=False):
      """
      Will plot the data as a histogram
      """
      if not a:
         f, a = plt.subplots()
   
      frq, edges = np.histogram(data)
      edges *= 180./np.pi
      #frq *= 100.
      bars = a.bar(edges[:-1], frq, width=np.diff(edges))
      return bars


def plot_angle_dist_from_file(xyz_filepath):
   """
   Will plot the angle distribution of the xyz coords for the long axis, short axis and perp axis of each molecule.

   Inputs:
      * xyz_filepath => The file pointing towards the data to be plotted

   Outputs:
      figure, axis and title.
   """
   # Actually calculate the data
   ats, crds = load_xyz.read_1_step_xyz(xyz_filepath)
   num_ats_per_mol = topology.get_num_ats_per_mol(crds)
   print("\n\nNumber of atoms per molecule = %i\n\n" % num_ats_per_mol)
   
   short_ax_ang, short_ax_vec = get_short_axes_angles(ats, crds)#, [[0, 4], [11, 1], [2, 9], [3, 13], [10, 13]], 18)
   long_ax_ang, long_ax_vec = get_long_ax_angles(ats, crds)#, [4, 13], 18)
   perp_ax_ang, perp_vec = get_perp_angles(ats, crds, short_ax_vec, long_ax_vec)#, 18)
   
   # set the style of the axes and the text color
   plt.rcParams['axes.edgecolor']='#333F4B'
   plt.rcParams['axes.linewidth']=0.8
   plt.rcParams['xtick.color']='#333F4B'
   plt.rcParams['ytick.color']='#333F4B'
   plt.rcParams['text.color']='#333F4B'
   
   f, a = plt.subplots(3)
   barsS = plot_hist_percent(short_ax_ang, a[0])
   barsL = plot_hist_percent(long_ax_ang, a[1])
   barsP = plot_hist_percent(perp_ax_ang, a[2])
   
   a[0].set_ylabel("Short Ax (%)", fontsize=18)
   a[1].set_ylabel("Long Ax (%)", fontsize=18)
   a[2].set_ylabel("Perpendicular Ax (%)", fontsize=18)
   for ax in a:
      ax.spines['bottom'].set_visible(False)
      ax.yaxis.grid(True, color='#aaaaaa')
      ax.xaxis.grid(False)
      ax.set_facecolor("#FAFAFA")
   a[2].spines['bottom'].set_visible(True)
   a[2].set_xlabel(r"Angle [$^{o}$]", fontsize=20)
   
   title = xyz_filename.replace("_", " ").replace(".xyz", "")
   title = title.replace("quenched", "quenching").replace(" quench ", " quenching ").replace(" quench.", " quenching ")
   title = title.replace("NPT", "").replace("NVT", "")
   
   #title = "Perfect Crystal"
   plt.suptitle(title, fontsize=22)
   plt.tight_layout()
   plt.subplots_adjust(top=0.939,
                       bottom=0.103,
                       left=0.064,
                       right=0.989,
                       hspace=0.238,
                       wspace=0.2)
   return f, a, title
