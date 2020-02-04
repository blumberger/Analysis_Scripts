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
    

def get_ax_angle(ats, crds, ax_ats, ats_per_mol=36):
    """
    Will get the angles that the axes defined by the atoms given of each molecule makes with that axis of the
    first molecule selected.
    
    Inputs:
        * atom_types => The atom types (array of strings giving element type)
        * crds => The atomic coordinates of shape (num_atoms, 3)
        * ax_ats => Atom indices on the long axis of the molecule (list of len 2 lists)
        * ats_per_mol => The number of atoms in a single molecule

    Outputs:
        * An array of shape num_ats
    """
    allMolCrds = topology.atoms_to_mols(crds, ats_per_mol)
    mol_COMs = topology.get_mol_COMs(crds, ats, ats_per_mol)
    mol_ind, _ = geom.find_center_atom(mol_COMs)
    ax_ats = np.array(ax_ats)
    
    all_disp_vecs = allMolCrds[:, ax_ats[:, 0]] - allMolCrds[:, ax_ats[:, 1]]
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

   
def plot_hist_percent(data, a=False, nbins=False, legend=""):
      """
      Will plot the data as a histogram
      """
      if not a:
         f, a = plt.subplots()
   
      if nbins: frq, edges = np.histogram(data, bins=nbins, density=True)
      else: frq, edges = np.histogram(data, density=True)
      edges *= 180./np.pi
      #frq *= 100.
      bars = a.bar(edges[:-1], frq, width=np.diff(edges), label=legend)
      return bars


def plot_angle_dist_from_file(xyz_filepath, a=False, nbins=False, legend_label=""):
   """
   Will plot the angle distribution of the xyz coords for the long axis, short axis and perp axis of each molecule.

   Inputs:
      * xyz_filepath => The file pointing towards the data to be plotted
      * a => 3 matplotlib axes on which to plot (top, middle, bottom)

   Outputs:
      figure, axis and title.
   """
   # Actually calculate the data
   ats, crds = load_xyz.read_1_step_xyz(xyz_filepath)
   num_ats_per_mol = 36 ##topology.get_num_ats_per_mol(crds)
   print("\n\nNumber of atoms per molecule = %i\n\n" % num_ats_per_mol)
   
   short_ax_ang, short_ax_vec = get_ax_angle(ats, crds, [[23, 24], [22, 25],
                                                         [21, 26], [20, 27],
                                                         [19, 18], [34, 10],
                                                         [9, 0], [8, 1], [7, 2],
                                                         [6, 3], [5, 4]], 36)
   long_ax_ang, long_ax_vec = get_ax_angle(ats, crds, [[23, 5], [24, 4], [22, 6], [25, 3]], 36)
   perp_ax_ang, perp_vec = get_perp_angles(ats, crds, short_ax_vec, long_ax_vec)#, 18)
   
   # set the style of the axes and the text color
   plt.rcParams['axes.edgecolor']='#333F4B'
   plt.rcParams['axes.linewidth']=0.8
   plt.rcParams['xtick.color']='#333F4B'
   plt.rcParams['ytick.color']='#333F4B'
   plt.rcParams['text.color']='#333F4B'
   
   if len(a) != 3:
      f, a = plt.subplots(3)
   else:
      f = False

   barsS = plot_hist_percent(short_ax_ang, a[0], nbins=nbins, legend=legend_label)
   barsL = plot_hist_percent(long_ax_ang, a[1], nbins=nbins)
   barsP = plot_hist_percent(perp_ax_ang, a[2], nbins=nbins)
   
   a[0].set_ylabel("Short Ax (%)\nDensity", fontsize=18)
   a[1].set_ylabel("Long Ax (%)\nDensity", fontsize=18)
   a[2].set_ylabel("Perpendicular Ax (%)\nDensity", fontsize=18)
   for ax in a:
      ax.spines['bottom'].set_visible(False)
      ax.yaxis.grid(True, color='#aaaaaa')
      ax.xaxis.grid(False)
      ax.set_facecolor("#FAFAFA")
   a[2].spines['bottom'].set_visible(True)
   a[2].set_xlabel(r"Angle [$^{o}$]", fontsize=20)
   
   plt.tight_layout()
   return f, a
