#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate the angular distributions of molecules in a system.
"""
import numpy as np
import json
import matplotlib.pyplot as plt

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
    _write_types = ('json', )
    required_metadata = ('atoms_per_molecule',)
    _defaults = {'number_bins': 'auto', 'histogram_density': True,
                 'short_axis_atoms': [[5, 6], [4, 7], [3, 8], [1, 28], [0, 18],
                                      [10, 19], [27, 20], [26, 21], [25, 22],
                                      [24, 23]],
                 'long_axis_atoms': [[24, 5], [23, 6]],
                }
    # Need these 3 attributes to create a new variable type
    data = {}
    metadata = {'file_type': 'json'}
    name = "Angular Distribution"
    with open(consts.PT_FILEPATH) as f: PT = json.load(f)

    def calc(self):
        """
        Will calculate the angular distribution of the molecular system.

        This will loop over each molecule, find the rotation (wrt long ax, short
        ax of central molecule) of the long and short axes of the molecule then
        create a histogram of this data.
        """
        self.get_xyz_data()
        ats_per_mol = self.Var.metadata['atoms_per_molecule']
        if 'long_axis_atoms' not in self.Var.metadata:
            long_ax_ats = self.metadata['long_axis_atoms']
        else: long_ax_ats = self.Var.metadata['long_axis_atoms']
        if 'short_axis_atoms' not in self.Var.metadata:
            short_ax_ats = self.metadata['short_axis_atoms']
        else: short_ax_ats = self.Var.metadata['short_axis_atoms']

        # Divide atomic coordinates into molecular coordinates
        all_at_crds = self.compute_data
        all_mol_crds = mol_utils.atoms_to_mols(all_at_crds, ats_per_mol,
                                               nstep=len(all_at_crds))

        # Loop over all steps
        for mol_crds in all_mol_crds:

            # Get the vectors that describe the axes
            self.long_vecs = self.get_all_ax_vecs(long_ax_ats, mol_crds)
            self.short_vecs = self.get_all_ax_vecs(short_ax_ats, mol_crds)
            long_mags = np.linalg.norm(self.long_vecs, axis=1)
            short_mags = np.linalg.norm(self.short_vecs, axis=1)

            # Get the angle of the long and short axis with every other molecule
            nmol = len(mol_crds)
            data_size = int(nmol * (nmol - 1) / 2)
            self.long_angles = np.zeros(data_size)
            self.short_angles = np.zeros(data_size)
            count = 0
            for mol1 in range(1, nmol):
                long_angles = self.get_angles_between_vecs(mol1, self.long_vecs,
                                                           long_mags)
                short_angles = self.get_angles_between_vecs(mol1, self.short_vecs,
                                                            short_mags)

                self.short_angles[count: count + len(short_angles)] = short_angles
                self.long_angles[count: count + len(short_angles)] = long_angles
                count += len(short_angles)

            # Convert to degress cuz they're nicer
            #self.long_angles  *= 180. / np.pi
            #self.short_angles *= 180. / np.pi

            self.long_counts, self.long_bin_edges = np.histogram(self.long_angles,
                                  density=True, bins=self.metadata['number_bins'])
            self.short_counts, self.short_bin_edges = np.histogram(self.short_angles,
                                  density=True, bins=self.metadata['number_bins'])

    def get_angles_between_vecs(self, vec_ind, vecs, mags):
        """
        Will get the angle a vector makes with many other vectors

        Inputs:
            * vec_ind <int> => The index of the vector to compare with
            * vecs <array> => All vectors
            * mags <array> => The magnitude of all the vectors
        Outputs:
            <array> All angles between vecA and other_vecs
        """
        vecA = vecs[vec_ind]
        magA = mags[vec_ind]
        vecs = vecs[:vec_ind]
        mags = mags[:vec_ind]

        dots = np.sum(vecA * vecs, axis=1)
        mags = magA * mags
        angles = dots/mags# * 180/np.pi

        return angles

    def get_all_ax_vecs(self, at_inds, mol_crds):
        """
        Will get the vector that the atoms with ID long_ac_ats of each mol form.

        Inputs:
            * at_inds <list<int>> => the indices of the atoms on the axis you want
            * mol_crds <list<float>> => The molecular coords in shape (nmol, nat per mol, 3)
        """
        if len(np.shape(at_inds)) == 1:
            at_inds = [at_inds]

        ats = mol_crds[:, at_inds, :]
        ats = np.mean(ats, axis=1)
        vecs = ats[:, 0] - ats[:, 1]
        return vecs


    def json_data(self):
        """
        Will return data in a form that the json writer can write.

        This function is just used for writing
        """
        data = {
                'short_ax_angles': self.short_angles.tolist(),
                'short_ax_vecs': self.short_vecs.tolist(),
                'short_ax_histogram': {
                                        'counts': self.short_counts.tolist(),
                                        'bin_edges': self.short_bin_edges.tolist(),
                                      },
                'long_ax_angles': self.long_angles.tolist(),
                'long_ax_vecs': self.long_vecs.tolist(),
                'long_ax_histogram':  {
                                        'counts': self.long_counts.tolist(),
                                        'bin_edges': self.long_bin_edges.tolist(),
                                      },
                }
        return data

    def get_angle_dist(self, mol_crds, at_inds):
        """
        Will get the vectors describing the displacement between 2 atoms.

        Inputs:
            * mol_crds <np.NDArray> => (nmol, nat_per_mol, 3) The molecular coordinates
            * at_inds <list<int>> => The atoms to get the displacement vec for.
        Outputs:
            (<np.NDArray>, <np.array>) The displacement from at1 to at2 for each mol
                                        and just the center mol.
        """
        # Get the center mol (the one to compare to)
        avg_mol_crds = np.mean(mol_crds, axis=1)
        center_ind, _, _ = geom.find_center_atom(avg_mol_crds)

        # First convert to array
        if len(np.shape(at_inds)) == 1:
            at_inds = [at_inds]
        at_inds = np.array(at_inds)

        # First get the vector describing an axis for all mols and the center one
        axis_vecs = mol_crds[:, at_inds[:, 0]] - mol_crds[:, at_inds[:, 1]]
        axis_vecs = np.mean(axis_vecs, axis=1)
        center_vec = axis_vecs[center_ind]

        # Remove center mol (it always make a 0 angle with itself)
        axis_vecs = axis_vecs[np.arange(len(axis_vecs)) != center_ind]

        # Now find what angle the vector on each mol makes with the center
        all_mags = np.linalg.norm(axis_vecs, axis=1) * np.linalg.norm(center_vec)
        all_mags += 1e-12      # Just to remove silly numerical errors in arccos
        all_dots = np.sum(axis_vecs * center_vec, axis=1)
        all_angles = np.arccos(all_dots / all_mags), axis_vecs

        return all_angles

    def plot_single(self, axis, bin_edges, counts, label=""):
        """
        Will plot 1 histogram on 1 axis according to bin_edges and counts.

        Inputs:
            * axes <plt.axis> OPTIONAL => The axis on which to plot.
        Outputs:
            (plt.bar, plt.axis) bar class and axis that the hist has been plotted on.
        """
        # Convert to degress
        #bin_edges = bin_edges
        #bars = axis.bar(bin_edges[:-1], counts, width=np.diff(bin_edges),
        #                label=label, alpha=0.7)
        min_bin, max_bin = min(bin_edges), max(bin_edges)
        axis.plot([-1, min_bin], [0, 0], 'k-')
        axis.plot([max_bin, 1], [0, 0], 'k-')
        line, = axis.plot(bin_edges[:-1], counts, 'k-')
        return line, axis

    def plot(self, axes=False, label=""):
        """
        Will plot a histogram of angles and return axis.

        If the axis argument is supplied then that axis will be used.

        Inputs:
            axes <array<plt.axis>> OPTIONAL => The axes on which to plot (len 2).
        Outputs:
            (plt.axis) Axes on which the histograms are plotted.
        """
        if axes is False:
            _, axes = plt.subplots(2)

        # Plot long axis angle distribution
        all_bars = []
        bars, ax = self.plot_single(axes[0], self.long_bin_edges,
                                    self.long_counts, label=label)
        all_bars.append(bars)
        axes[0] = ax

        # Plot short axis angle distribution
        all_bars = []
        bars, ax = self.plot_single(axes[1], self.short_bin_edges,
                                    self.short_counts, label=label)
        all_bars.append(bars)
        axes[1] = ax

        # Make it pretty
        axes[0].set_ylabel("Long Ax Density", fontsize=24)
        axes[1].set_ylabel("Short Ax Density", fontsize=26)
        axes[1].set_xlabel(r"cos($\theta$)", fontsize=26)

        if label:
            axes[0].legend(fontsize=16)
        for i, ax in enumerate(axes):
            ax.set_xlim([-1, 1])
        plt.tight_layout()

        return axes
