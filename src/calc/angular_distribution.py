#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to calculate the angular distributions of molecules in a system.
"""
import numpy as np
import json
import pandas as pd
import matplotlib.pyplot as plt

# Import calculating functions
from src.calc import general_calc as gen_calc
from src.calc import geometry as geom
from src.calc import molecule_utils as mol_utils

from src.data import consts

# import type checking functions
from src.system import type_checking as type_check

class Angular_Dist(gen_calc.Calc_Type):
    """
    Will calculate the angular distributions of the molecules in a system.

    Inputs:
        * Variable <Variable> => An instance of the Variable class

    Important Attributes:
        * required_metadata <tuple> => Any keys that are required in the metadata dictionary.
        * required_calc <tuple> => Any values that need calculating to calculate this value.
        * data <*> => The data that has been calculated.
    """
    _write_types = ('json', 'csv',)
    required_metadata = ('atoms_per_molecule', 'plot_angular_distribution')
    _defaults = {'number_bins': 'auto', 'histogram_density': True,
                 'plot_angular_distribution': False, 'nearest_neighbour_dist': float("inf"),
                }
    # Need these 3 attributes to create a new variable type
    data = {}
    metadata = {'file_type': 'json'}
    required_data_types = ('pos',)
    name = "Angular Distribution"
    with open(consts.PT_FILEPATH) as f: PT = json.load(f)

    def _calc_(self):
        """
        Will calculate the angular distribution of the molecular system.

        This will loop over each molecule, find the rotation (wrt long ax, short
        ax of central molecule) of the long and short axes of the molecule then
        create a histogram of this data.
        """
        self.Var['coordinate_wrapping'] = 'unwrapped'
        self.xyz_data = self.Var.data.get_xyz_data()

        ats_per_mol = self.Var.metadata['atoms_per_molecule']
        if 'long_axis_atoms' not in self.Var.metadata:
            long_ax_ats = self.metadata['long_axis_atoms']
        else: long_ax_ats = self.Var.metadata['long_axis_atoms']
        if 'short_axis_atoms' not in self.Var.metadata:
            short_ax_ats = self.metadata['short_axis_atoms']
        else: short_ax_ats = self.Var.metadata['short_axis_atoms']

        # Divide atomic coordinates into molecular coordinates


        # This line means we will only work with the first file that has xyz data.
        # This will need changing to make it more general.
        all_at_crds = self.xyz_data[0]

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
            self.long_angles = []
            self.short_angles = []
            count = 0
            centroids = np.mean(mol_crds, axis=1)
            inds = np.arange(nmol)

            for mol1 in range(nmol):
                disp = centroids - centroids[mol1]
                dist = np.linalg.norm(disp, axis=1)
                close_inds = inds[dist < self.metadata['nearest_neighbour_dist']]

                long_angles = self.get_angles_between_vecs(self.long_vecs[mol1], self.long_vecs[close_inds],
                                                           long_mags[mol1], long_mags[close_inds])
                short_angles = self.get_angles_between_vecs(self.short_vecs[mol1], self.short_vecs[close_inds],
                                                            short_mags[mol1], short_mags[close_inds])

                self.short_angles.extend(short_angles)
                self.long_angles.extend(long_angles)
                count += len(close_inds)

            self.long_counts, self.long_bin_edges = np.histogram(self.long_angles,
                                  density=True, bins=self.metadata['number_bins'])
            self.short_counts, self.short_bin_edges = np.histogram(self.short_angles,
                                  density=True, bins=self.metadata['number_bins'])

        if self.metadata['plot_angular_distribution'] is True:
            self.plot()
            plt.show()

    def get_angles_between_vecs(self, vec, vecs, mag, mags):
        """
        Will get the angle a vector makes with many other vectors

        Inputs:
            * vec_ind <int> => The index of the vector to compare with
            * vecs <array> => All vectors
            * mags <array> => The magnitude of all the vectors
        Outputs:
            <array> All angles between vecA and other_vecs
        """
        dots = np.sum(vec * vecs, axis=1)
        mags = mag * mags
        angles = dots/mags #) * 180/np.pi

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

        # The mean averages over multiple pairs of atoms.
        #  i.e. if [[3, 5], [10, 13]] was given then we would get the
        #       average the vector from 3 -> 5 and 10 -> 13.
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

    def get_csv_data(self):
        """
        Will return the data in a form the csv writer can use.

        This is only used for writing.
        """
        short_len = len(self.short_counts)
        long_len = len(self.long_counts)
        return pd.DataFrame({
                'short_ax_count': self.short_counts.tolist(),
                'short_ax_edges': self.short_bin_edges.tolist()[:short_len],
                'long_ax_count': self.long_counts.tolist(),
                'long_ax_edges': self.long_bin_edges.tolist()[:long_len],
                })

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
        line, = axis.plot(bin_edges[:-1], counts, '-', label=label)

        if min_bin > -1: axis.plot([-1, min_bin], [0, 0], color=line.get_color())
        if max_bin < 1:  axis.plot([max_bin, 1], [0, 0], color=line.get_color())

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
            axes[1].legend(fontsize=16)

        plt.legend()

        #for i, ax in enumerate(axes):
        #    ax.set_xlim([-1.05, 1.05])
        plt.tight_layout()

        return axes
