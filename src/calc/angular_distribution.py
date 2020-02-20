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
    required_metadata = ('long_axis_atoms', 'short_axis_atoms',
                         'atoms_per_molecule')
    _defaults = {'number_bins': 'auto', 'histogram_density': True}
    # Need these 3 attributes to create a new variable type
    data = {}
    metadata = {'file_type': 'json'}
    name = "Angular Distribution"
    with open(consts.PT_FILEPATH) as f: PT = json.load(f)

    def get_data(self):
        """
        Will get the data to use from the inputted class.
        """
        if 'xyz_data' in dir(self.Var.data):
            self.data = self.Var.data.xyz_data
        elif 'csv_data' in dir(self.Var.data):
            self.data = self.Var.data.csv_data[['x', 'y', 'z']].to_numpy()
            self.data = np.array([self.data])

    def calc(self):
        """
        Will calculate the angular distribution of the molecular system.

        This will loop over each molecule, find the rotation (wrt long ax, short
        ax of central molecule) of the long and short axes of the molecule then
        create a histogram of this data.
        """
        self.get_data()
        ats_per_mol = self.Var.metadata['atoms_per_molecule']
        long_ax_ats = self.Var.metadata['long_axis_atoms']
        short_ax_ats = self.Var.metadata['short_axis_atoms']
        all_at_crds = self.data

        self.long_ax_angles, self.long_ax_vecs = [], []
        self.long_ax_bin_edges, self.long_ax_counts = [], []
        self.short_ax_angles, self.short_ax_vecs = [], []
        self.short_ax_bin_edges, self.short_ax_counts = [], []
        for at_crds in all_at_crds:
            # Divide atomic coordinates into molecular coordinates
            mol_crds = mol_utils.atoms_to_mols(at_crds, ats_per_mol,
                                                    nstep=len(at_crds))

            # Get the long ax angles
            _ = self.get_angle_dist(mol_crds, long_ax_ats)
            self.long_ax_angles.append(_[0])
            self.long_ax_vecs.append(_[1])
            hist = np.histogram(_[0], bins=self.metadata['number_bins'],
                               density=self.metadata['histogram_density'])
            self.long_ax_counts.append(hist[0])
            self.long_ax_bin_edges.append(hist[1])

            # Get the long ax angles
            _ = self.get_angle_dist(mol_crds, short_ax_ats)
            self.short_ax_angles.append(_[0])
            self.short_ax_vecs.append(_[1])
            hist = np.histogram(_[0], bins=self.metadata['number_bins'],
                               density=self.metadata['histogram_density'])
            self.short_ax_counts.append(hist[0])
            self.short_ax_bin_edges.append(hist[1])

        self.short_ax_angles = np.array(self.short_ax_angles)
        self.short_ax_vecs = np.array(self.short_ax_vecs)
        self.short_ax_bin_edges = np.array(self.short_ax_bin_edges)
        self.short_ax_counts = np.array(self.short_ax_counts)

        self.long_ax_angles = np.array(self.long_ax_angles)
        self.long_ax_vecs = np.array(self.long_ax_vecs)
        self.long_ax_bin_edges = np.array(self.long_ax_bin_edges)
        self.long_ax_counts = np.array(self.long_ax_counts)


    def json_data(self):
        """
        Will return data in a form that the json writer can write.

        This function is just used for writing
        """
        data = {
                'short_ax_angles': self.short_ax_angles.tolist(),
                'short_ax_vecs': self.short_ax_vecs.tolist(),
                'short_ax_histogram': {
                                        'counts': self.short_ax_counts.tolist(),
                                        'bin_edges': self.short_ax_bin_edges.tolist(),
                                      },
                'long_ax_angles': self.long_ax_angles.tolist(),
                'long_ax_vecs': self.long_ax_vecs.tolist(),
                'long_ax_histogram':  {
                                        'counts': self.long_ax_counts.tolist(),
                                        'bin_edges': self.long_ax_bin_edges.tolist(),
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
        center_ind, center = geom.find_center_atom(avg_mol_crds)

        at1, at2 = at_inds[0], at_inds[1]

        # First get the vector describing an axis for all mols and the center one
        axis_vecs = mol_crds[:, at1] - mol_crds[:, at2]
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
        edges = bin_edges * 180 / np.pi

        bars = axis.bar(edges[:-1], counts, width=np.diff(edges), label=label)
        return bars, axis

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
        bars, ax = self.plot_single(axes[0], self.long_ax_bin_edges[0],
                                    self.long_ax_counts[0], label=label)
        all_bars.append(bars)
        axes[0] = ax

        # Plot short axis angle distribution
        all_bars = []
        bars, ax = self.plot_single(axes[1], self.short_ax_bin_edges[0],
                                    self.short_ax_counts[0], label=label)
        all_bars.append(bars)
        axes[1] = ax

        # Make it pretty
        axes[0].set_ylabel("Long Ax Density")
        axes[1].set_ylabel("Short Ax Density")
        axes[1].set_xlabel(r"Angle [$^o$]")

        for ax in axes:
            if label:
                ax.legend()
            ax.set_xlim([0, 180])

        return axes
