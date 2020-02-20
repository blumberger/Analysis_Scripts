#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module holds some functions to calculate certain geometric properties such
as angle between 2 vectors and center of masses etc...
"""

import numpy as np


#### Define Consts #####
PI_4_3 = 4./3. * np.pi #
PI_4   = 4. * np.pi    #
########################


def volume_sphere(r):
    """
    Returns the volume of a sphere
    """
    return PI_4_3 * r**3


def volume_concentric_spheres(R1, R2):
    """
    Returns the volume between 2 concentric spheres
    """
    if (R1 > R2):
        return PI_4_3 * (R1**3 - R2**3)
    else:
        return PI_4_3 * (R2**3 - R1**3)


def volume_differential_shell(r, dr):
    """
    Should be the same as volume_concentric_spheres for vanishing dr.
    """
    return PI_4 * (r**2) * dr


def beginning_tests(xyz, mass=False):
    """
    Carries out some basic tests at the beginning of a function to check
    if the input array looks right.
    """
    if mass and len(xyz) != len(mass):
        raise SystemExit("Mass array and xyz array not the same length")
    elif len(xyz[0]) != 3:
        raise SystemExit("This only works with 3D coords given in shape (num_crds, 3)")
    elif type(xyz) != type(np.array(1)):
        xyz = np.array(xyz)
        if mass:
            mass = np.array(mass)
    return xyz, mass


def get_COM(xyz, mass):
    """
    Will get the 3D center of mass of an array of points (xyz) given their
    masses (mass).

    Inputs:
        * xyz   => 2D array of atomic coords of shape <num coords, 3>
        * mass  => 1D array of masses same as len(xyz).

    Returns:
        * 3D array with [x,y,z] of center of mass
    """
    xyz, mass = beginning_tests(xyz, mass)

    tot_mass = sum(mass)

    xmean = sum(xyz[:,0]*mass) / tot_mass
    ymean = sum(xyz[:,1]*mass) / tot_mass
    zmean = sum(xyz[:,2]*mass) / tot_mass

    return [xmean, ymean, zmean]


def find_center_atom(crds):
    """
    Will find the atom that is closest to the arthimetic center

    Inputs:
        * 2D crds array of shape (num_atoms, 3)

    Outpus:
        * index of central most atom
        * xyz of central most atom
    """
    crds, _ = beginning_tests(crds)

    # Center point = (max - min) / 2
    center_point = (np.max(crds, axis=0) - np.min(crds, axis=0))/2.

    dist_from_center = np.linalg.norm(crds - center_point, axis=1)
    min_at_ind = np.argmin(dist_from_center)
    center_xyz = crds[min_at_ind]

    return min_at_ind, center_xyz, dist_from_center


def get_vector_from_cds(at1, at2):
    """
    Will get the displacement vector from at1 to at2.

    Inputs:
        * at1 => 1D vector of length 3.
        * at2 => 1D vector of length 3.

    Outputs:
        * 3D vector pointing from atom1 to atom2.
    """
    at1, at2 = np.array(at1), np.array(at2)
    return at2 - at1


def get_angle_between_2_vecs(vec1, vec2):
    """
    Will find the angle between 2 vectors.

    Inputs:
        * vec1 => 1D vector (np.array)
        * vec2 => 1D vector (np.array)

    Outputs:
        Single float with angle between the 2 vectors
    """
    if len(vec1) != len(vec2):
        raise SystemExit("Vectors must be same dimension!")

    dot_prod = np.dot(vec1, vec2)
    norm_mult = np.linalg.norm(vec1) * np.linalg.norm(vec2)
    return np.arccos(dot_prod / norm_mult)
