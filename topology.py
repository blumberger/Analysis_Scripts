import numpy as np
import consts
import geom

import collections as coll


def atoms_to_mols(crds, num_ats_in_mol, cart_dims=3):
    """
    Will reshape array to divide up the atoms into arrays with molecules, i.e:

        reshape array from (num_ats, cart_dims) to (num_mols, ats_in_mol, cart_dims)

    Inputs:
        * crds => coordinates in shape (num_ats, cart_dims)
        * num_ats_in_mol => how many atoms in a single molecule

    Outputs:
        * array of shape (num_mols, ats_in_mol, cart_dims)
    """
    # Could test out that walrus operator here!
    if (len(crds) / num_ats_in_mol) != int(len(crds) / num_ats_in_mol):
        print("Please double check that you have entered the correct number of atoms per molecule!")
        raise SystemExit("The number of atoms per molecule does give an integer number of molecules!")

    num_mols = int(len(crds) // num_ats_in_mol)
    return np.reshape(crds, (num_mols, num_ats_in_mol, cart_dims))


def get_mol_COMs(at_crds, ats, num_ats_in_mol):
    """
    Will get the center of masses of the molecules within in the atomic
    coordinates array.

    Inputs:
        * at_crds => the atomic coordinates
        * ats => the atom types
    Outputs:
        * An array containing all the xyz coordinates of the COM of the molecule.
    """
    at_masses = {'C': consts.C_mass, 'H': consts.H_mass}

    ats_per_mol = atoms_to_mols(at_crds, num_ats_in_mol)
    at_types = atoms_to_mols(ats, num_ats_in_mol, 1)
    masses = [at_masses[i] for i in at_types[0, :, 0]]

    all_COM = [geom.get_COM(i, masses) for i in ats_per_mol]
    return all_COM


def get_long_axis_vec(molCrds, ats=[5, 23]):
    """
    Will get the vector describing the displacement from atom 6 to atom 24.
    In the geomety this was coded for these atoms sit on the long axis. If
    this isn't the case please change the ats list.

    Inputs:
        * molCrds => The coords of the atoms in the molecule.
        * ats => The indices of the atoms that sit on the long axis.

    Outputs:
        A 1D vector of length 3.
    """
    if len(ats) != 2:
        raise SystemExit("Only 2 input atoms are allowed")

    vec1 = np.array(molCrds[ats[0]])
    vec2 = np.array(molCrds[ats[1]])

    return vec1 - vec2


def get_long_axis_rotation_about(molCrds, vec=[1, 0, 0], long_ax_ats=[5, 23]):
    """
    Will get the rotation of the long axis of the molecule about a given vector.

    Inputs:
        * molCrds => 2D array of molecular coordinates of shape (num_atoms, 3)
        * vec => The vector to check the rotation about.
        * long_ax_ats => The indices of the atoms sitting on the long axis.
    
    Outputs:
        A float with the rotation angle of the long axis from the vector given.
    """
    rot_vec = get_long_axis_vec(molCrds, long_ax_ats)
    return geom.get_angle_between_2_vecs(rot_vec, vec)


def get_num_ats_per_mol(atCrds):
   """
   Will estimate the number of atoms per molecule by checking when the gradient with respect to atom number of the x coordinate spikes.

   Inputs:
      * atCrds => The atomic coordinates as an array of shape (num_at, 3)

   Outputs:
      An integer
   """
   x_crds = atCrds[:, 0]
   grad = np.gradient(x_crds)
   grad_tol = np.std(grad) * 2
   at_inds = np.arange(len(x_crds))[np.abs(grad) > grad_tol]

   counts = coll.Counter(np.diff(at_inds))
   counts = {i: counts[i] for i in counts if i > 5}

   return min(counts.keys()) + 1
