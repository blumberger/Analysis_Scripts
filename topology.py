import numpy as np
import consts
import geom


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
