import load_xyz
import geom
import topology
import consts

import numpy as np
import matplotlib.pyplot as plt
import os


folder = "/home/kev-boi-2k18/Documents/PhD/AMORPHOUS_PENTACENE/10ns_quenched_NVT_equilibriation"
filename = "10nsNPT_quench.xyz"

xyz_fpath = os.path.join(folder, filename)
if not os.path.isfile(xyz_fpath):
    raise SystemExit("Can't find file %s" % xyz_fpath)


def calc_RDF(crds, max_dist=False, dr=False, origin=False):
    """
    Will calculate the radial distribution function of a coordinate array with
    origin at origin.

    Inputs:
        * crds = 2D array of shape (num_atoms, 3)
        * origin = 3 element array with xyz of origin
        * dr = with used for the concentric spheres.
    Output:
        * an array with the radial distribution as a function of radius
          from the origin.
    """
    if not origin:
        at_ind, origin = geom.find_center_atom(crds)

    if np.shape(origin) != (3,):
        raise SystemExit("Shape of the origin needs to be (3,)")
    
    dist_from_origin = np.linalg.norm(crds - origin, axis=1)

    if not max_dist:
        max_dist = max(dist_from_origin)
    if not dr:
        dr = max_dist / 100.    

    spaces = np.arange(0, max_dist, dr)
    rdf = [0]
    r = [0, ]
    for i in range(1, len(spaces) - 1):
        start_R, end_R = spaces[i], spaces[i+1]
        
        num_atoms = sum((end_R > dist_from_origin) & (dist_from_origin >= start_R))
        volume = geom.volume_concentric_spheres(start_R, end_R)

        rdf.append(num_atoms / volume)
        r.append((start_R + end_R) / 2)

    return r, rdf


def plot_RDF_from_file(xyz_file, num_ats_per_mol=36, max_dist=False, dr=False, origin=False):
    """
    Will load the xyz_coords from an xyz file and will plot the RDF vs R
    for the COM of each molecule in the file.
    """
    # Load coords
    if not os.path.isfile(xyz_file):
        raise SystemExit("Can't find file %s" % xyz_file)

    ats, crds = load_xyz.read_1_step_xyz(xyz_file)
    mol_COMs = topology.get_mol_COMs(crds, ats, num_ats_per_mol)
    
    r, rdf = calc_RDF(mol_COMs, dr=dr, max_dist=max_dist, origin=origin)

    f, a = plt.subplots()
    a.plot(r, rdf, 'k--')
    a.set_xlabel(r"R [$\AA$]", fontsize=18)
    a.set_ylabel(r"g(r)", fontsize=18)
    a.set_title("File = %s" % xyz_file)
    plt.show()
    
    

    
    
        
