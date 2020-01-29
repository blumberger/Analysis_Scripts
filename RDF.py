import load_xyz
import geom
import topology
import consts

import numpy as np
import matplotlib.pyplot as plt
import os


#folder = "/home/kev-boi-2k18/Documents/PhD/AMORPHOUS_PENTACENE/10ns_quenched_NVT_equilibriation"`
#filename = "10nsNPT_quench.xyz"
#
#xyz_fpath = os.path.join(folder, filename)
#if not os.path.isfile(xyz_fpath):
#    raise SystemExit("Can't find file %s" % xyz_fpath)


def calc_RDF(crds, max_dist=False, dr=False, origin=False, nbins="2sqrt(N)"):
    """
    Will calculate the radial distribution function of a coordinate array with
    origin at origin.

    Inputs:
        * crds = 2D array of shape (num_atoms, 3)
        * origin = 3 element array with xyz of origin
        * dr = with used for the concentric spheres.
        * nbins = how many divisions of the total volume to use in calculating g(r)
    Output:
        * an array with the radial distribution as a function of radius
          from the origin.
    """
    if nbins and dr:
        print("Warning ignoring parameter nbins as dr has been set")
    if not origin:
        at_ind, origin = geom.find_center_atom(crds)

    if np.shape(origin) != (3,):
        raise SystemExit("Shape of the origin needs to be (3,)")
    
    dist_from_origin = np.linalg.norm(crds - origin, axis=1)

    if not max_dist:
        max_dist = max(dist_from_origin) * 0.7

    # Tot num of particles within a sphere
    N = float(sum(dist_from_origin <= max_dist))
    # Number density of particles within a sphere
    rho = N / geom.volume_sphere(max_dist)
    if not dr and nbins == "2sqrt(N)":
        dr = max_dist / (2*np.sqrt(N))
    elif type(nbins) == int:
        dr = max_dist / float(nbins)
    elif type(nbins) == float:
        dr = max_dist / float(int(nbins))
        
    spaces = np.arange(0, max_dist, dr)
    rdf = [0]
    r = [0, ]
    for i in range(1, len(spaces) - 1):
        start_R, end_R = spaces[i], spaces[i+1]
        
        num_atoms = sum((end_R >= dist_from_origin) & (dist_from_origin > start_R))
        volume = geom.volume_differential_shell(start_R, dr)
        norm_const = volume / rho

        rdf.append(num_atoms / (volume*rho))
        r.append((start_R + end_R) / 2)

    return r, rdf


def plot_RDF_from_file(xyz_file, num_ats_per_mol=False, max_dist=False, dr=False,
                       origin=False, nbins="2sqrt(N)"):
    """
    Will load the xyz_coords from an xyz file and will plot the RDF vs R
    for the COM of each molecule in the file.
    """
    # Load coords
    if not os.path.isfile(xyz_file):
        raise SystemExit("Can't find file %s" % xyz_file)


    ats, crds = load_xyz.read_1_step_xyz(xyz_file)
    if not num_ats_per_mol: num_ats_per_mol = topology.get_num_ats_per_mol(crds)

    mol_COMs = topology.get_mol_COMs(crds, ats, num_ats_per_mol)
    
    r, rdf = calc_RDF(mol_COMs, dr=dr, max_dist=max_dist, origin=origin, nbins=nbins)

    f, a = plt.subplots()
    a.plot(r, rdf, 'k')
    a.set_xlabel(r"R [$\AA$]", fontsize=20)
    a.set_ylabel(r"g(r)", fontsize=20)

    title = xyz_file[xyz_file.rfind('/')+1:]
    title = title.replace("_", " ").replace(".xyz", "") 
    title = title.replace("quenched", "quenching").replace(" quench ", " quenching ").replace(" quench.", " quenching ")
    title = title.replace("NPT", "").replace("NVT", "") 

    a.set_title(title, fontsize=22)
    plt.tight_layout()

    return f, a, title
    
