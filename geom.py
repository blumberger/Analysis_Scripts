import numpy as np


#### Define Consts #####
pi_4_3 = 4./3. * np.pi #
########################





def volume_sphere(r):
    """
    Returns the volume of a sphere
    """
    return pi_4_3 * r**3


def volume_concentric_spheres(R1, R2):
    """
    Returns the volume between 2 concentric spheres
    """
    if (R1 > R2):
        return pi_4_3 * (R1**3 - R2**3)
    else:
        return pi_4_3 * (R2**3 - R1**3)


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

    center_point = np.mean(crds, axis=0)
    dist_from_mean = np.linalg.norm(crds - center_point, axis=1)    
    min_at_ind = np.argmin(dist_from_mean)
    center_xyz = crds[min_at_ind]

    return min_at_ind, center_xyz
