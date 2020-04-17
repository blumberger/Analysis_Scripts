#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contains various functions to perform/calculate geometric operations/properties.

These include things like: calculate rotation angle, get rotation matricies
"""

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import DBSCAN



def rotate_crds(xyz, rot_mat):
    """
    Will rotate some coordinates using numpy.apply_along_axis.

    Inputs:
        * xyz <array> => An array that is of shape that has 1 axis with a length that is the same
                         as the length of the rotation matrix. This will be the coord axis.
        * rot_mat <array> => The rotation matrix to use to rotate the coordinates
    Outputs:
        <array> xyz coordinates of rotated system.
    """
    crds_dim = len(rot_mat)

    # First find the coordinate axis
    coord_axis = _get_crd_ax(xyz, crds_dim)
    rotate_vec = lambda vec: np.matmul(rot_mat, vec)

    return np.apply_along_axis(rotate_vec, coord_axis, xyz)

def angle_between_vecs(vec1, vec2):
    """
    Will find the angle between 2 vectors.

    This uses np.arccos( (a . b) / |a||b| )

    Inputs:
      * vec1 <list> => The first vector
      * vec2 <list> => The second vector
    Outputs:
      <float> angle.
    """
    dot_prod = np.dot(vec1, vec2)
    mags = np.linalg.norm(vec1) * np.linalg.norm(vec2)

    return np.arccos(dot_prod / mags)

def map_vec1_to_unit(vec, unit=[1.0, 0.0, 0.0]): 
    """
    Will map 1 vector to another vector and return the rotation matrix.

    This works by finding the perp vector to the input vec and unit (using np.cross).
    The angle is then found between the input vec and unit.
    The rotation matrix is then constructed to rotate about the calculated axis by
     the calculated angle.

    Inputs: 
        * vec <list> => The first vector to rotate
        * unit <list> => The vector to map 'vec' onto
    Outputs:
        <array<array>> The rotation matrix.

    Taken from: https://stackoverflow.com/questions/43507491/imprecision-with-rotation-matrix-to-align-a-vector-to-an-axis
    """
    # Normalize vector length 
    vec /= np.linalg.norm(vec) 
    unit /= np.linalg.norm(unit)

    # Get axis 
    uvw = np.cross(vec, unit) 

    # compute trig values - no need to go through arccos and back 
    rcos = np.dot(vec, unit) 
    rsin = np.linalg.norm(uvw) 

    #normalize and unpack axis 
    if not np.isclose(rsin, 0): 
       uvw /= rsin 
    u, v, w = uvw 

    # Compute rotation matrix - re-expressed to show structure 
    return ( 
       rcos * np.eye(3) + 
       rsin * np.array([ 
           [ 0, -w,  v], 
           [ w,  0, -u], 
           [-v,  u,  0] 
       ]) + 
       (1.0 - rcos) * uvw[:,None] * uvw[None,:] 
    )

def get_system_size_info(xyz):
    """
    Will get some basic geometric info about the inputted system.

    Inputs:
        * xyz <array> => The array containing xyz data
    Outputs:
        <dict> A dictionary containing xmin, xmax, ymin, ymax, zmin, zmax, xlen, ylen, zlen.
    """
    sys_info = {}

    # Get which axis contains coordinate info
    coord_axis = _get_crd_ax(xyz, 3)

    # Get the coordinate info
    x, y, z = np.take(xyz, 0, coord_axis), np.take(xyz, 1, coord_axis), np.take(xyz, 2, coord_axis)

    # Calculate some basic properties
    sys_info['xmin'], sys_info['xmax'] = np.min(x), np.max(x)
    sys_info['ymin'], sys_info['ymax'] = np.min(y), np.max(y)
    sys_info['zmin'], sys_info['zmax'] = np.min(z), np.max(z)

    sys_info['xlen'] = sys_info['xmax'] - sys_info['xmin']
    sys_info['ylen'] = sys_info['ymax'] - sys_info['ymin']
    sys_info['zlen'] = sys_info['zmax'] - sys_info['zmin']

    return sys_info

def _get_crd_ax(xyz, dim=3):
    """
    Will get which axis is the coordinate axis

    Inputs:
        * xyz <array> => The data to find the coordinate axis for
        * dim <int>   => The number of elements in the coordinate axis
    Outputs:
        <int> The coordinate axis
    """
    coord_axis_lens = np.shape(xyz)

    if dim not in coord_axis_lens:
        print("The xyz coordinates can't be found in the array you've inputted.")
        print(f"The inputted array has shape {coord_axis_lens}.\n\n")
        raise SystemExit("Can't rotate coordinates.")
    elif coord_axis_lens.count(dim) > 1:
        print("Too many possible axes that could be coordinates.")
        print(f"The inputted array has shape {coord_axis_lens}.\n\n")
        raise SystemExit("Can't rotate coordinates.")

    return coord_axis_lens.index(dim)


def cluster_1D_points(points_1D, eps=0.02, min_samples=1):
    """
    Will return an array of clustered 1D points.

    This uses sklearn's implementation of DBSCAN to cluster points.

    Inputs:
        * points_1D <array> => 1D array of data
        * eps <float> OPTIONAL => Space between to define clusters (default 0.02)
        * min_samples <int> OPTIONAL => min num of points in cluster (default 1)
    Outputs
        <DBSCAN class>, <list> The DBSCAN class and the avg point in each cluster
    """
    if len(points_1D) == 0:
        raise SystemError("1D Point clustering error, input list has length 0.")
    if type(points_1D) == list:
        points_1D = np.array(points_1D)
    if isinstance(points_1D, (pd.DataFrame, pd.Series)) :
        points_1D = points_1D.to_numpy()

    scaled_data = StandardScaler().fit_transform(points_1D.reshape(-1, 1))
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(scaled_data)
    clustered_means = [np.mean(points_1D[db.labels_ == i]) for i in set(db.labels_)]
    return db, clustered_means

def find_local_min(y, start_ind, allow_boundaries=False):
    """
    Apply a steepest descent like algorithm to get a local minima.

    This is a very quick and dirty algorithm that steps by one in the index of
    values until the value of y stops decreasing.

    This isn't built for speed -if a faster algorithm is required look at conjugate
    gradient methods in scipy.

    Inputs:
        * x <arr> => 1D array of x positions
        * y <arr> => 1D array of y positions
        * start_ind <int> => where to start looking for the local minima
        * allow_boundaries <bool> => Whether to allow the ends/boundaries of the data or 
                                     only 'true' minima.

    Outputs:
        <ind>, <float> The index of the local minima and its value.
    """
    # Check we are looking at indices within the range of the data.
    index = np.arange(len(y))
    if not (start_ind in index and start_ind - 1 in index and start_ind + 1 in index):
        return False, (False, False)

    # Find whether to go left or right
    cond = True
    if y[start_ind - 1] < y[start_ind]:  step = -1
    elif y[start_ind + 1] < y[start_ind]:  step = +1
    else: cond = False

    # Keep stepping left or right until the value stops decreasing
    while cond:
        start_ind += step
        if not (start_ind in index and start_ind + step in index):
            break

        if y[start_ind + step] >= y[start_ind]:
            cond = False

    # if we hit a boundary (not a 'proper' minima) then return False
    if cond and allow_boundaries is False:
        return False, (False, False)

    return start_ind, y[start_ind]