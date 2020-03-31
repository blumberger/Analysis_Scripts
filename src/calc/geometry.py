#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contains various functions to perform/calculate geometric operations/properties.

These include things like: calculate rotation angle, get rotation matricies
"""

import numpy as np


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

