"""
   Helper functions that are used throughout the package are defined here
"""

import numpy as np

def molecule_leaflet_orientation(atompos1: np.array, atompos2: np.array,
      axis=np.array([0.0, 0.0, 1.0] ) ) -> int:
    ''' Takes two positions and returns 0 or 1
        if molecule is oriented up or downwards respectively
    '''
    new_coords = atompos1 - atompos2
    cos = np.dot(new_coords, axis) / np.linalg.norm(new_coords)
    return  0 if cos <= 0 else 1

def create_box_mirrors_2D(points: np.array, box: np.array) -> np.array:
    ''' A set of points is copied eightfold using dimensions in box '''
    outp_ar = []
    if not np.all(box != 0):
        raise ValueError("Box edge length must not be 0")
    img_vectors = [
        (0,  1), (1,  0), (1,  1), (0, -1),
        (-1, 0), (-1, -1), (1, -1), (-1, 1),
        (0, 0),
        ]
    for vec in img_vectors:
        newb = vec * box
        newp = points + newb
        outp_ar.append(newp)
    return np.array([p for ar in outp_ar for p in ar])

def angle_to_axis(vec: np.array, axis=np.array([0,0,1])) -> float:
    ''' Calculates angle of vector to axis (default: z axis)
        Returns angle in degree
    '''
    return np.arccos(np.dot(vec, axis)/np.linalg.norm(vec)) * (180/np.pi)

def dist_helper(v1, v2, box):
    img_vectors = [
        (0,  1),
        (1,  0),
        (1,  1),
        (0, -1),
        (-1, 0),
        (-1, -1),
        (1, -1),
        (-1, 1),
        (0, 0),
        ]
    ns = []
    for img in img_vectors:
        img = (*img, 0)
        newv = v2 + (box[:3]*img)
        ns.append(np.linalg.norm(v1-newv))
    return np.array(ns).min()

def rotate_2d(vec, angle_rad):
    ''' Rotate vec by angle_rad (anti clockwise) '''
    mat = np.matrix([[np.cos(angle_rad), np.sin(angle_rad) ], [-1*np.sin(angle_rad), np.cos(angle_rad) ]])
    nvec = np.array( mat.dot(vec) )[0]
    return nvec

def angle_clockwise(A, B):
    ''' Return the inner angle of two vectors '''
    inner_ang = np.dot(A, B) / np.linalg.norm(A) * np.linalg.norm(B)
    det = np.linalg.det(np.matrix([A, B]))
    if det < 0:
        return inner_ang
    else:
        return 2*np.pi - inner_ang
