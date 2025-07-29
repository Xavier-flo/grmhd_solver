"""
src/grid.py

Generate a 2D polar grid (r, φ) for the GRMHD solver.
"""

import numpy as np

def create_polar_grid(r_min, r_max, nr, nphi, endpoint=True):
    """
    Create arrays of r and φ coordinates on a 2D polar mesh.

    Parameters
    ----------
    r_min : float
        Inner radius of the grid.
    r_max : float
        Outer radius of the grid.
    nr : int
        Number of radial zones.
    nphi : int
        Number of azimuthal zones.
    endpoint : bool
        If True, include r_max and 2π in the grids.

    Returns
    -------
    R : ndarray, shape (nr, nphi)
        Radial coordinate at each cell center.
    Phi : ndarray, shape (nr, nphi)
        Azimuthal coordinate at each cell center.
    dr : float
        Radial zone width.
    dphi : float
        Azimuthal zone width.
    """
    r = np.linspace(r_min, r_max, nr, endpoint=endpoint)
    phi = np.linspace(0, 2*np.pi, nphi, endpoint=endpoint)
    R, Phi = np.meshgrid(r, phi, indexing='ij')
    dr = (r_max - r_min) / (nr - (1 if endpoint else 0))
    dphi = (2*np.pi) / (nphi - (1 if endpoint else 0))
    return R, Phi, dr, dphi
