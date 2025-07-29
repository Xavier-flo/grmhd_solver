"""
src/flux.py

Compute HLL‐type numerical fluxes for a 1D hyperbolic system.
"""

import numpy as np

def hll_flux(uL, uR, fL, fR, sL, sR):
    """
    HLL flux between left and right states.

    Parameters
    ----------
    uL, uR : array-like
        Conserved variable vectors on the left and right.
    fL, fR : array-like
        Physical flux vectors at left and right.
    sL, sR : float
        Minimum and maximum signal speeds (estimates).

    Returns
    -------
    flux : ndarray
        HLL‐approximate flux at the interface.
    """
    uL = np.asarray(uL)
    uR = np.asarray(uR)
    fL = np.asarray(fL)
    fR = np.asarray(fR)

    # If all waves move right
    if sL >= 0:
        return fL
    # If all waves move left
    if sR <= 0:
        return fR

    # HLL formula
    return (sR*fL - sL*fR + sL*sR*(uR - uL)) / (sR - sL)
