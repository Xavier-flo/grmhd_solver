"""
src/metric.py

Schwarzschild metric components in the equatorial plane (θ = π/2),
using geometric units G = c = 1.
Coordinates: (t, r, φ).
"""

import numpy as np

def schwarzschild_metric(r, M=1.0):
    """
    Compute diagonal metric components at radius r.
    
    Parameters
    ----------
    r : float or array-like
        Radial coordinate(s).
    M : float
        Black hole mass (default 1.0).
    
    Returns
    -------
    g : dict
        {'g_tt': g_{tt}, 'g_rr': g_{rr}, 'g_phiphi': g_{φφ}'}
    """
    r = np.asarray(r)
    f = 1.0 - 2.0*M / r
    return {
        'g_tt':    -f,
        'g_rr':     1.0 / f,
        'g_phiphi': r**2
    }
