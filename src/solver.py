"""
src/solver.py

Main time‐stepping driver for the 2D ideal GRMHD toy solver.
"""

import numpy as np
import h5py

from src.metric import schwarzschild_metric
from src.grid import create_polar_grid
from src.flux import hll_flux

def initialize_U(R, Phi):
    """
    Set up initial conserved variables U on the grid.
    Here U = [rho, rho*v_r, rho*v_phi, B_r, B_phi] for a simple test.
    """
    nr, nphi = R.shape
    U = np.zeros((nr, nphi, 5))
    # Example: uniform density with a small radial velocity perturbation
    U[..., 0] = 1.0                     # rho
    U[..., 1] = 0.01 * np.sin(Phi)      # rho*v_r
    U[..., 2] = 0.0                     # rho*v_phi
    U[..., 3] = 0.001                   # B_r
    U[..., 4] = 0.0                     # B_phi
    return U

def step(U, R, dr, dphi, dt):
    """
    Single time‐step forward (1st order Euler).
    Compute fluxes in r and phi and update U.
    """
    nr, nphi, nvar = U.shape
    U_new = U.copy()
    
    # Loop over each cell interface in r‐direction
    for i in range(1, nr):
        for j in range(nphi):
            # Left/Right states
            uL = U[i-1,j];  uR = U[i,j]
            # Compute local fluxes f(U)
            fL = uL  # placeholder; replace with physical flux function
            fR = uR
            # Estimate signal speeds (e.g. max eigenvalues)
            sL, sR = -1.0, +1.0
            F = hll_flux(uL, uR, fL, fR, sL, sR)
            U_new[i-1,j] -= dt/dr * F
            U_new[i,j]   += dt/dr * F
    
    # Similar loop in phi‐direction
    for i in range(nr):
        for j in range(nphi-1):
            uL = U[i,j];  uR = U[i,j+1]
            fL, fR = uL, uR
            sL, sR = -1.0, +1.0
            F = hll_flux(uL, uR, fL, fR, sL, sR)
            U_new[i,j]   -= dt/dphi * F
            U_new[i,j+1] += dt/dphi * F
    
    return U_new

def run_simulation(r_min, r_max, nr, nphi, t_end, cfl=0.4, output_file='out.h5'):
    """
    Full simulation driver.
    """
    # Build grid
    R, Phi, dr, dphi = create_polar_grid(r_min, r_max, nr, nphi)
    # Initialize U
    U = initialize_U(R, Phi)
    t = 0.0; step_num = 0
    
    # HDF5 output
    f5 = h5py.File(output_file, 'w')
    dset = f5.create_dataset('U0', data=U)
    
    while t < t_end:
        # Simple CFL time‐step estimate
        dt = cfl * min(dr, dphi)
        if t + dt > t_end:
            dt = t_end - t
        
        U = step(U, R, dr, dphi, dt)
        t += dt; step_num += 1
        
        # Save every 10 steps
        if step_num % 10 == 0:
            dset = f5.create_dataset(f'U{step_num}', data=U)
    
    f5.close()

if __name__ == '__main__':
    # Example run parameters
    run_simulation(r_min=2.5, r_max=20.0, nr=64, nphi=128, t_end=1.0)
