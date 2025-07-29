# 2D Ideal GRMHD Accretion‑Disk Solver

## Overview
A toy GRMHD code solving the ideal MHD equations in a Schwarzschild metric on a 2D polar grid.

## Setup
\`\`\`bash
conda env create -f environment.yml
conda activate grmhd
\`\`\`

## Running
- **Local test:** \`python scripts/run_local.py --grid-size 32 32 --timesteps 100\`
- **HPC:** \`sbatch scripts/submit.slurm\`

## Structure
- \`src/\`: physics modules  
- \`notebooks/\`: interactive demos  
- \`scripts/\`: launch scripts  
- \`tests/\`: unit tests  
- \`docs/\`: design notes  
