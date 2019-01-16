## /notebooks

This contains the Jupyter notebook files, and converted Python scripts, used for the analysis of the neutron reflectometry data and molecular dynamics simulations.

When the `Snakefile` in the top level directory is run `nbconvert` generates the `.py` files from the existing `.ipynb` notebooks. These Python scripts are then run for each surface pressure and the molecular dynamics simulations analysed.

#### Contents

- `simulation/analysis.*` - analyses the molecular dynamics simulations to produce the resultant neutron reflectometry
- `simulation/chain_tilt.*` - calculates some parameters of interest from the molecular dynamics simulation
- `simulation/density_plot.*` - plots the number density for the slipids simulation at 30 mN/m
- `simulation/martiniorder.*` - plots the scattering length density for the martini simulation at different layer thicknesses
- `simulation/plot.*` - plots the reflectometry from simulation and scattering length density profiles
- `traditional/analysis.*` - analyses the neutron reflectometry using the chemically-consistent analysis method
- `traditional/plot.*` - plots the reflectometry from traditional analyses and scattering length density profiles
