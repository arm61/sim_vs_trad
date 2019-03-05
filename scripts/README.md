## /scripts

This contains the Python scripts, used for the analysis of the neutron reflectometry data and molecular dynamics simulations.

#### Contents

- `simulation/martini_order.py` - plotting the difference between the MARTINI simulation scattering length density under different conditions
- `simulation/md_simulation.py` - a refnx class to find the reflectometry and SLD profile directly from a molecular dynamics simulation
- `simulation/ref_help.py` - some reflectometry helper scripts
- `simulation/roughness.py` - quantification of roughness from MD simulation
- `simulation/short_sim_analysis.py` - the implementation of the md_simulation class mentioned above for the specific work but only for the first 5 ns
- `simulation/sim_analysis.py` - the implementation of the md_simulation class mentioned above for the specific work
- `simulation/sim_lengths.py` - functions to get the scattering length for the md_simulation
- `simulation/tail_length.py` - determine the length of the carbon tails
- `simulation/water_plot.py` - plot the intrinsic surface approach for water penetration
- `simulation/waters.py` - implement the intrinsic surface approach to find water penetration
- `simulation/wph.py` - calculate the number of waters per head group
- `traditional/cc_plot.py` - plot the chemically-consistent model
- `traditional/chemically_consistent.py` - implementation of the mol_vol class
- `traditional/mol_vol.py` - a refnx class to implement the chemically-consistent model
- `traditional/ref_help.py` - some reflectometry helper scripts
