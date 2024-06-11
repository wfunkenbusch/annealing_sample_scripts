# annealing_sample_scripts

Sample scripts for running simulations for paper "Shear annealing of a self-interacting sheet" (submitted June 2024). All simulations were run using HOOMD-blue 2.4 and Python 3.6.9 using a custom package from Silmore, K. S., Strano, M. S., & Swan, J. W. (2021). *Soft Matter* which was adapted from Fiore, A. M., *et al.* (2017). *The Journal of Chemical Physics*.

Scripts:

generate_hexagonal_sheet.jl - Julia file for generating a flat sheet in the flow-vorticity plane using a custom package from Silmore, K. S., Strano, M. S., & Swan, J. W. (2021). *Soft Matter* which can be found at https://github.com/ksil/KMesh.jl. Tested with Julia 1.1 and gsd 2.9.0.

annealing_simulation.py - sample annealing simulation with $2000\dot{\gamma}_0t$ for each constant shear portion.

protocol_simulation.py - sample protocol simulation with $2000$\dot{\gamma}t$ for each constant shear portion.
