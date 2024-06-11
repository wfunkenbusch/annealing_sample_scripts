from hoomd import *
import hoomd.md
from hoomd.md import *
import numpy as np
import hoomd.RPY
from hoomd import deprecated
import itertools
import os
from math import pi
from numpy import cos, sin

hoomd.context.initialize()

# assume hexagonal sheet (NOT a ribbon)
L = 39

thetas = [5] # angle about vorticity axis
phis = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90] # angle about flow axis
kbond_vals = np.array([1000]) * ((2*L - 1)/58) # Bond strength ~ shear
shear_vals = [1.00] # shear strength
kTeps_vals = [0.00]  # thermal energy/bending rigidity
ljf_vals = np.array([2])
protocols = [[[3, 0.03], [0.3, 0.03], [3, 0.03]], # Protocol 1 with K = 0.03
		  [[0.03, 0.03], [0.3, 0.03], [3, 0.03]]] # Protocol 2 with K = 0.03
		  # each row is a list of [lje, K] pairs. The simulation will shear at the first pair for 2000 shear cycles, 
		  # move to the next pair over 1000 cycles, 
		  # stay at that pair for 1000 cycles, and so on
runs = [1]
phi_old = phis[0]

# ONLY WORKS FOR kT = 0 RIGHT NOW
for param in itertools.product(kbond_vals, shear_vals, thetas, phis, kTeps_vals, ljf_vals, protocols, runs):
	# initialize hoomd context
	hoomd.context.initialize()

	kbond = param[0]
	shear = param[1]
	theta = param[2]
	phi = param[3]
	kTeps = param[4]
	ljf = param[5]
	protocol = param[6]
	runi = param[7]

	# update protocol number for file naming purposes
	if phi != phi_old:
		protocol_number = -1
	phi_old = phi

	protocol_number += 1

	radius = 1.0
	sigma = 2*np.sqrt(6)/3*radius*ljf
	kT = 0.00

	dt = 2e-4 # time step
	tf = 2000 # strain cycles at each [lje, K] pair

	dump_steps = int(1/dt/4)# dump every 0.5 diffusion times
	log_steps = int(1/dt/10)

	rest_steps = int(np.ceil(tf / dt)) # number of steps "resting" at each [lje, K] pair

	seed = hash((kbond, shear, theta, phi, kTeps, ljf, protocol_number, runi)) % (2**32 - 1)
	print('Seed is {}'.format(seed))

	# create folder for results
	foldername = '{}x{}_theta{}_phi{}_run{}'.format(L, L, theta, phi, runi)
	if ( not os.path.isdir(foldername) ):
		os.mkdir(foldername)
	filename = '{}/kbond{:.1f}_shear{:.2f}_ljf{:.4f}_protocol{}'.format(foldername, kbond, shear, ljf, protocol_number)
	print('Filename: {}'.format(filename))

	# set up system
	snap = hoomd.data.gsd_snapshot('hex_20x20.gsd')

	# rotate the sheet
	Rz = np.array([[cos(theta*pi/180), sin(theta*pi/180), 0], [-sin(theta*pi/180), cos(theta*pi/180), 0], [0, 0, 1]])
	Rx = np.array([[1, 0, 0], [0, cos(-phi*pi/180), sin(-phi*pi/180)], [0, -sin(-phi*pi/180), cos(-phi*pi/180)]])

	pos = snap.particles.position
	snap.particles.position[:] = np.dot(pos, Rx.T @ Rz.T)
	system = hoomd.init.read_snapshot(snap)

	all = group.all()
	integrate.mode_standard(dt=dt)

	# harmonic bond potential
	harmonic = bond.harmonic(name='bnd')
	harmonic.bond_coeff.set('bnd', k=kbond, r0=2*radius)

	# set up logging
	dump = hoomd.dump.gsd(filename='{}.gsd'.format(filename), period = dump_steps, group=all, overwrite=True)
	log = hoomd.analyze.log(filename='{}.dat'.format(filename), quantities=['bond_harmonic_energy', 'dihedral_harmonic_energy', 'pair_lj_energy', 'pressure_xx', 'pressure_xy', 'pressure_xz', 'pressure_yy', 'pressure_yz', 'pressure_zz'], period = log_steps, overwrite=True)

	# hard sphere potential
	def hs_potential_stokes(r, rmin, rmax, coeff1, coeff2, coeff3):
		V = 0.0
		F = -coeff1 * 8 * coeff2*coeff3 * r**3 * ( r - coeff2 - coeff3 )  / (3*r**4 + 6*(coeff2-coeff3)**2 * r**2 - (coeff2-coeff3)**4)
		return (V, F)

	nl = nlist.cell()
	table_hs = pair.table(width=1000, nlist=nl)
	table_hs.pair_coeff.set('A', 'A', func=hs_potential_stokes, rmin=0.001, rmax=2*radius, coeff=dict(coeff1=1/dt, coeff2=radius, coeff3=radius))

	# dihedral and LJ potential initialization
	dhdl = dihedral.harmonic()
	lj = pair.lj(r_cut = 2.5*sigma, nlist = nl) # Standard LJ cutoff

	# integrator setup
	bd = hoomd.RPY.integrate.RPY(group=all, T = kT, seed=seed)
	steady_shear = hoomd.RPY.shear_function.steady(dt = dt, shear_rate = shear)
	bd.set_params(function_form = steady_shear)
	box_resize = hoomd.update.box_resize(xy = hoomd.RPY.variant.shear_variant(function_form = steady_shear, total_timestep = rest_steps*len(protocol) + 1), scale_particles=False)

	# integration at each [lje, K] pair for given protocol
	for i in range(len(protocol)):
		[lje, keps] = protocol[i]

		lje = lje
		kappa = keps * lje * sigma**2 / radius**2  / 2 / np.sqrt(3) # keps gives k/eps tilde sigma^2

		# update dihedral potential
		dhdl.dihedral_coeff.set('dhdl', k=2*kappa, d=1, n=1) # multiply by 2 since it's defined with 1/2 in hoomd!
	
		# update LJ potential
		lj.pair_coeff.set('A', 'A', epsilon = lje, sigma = sigma, r_on = 2*radius)

		# integrate
		run(rest_steps)

	run(1)
