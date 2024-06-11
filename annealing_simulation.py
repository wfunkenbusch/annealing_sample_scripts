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

L = 39

thetas = [5] # angle about vorticity axis
phis = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90] # angle about flow axis
kbond_vals = np.array([1000]) * ((2*L - 1)/58) # Bond strength ~ shear
keps_vals = np.array([0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0]) # K
shear_vals = [1.00] # shear strength
kTeps_vals = [0.00] # thermal energy/bending rigidity
lje_vals = np.array([0.1]) # LJ strength between beads
ljf_vals = np.array([2]) # LJ range as a multiple of 2 sqrt(6)/3
tq_vals = np.array([10.0, 100.0, 1000.0, 10000.0]) # quench times
runs = [1]

for param in itertools.product(kbond_vals, shear_vals, keps_vals, thetas, phis, kTeps_vals, lje_vals, ljf_vals, tq_vals, runs):
	# initialize hoomd context
	hoomd.context.initialize()

	kbond = param[0]
	shear = param[1]
	keps = param[2]
	theta = param[3]
	phi = param[4]
	kTeps = param[5]
	lje = param[6]
	ljf = param[7]
	tq = param[8]
	runi = param[9]

	radius = 1.0 # bead radius
	sigma = 2*np.sqrt(6)/3*radius*ljf # LJ range
	lje = lje
	kT = kTeps * lje * sigma**2 / radius**2 / 2 / np.sqrt(3) # kTeps gives kT/eps tilde sigma^2
	kappa = keps * lje * sigma**2 / radius**2  / 2 / np.sqrt(3) # keps gives k/eps tilde sigma^2

	dt = 2e-4 # time step
	tf = tq
	dshear = 0.001 # decrease in 1000 steps
	nsteps = int(dshear / (shear * dt / tf)) # number of steps at each shear rate
	rest_steps = int(np.ceil(2000/dt)) # time to initialize and rest

	seed = hash((kbond, kappa, shear, theta, phi, kT, runi, lje, tq, ljf)) % (2**32 - 1)
	print('Seed is {}'.format(seed))

	dump_steps = int(1/dt/4) # dump every 0.5 diffusion times
	log_steps = int(1/dt/10)

	# create folder for results
	foldername = '{}x{}_theta{}_phi{}_run{}'.format(L, L, theta, phi, runi)
	if ( not os.path.isdir(foldername) ):
		os.mkdir(foldername)
	filename = '{}/tq{:.2f}_kbond{:.1f}_keps{:.4f}_shear{:.2f}_kT{:.4f}_lje{:.4f}_ljf{:.4f}'.format(foldername, tq, kbond, keps, shear, kT, lje, ljf)
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

	# dihedral potential
	dhdl = dihedral.harmonic()
	dhdl.dihedral_coeff.set('dhdl', k=2*kappa, d=1, n=1) # multiply by 2 since it's defined with 1/2 in hoomd!

	# hard sphere potential
	def hs_potential_stokes(r, rmin, rmax, coeff1, coeff2, coeff3):
		V = 0.0
		F = -coeff1 * 8 * coeff2*coeff3 * r**3 * ( r - coeff2 - coeff3 )  / (3*r**4 + 6*(coeff2-coeff3)**2 * r**2 - (coeff2-coeff3)**4)
		return (V, F)

	nl = nlist.cell()
	table_hs = pair.table(width=1000, nlist=nl)
	table_hs.pair_coeff.set('A', 'A', func=hs_potential_stokes, rmin=0.001, rmax=2*radius, coeff=dict(coeff1=1/dt, coeff2=radius, coeff3=radius))
	
	# LJ potential
	lj = pair.lj(r_cut = 2.5*sigma, nlist = nl) # Standard LJ cutoff
	lj.pair_coeff.set('A', 'A', epsilon = lje, sigma = sigma, r_on = 2*radius)

	# set up logging
	dump = hoomd.dump.gsd(filename='{}.gsd'.format(filename), period = dump_steps, group=all, overwrite=True)
	log = hoomd.analyze.log(filename='{}.dat'.format(filename), quantities=['bond_harmonic_energy', 'dihedral_harmonic_energy', 'pair_lj_energy', 'pressure_xx', 'pressure_xy', 'pressure_xz', 'pressure_yy', 'pressure_yz', 'pressure_zz'], period = log_steps, overwrite=True)
	
	# integrator setup
	bd = hoomd.RPY.integrate.RPY(group = all, T = kT, seed = seed)
	steady_shear = hoomd.RPY.shear_function.steady(dt = dt, shear_rate = shear)
	bd.set_params(function_form = steady_shear)
	box_resize = hoomd.update.box_resize(xy = hoomd.RPY.variant.shear_variant(function_form = steady_shear, total_timestep = rest_steps), scale_particles = False)
	run(rest_steps) # initialization
	
	box_resize.disable()

	# integrate while shearing the box with the specified shear rate
	for i in range(int(1/dshear)):
		steady_shear = hoomd.RPY.shear_function.steady(dt = dt, shear_rate = shear)
		bd.set_params(function_form = steady_shear)
		box_resize = hoomd.update.box_resize(xy = hoomd.RPY.variant.shear_variant(function_form = steady_shear, total_timestep = nsteps), scale_particles = False)

		run(nsteps)

		box_resize.disable()
		shear += -dshear # update shear rate

	# no shear simulation
	steady_shear = hoomd.RPY.shear_function.steady(dt = dt, shear_rate = shear)
	bd.set_params(function_form = steady_shear)
	box_resize = hoomd.update.box_resize(xy = hoomd.RPY.variant.shear_variant(function_form = steady_shear, total_timestep = 1), scale_particles = False)

	run(rest_steps + 1) # resting

