using PyCall
using KMesh

function generate_hexagonal_sheet(N1::Int64, N2::Int64, b::Float64 = 2.0)
	#= generate a hexagonal sheet lying in the x-z plane with with the number of
	edge beads equal to N1 and N2 along the x and x-z edges, respectively.
	The sheet is centered about the origin.

	INPUT
	N1 - number of beads along the x edge
	N2 - number of beads along the x-z edge (if equal to N1, this would be a perfect hexagon)
	b - length of edges

	OUTPUT
	p - 3 x Nbeads matrix containing the coordinates of vertices
	mi - MeshInfo containing triangulation info
	=#

	if N1 < 1 || N2 < 1
		throw(DomainError(min(N1,N2), "N1 and N2 must be greater than 0"))
	end

	# should rotate the hexagon if N2 is greater than N1
	should_rotate = false
	if N2 > N1
		swp = N1
		N1 = N2
		N2 = swp
		should_rotate = true
	end

	# create matrices for vertices and triangle connectivity
	Ncenter = N1 + N2 - 1
	Nbeads = (N1 + Ncenter)*N2 - Ncenter
	p = zeros(3, Nbeads)

	Ntri = (2*N1-1 + 2*N1-1 + 2*(N2 - 2))*(N2 - 1)
	t = zeros(Int64, 3, Ntri)

	# create center row
	leftx = -b*(Ncenter-1)/2
	p[1,1:Ncenter] .= b.*(0:(Ncenter-1)) .+ leftx

	# create other rows and add triangles
	cur = Ncenter + 1
	curt = 1
	last_bottom_ind = 1	# the beginning index of the last top and bottom rows added
	last_top_ind = 1

	for j = 1:(N2 - 1)
		Nbeads = Ncenter - j	# number of beads in this row
		row_leftx = leftx + b*j/2

		cur_bottom_ind = cur
		cur_top_ind = cur + Nbeads

		# calculate positions for bottom and top rows
		p[1, cur_bottom_ind:(cur_bottom_ind + Nbeads-1)] .= b.*(0:(Nbeads-1)) .+ row_leftx
		p[3, cur_bottom_ind:(cur_bottom_ind + Nbeads-1)] .= -b*j*sqrt(3)/2

		p[1, cur_top_ind:(cur_top_ind + Nbeads - 1)] .= b.*(0:(Nbeads-1)) .+ row_leftx
		p[3, cur_top_ind:(cur_top_ind + Nbeads - 1)] .= b*j*sqrt(3)/2

		# add triangles for bottom row
		for k = 0:Nbeads-2
			t[:,curt] .= [last_bottom_ind+k, last_bottom_ind+k+1, cur_bottom_ind+k]
			t[:,curt+1] .= [last_bottom_ind+k+1, cur_bottom_ind+k, cur_bottom_ind+k+1]
			curt += 2
		end

		# last triangle for bottom row
		t[:,curt] .= [last_bottom_ind+Nbeads-1, last_bottom_ind+Nbeads, cur_bottom_ind+Nbeads-1]
		curt += 1

		# add triangles for top row
		for k = 0:Nbeads-2
			t[:,curt] .= [last_top_ind+k, last_top_ind+k+1, cur_top_ind+k]
			t[:,curt+1] .= [last_top_ind+k+1, cur_top_ind+k, cur_top_ind+k+1]
			curt += 2
		end

		# last triangle for top row
		t[:,curt] .= [last_top_ind+Nbeads-1, last_top_ind+Nbeads, cur_top_ind+Nbeads-1]
		curt += 1

		# update counters
		cur = cur + 2*Nbeads
		last_bottom_ind = cur_bottom_ind
		last_top_ind = cur_top_ind
	end


	# rotate hexagon in x-z plane by 90 degrees
	if should_rotate
		@inbounds for j = 1:size(p, 2)
			swp = p[1,j];
			p[1, j] = p[3, j];
			p[3, j] = swap;
		end

	end

	mi = MeshInfo(t)
	order_triangles!(mi.t, mi.tnbor)
	return p, mi
end

function save_hoomd_gsd(p::Array{Float64,2}, mi::MeshInfo, filename::String, a::Float64 = 1.0)
	#= saves the gsd mesh into a hoomd gsd file with box size greater than twice
	the greatest diagonal in the object. Assumes the mesh is centered at the origin

	INPUT
	p - 3 x N array of bead positions
	mi - MeshInfo object containing bond connectivity
	filename - filename for gsd file
	[a] - bead radius (default is 1.0)
	=#

	# find smallest radius that contains mesh
	maxR = reduce(max, sqrt(p[i]^2 + p[i+1]^2 + p[i+2]^2) for i in 1:3:length(p))

	L = 4*maxR + 2*a

	# ------------- create hoomd snapshot ---------------------
	gsdhoomd = pyimport("gsd.hoomd")

	f = gsdhoomd.open(name=filename, mode="wb")

	snap = gsdhoomd.Snapshot()
	snap.configuration.box = [L, L, L, 0.0, 0.0, 0.0]

	# particles
	Np = size(p, 2)
	snap.particles.N = Np
	snap.particles.types = ["A"]
	snap.particles.typeid = repeat([0], Np)
	snap.particles.position = transpose(p)
	snap.particles.diameter = repeat([2*a], Np)

	# bonds
	Nbonds = sum(mi.nbor_count) ÷ 2
	bonds = zeros(Int64, Nbonds, 2)
	cur = 1
	for j = 1:size(mi.nbor, 2)
		for i = 1:mi.nbor_count[j]
			if mi.nbor[i,j] > j
				bonds[cur, :] = [j, mi.nbor[i,j]]
				cur += 1
			end
		end
	end

	snap.bonds.N = Nbonds
	snap.bonds.types = ["bnd"]
	snap.bonds.typeid = repeat([0], Nbonds)
	snap.bonds.group = bonds .- 1 	# for 0 indexing!

	# dihedrals
	Ndhdls = sum(i > 0 ? 1 : 0 for i in vec(mi.tnbor)) ÷ 2
	dhdls = zeros(Int64, Ndhdls, 4)
	cur = 1
	for t1 = 1:size(mi.tnbor, 2)
		for t2 in mi.tnbor[:, t1]
			if t2 > t1 || t2 == 0
				continue
			end

			i1, i2, i3 = mi.t[:, t1]
			j1, j2, j3 = mi.t[:, t2]

			# ensure i1 contains index of point not in triangle 2
			if i1 != j1 && i1 != j2 && i1 != j3 && i1 > 0
			elseif i2 != j1 && i2 != j2 && i2 != j3 && i2 > 0
				swp = i1;
				i1 = i2;
				i2 = swp;
			elseif i3 != j1 && i3 != j2 && i3 != j3 && i3 > 0
				swp = i1;
				i1 = i3;
				i3 = swp;
			end

			# ensure j1 contains index of point not in triangle 1
			if j1 != i1 && j1 != i2 && j1 != i3 && j1 > 0
			elseif j2 != i1 && j2 != i2 && j2 != i3 && j2 > 0
				swp = j1;
				j1 = j2;
				j2 = swp;
			elseif j3 != i1 && j3 != i2 && j3 != i3 && j3 > 0
				swp = j1;
				j1 = j3;
				j3 = swp;
			end

			dhdls[cur, :] = [i1, i2, i3, j1]
			cur += 1
		end
	end

	snap.dihedrals.N = Ndhdls
	snap.dihedrals.types = ["dhdl"]
	snap.dihedrals.typeid = repeat([0], Ndhdls)
	snap.dihedrals.group = dhdls .- 1		# for 0 indexing!

	snap.validate()
	f.append(snap)
end
