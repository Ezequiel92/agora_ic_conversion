### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 568ed730-9319-11ec-09ca-b1e962cadf22
using DelimitedFiles, Unitful, UnitfulAstro, HDF5, DataFrames

# ╔═╡ 5ca972c9-548a-445e-ba2e-e38151f7bedf
#####################################################################################
# Base units
#####################################################################################

begin
	const ILLUSTRIS_L_UNIT = 3.085678e21u"cm"
	const ILLUSTRIS_M_UNIT = 1.989e43u"g"
	const ILLUSTRIS_V_UNIT = 1.0e5u"cm*s^-1"
end;

# ╔═╡ d3602da9-8263-4638-978e-bb7da9998885
#####################################################################################
# Dimensions of specific energy
#####################################################################################

@derived_dimension SpecificEnergy Unitful.𝐋^2 * Unitful.𝐓^-2 true;

# ╔═╡ b20de40e-21fe-4d4e-b0dc-fcb0d4ead59e
#####################################################################################
# As an example we use the low resolution ICs from the AGORA project site
# https://sites.google.com/site/santacruzcomparisonproject/data
#####################################################################################

begin
    const out_file = "./output/ic_low"
    const SIM_COSMO = false
end;

# ╔═╡ dd16a119-4f6e-4510-ac4c-17e4ab2fdf80
"""
Unit conversion struct.

# Fields

  - `x_cgs::Unitful.Length`: Length, from internal units to ``\\mathrm{cm}``.
  - `x_cosmo::Unitful.Length`: Length, from internal units to ``\\mathrm{kpc}``.
  - `x_comoving::Unitful.Length`: Length, from internal units to ``\\mathrm{ckpc}``.
  - `v_cgs::Unitful.Velocity`: Velocity, from internal units to ``\\mathrm{cm \\, s^{-1}}``.
  - `v_cosmo::Unitful.Velocity`: Velocity, from internal units to ``\\mathrm{km \\, s^{-1}}``.
  - `m_cgs::Unitful.Mass`: Mass, from internal units to ``\\mathrm{g}``.
  - `m_cosmo::Unitful.Mass`: Mass, from internal units to ``\\mathrm{M_\\odot}``.
  - `t_cgs::Unitful.Time`: Time, from internal units to ``\\mathrm{s}``.
  - `t_cosmo::Unitful.Time`: Time, from internal units to ``\\mathrm{Myr}``.
  - `U_cgs::Unitful.Energy`: Specific energy, from internal units to ``\\mathrm{erg \\, g^{-1}}``.
  - `rho_cgs::Unitful.Density`: Density, from internal units to ``\\mathrm{g \\, cm^{-3}}``.
  - `P_Pa::Unitful.Pressure`: Pressure, from internal units to ``\\mathrm{Pa}``.
"""
struct InternalUnits

    x_cgs::Unitful.Length      # Length, from internal units to cm
    x_cosmo::Unitful.Length    # Length, from internal units to kpc
    x_comoving::Unitful.Length # Length, from internal units to ckpc

    v_cgs::Unitful.Velocity    # Velocity, from internal units to cm * s^-1
    v_cosmo::Unitful.Velocity  # Velocity, from internal units to km * s^-1

    m_cgs::Unitful.Mass        # Mass, from internal units to g
    m_cosmo::Unitful.Mass      # Mass, from internal units to M⊙

    t_cgs::Unitful.Time        # Time, from internal units to s
    t_cosmo::Unitful.Time      # Time, from internal units to Myr

    U_cgs::SpecificEnergy      # Specific energy, from internal units to erg * g^-1

    rho_cgs::Unitful.Density   # Density, from internal units to g * cm^-3

    P_Pa::Unitful.Pressure     # Pressure, from internal units to Pa

    """
        InternalUnits(; <keyword arguments>)

    Constructor for `InternalUnits`.

    # Arguments

      - `l_unit::Unitful.Length=ILLUSTRIS_L_UNIT`: Code parameter `UnitLength_in_cm`.
      - `m_unit::Unitful.Mass=ILLUSTRIS_M_UNIT`: Code parameter `UnitMass_in_g`.
      - `v_unit::Unitful.Velocity=ILLUSTRIS_V_UNIT`: Code parameter `UnitVelocity_in_cm_per_s`.
      - `a0::Float64=1.0`: Cosmological scale factor of the simulation.
      - `h0::Float64=1.0`: Hubble constant as "little h".
    """
    function InternalUnits(;
        l_unit::Unitful.Length=ILLUSTRIS_L_UNIT,
        m_unit::Unitful.Mass=ILLUSTRIS_M_UNIT,
        v_unit::Unitful.Velocity=ILLUSTRIS_V_UNIT,
        a0::Float64=1.0,
        h0::Float64=1.0,
    )

        #############################################################################
        # Base units
        #############################################################################

        x_cgs = l_unit * a0 / h0
        x_cosmo = x_cgs |> u"kpc"
        x_comoving = l_unit / h0 |> u"kpc"

        v_cgs = v_unit * sqrt(a0)
        v_cosmo = v_cgs |> u"km*s^-1"

        m_cgs = m_unit / h0
        m_cosmo = m_cgs |> u"Msun"

        #############################################################################
        # Derived units
        #############################################################################

        # Only used in non-cosmological simulations
        t_cgs = x_cgs / v_cgs
        t_cosmo = t_cgs |> u"Myr"

        U_cgs = v_unit^2 |> u"erg*g^-1"

        rho_cgs = m_cgs * x_cgs^-3

        # Thermal pressure (it uses v_unit^2 instead of v_cgs^2, 
		# which would add an extra factor of a0)
        P_Pa = v_unit^2 * m_cgs * x_cgs^-3 |> u"Pa"

        new(
            x_cgs,
            x_cosmo,
            x_comoving,
            v_cgs,
            v_cosmo,
            m_cgs,
            m_cosmo,
            t_cgs,
            t_cosmo,
            U_cgs,
            rho_cgs,
            P_Pa,
        )

    end

end;

# ╔═╡ a8db93bd-b83d-4236-ab15-a51604199ed6
"""
Data in the "Header" group of a HDF5 snapshot file.

# Fields

  - `npart::Vector{Int32}`: Number of particles (of each type) included in this file chunk.
  - `massarr::Vector{Float64}`: Masses of particle types which have a constant mass.
  - `time::Float64`: The physical time/scale factor.
  - `z::Float64`: Redshift of the simulation.
  - `flag_sfr::Int32`: 1 if the simulation was run with star formation, else 0.
  - `flag_feedback::Int32`: 1 if the simulation was run with stellar feedback, else 0.
  - `nall::Vector{UInt32}`: Total number of particles (of each type) for this snapshot.
  - `flag_cooling::Int32`: 1 if the simulation was run with cooling, else 0.
  - `num_files::Int32`: Number of file chunks per snapshot.	
  - `omega_0::Float64`: The cosmological density parameter for matter.	
  - `boxsize::Float64`: Total size of the simulation box.
  - `omega_l::Float64`: The cosmological density parameter for the cosmological constant.
  - `h0::Float64`: Hubble parameter.
  - `flag_stellarage::Int32`: 1 if the simulation was run with stellar age, else 0.
  - `flag_metals::Int32`: 1 if the simulation was run with metals, else 0.
  - `npartTotalHighWord::Vector{UInt32}`: If Npart > 1584^3 (> 2^32) this contains a high bit: ntotal = npartTotalHighWord * 2^32 + nall.
  - `flag_entropy_instead_u::Int32`: 1 if the snapshot U field contains entropy instead of internal energy, else 0.
  - `flag_doubleprecision::Int32`: 1 if the snapshot is in double precision, else 0.
  - `flag_ic_info::Int32`: 1 if the initial snapshot file contains an info block, else 0.
  - `lpt_scalingfactor::Float32`: Factor to use second order IC generation.
  - `fill::Vector{Int32}`: The HEAD block needs to be filled with zeros to have a size of 256 bytes.
"""
@kwdef mutable struct SnapshotHeader
	npart::Vector{Int32}
	massarr::Vector{Float64}
	time::Float64	
	z::Float64	
	flag_sfr::Int32	
	flag_feedback::Int32	
	nall::Vector{UInt32}
	flag_cooling::Int32	
	num_files::Int32	
	omega_0::Float64	
	boxsize::Float64	
	omega_l::Float64	
	h0::Float64	
	flag_stellarage::Int32
	flag_metals::Int32	
	npartTotalHighWord::Vector{UInt32}	
	flag_entropy_instead_u::Int32	
	flag_doubleprecision::Int32	
	flag_ic_info::Int32	
	lpt_scalingfactor::Float32	
	fill::Vector{Int32}	
end;

# ╔═╡ 6647922f-bcc9-4438-a454-d8dba7a2d103
#####################################################################################
# Read IC files which are in the following format:
#
# Velocity: km/s
# Mass:     10^9 Msun
# Length:   kpc
#
# Gas particle     (gas.dat):   x, y, z, vx, vy, vz, mgas, u_gas
# Dark matter halo (halo.dat):  x, y, z, vx, vy, vz, mdark
# Stellar disk     (disk.dat):  x, y, z, vx, vy, vz, mdisk
# Stellar bulge    (bulge.dat): x, y, z, vx, vy, vz, mbulge
#####################################################################################

begin
    rawIC_type0 = readdlm("./AGORA_ICs/LOW/gas.dat")
    rawIC_type1 = readdlm("./AGORA_ICs/LOW/halo.dat")
    rawIC_type2 = readdlm("./AGORA_ICs/LOW/disk.dat")
    rawIC_type3 = readdlm("./AGORA_ICs/LOW/bulge.dat")

    s0 = size(rawIC_type0, 1)
    s1 = size(rawIC_type1, 1)
    s2 = size(rawIC_type2, 1)
    s3 = size(rawIC_type3, 1)
end;

# ╔═╡ f3cc08bd-b2c6-42a2-a757-5ef4e7c685e9
#####################################################################################
# Header
#####################################################################################

header = SnapshotHeader(
    npart                  = Int32[s0, s1, s2, s3, 0, 0],
    massarr                = [
        rawIC_type0[1, 7],
        rawIC_type1[1, 7],
        rawIC_type2[1, 7],
        rawIC_type3[1, 7],
        0.0,
        0.0,
    ] .* 10^9*u"Msun" ./ (ILLUSTRIS_M_UNIT |> u"Msun"),
    time                   = 0.0,
    z                      = 0.0,
	flag_sfr               = convert(Int32, 1),
	flag_feedback          = convert(Int32, 1),
	nall                   = UInt32[s0, s1, s2, s3, 0, 0],
	flag_cooling           = convert(Int32, 1),
	num_files              = convert(Int32, 1),
	omega_0                = 0.0,
	boxsize                = 0.0,
	omega_l                = 0.0,
	h0                     = 1.0,
	flag_stellarage        = convert(Int32, 1),
	flag_metals            = convert(Int32, 1),
	npartTotalHighWord     = UInt32[0, 0, 0, 0, 0, 0],
	flag_entropy_instead_u = convert(Int32, 0),
	flag_doubleprecision   = convert(Int32, 0),
	flag_ic_info           = convert(Int32, 0),
	lpt_scalingfactor      = 0.0f0,
	fill                   = Int32[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
);

# ╔═╡ a0f96fb1-8fda-410f-8861-563803f3798e
#####################################################################################
# Unit struct
#####################################################################################

IU = InternalUnits(; a0=SIM_COSMO ? header.time : 1.0, h0=header.h0);

# ╔═╡ 108bb88d-0d89-4601-a2b2-3af50fcd1f3b
#####################################################################################
# Write the ICs in HDF5 format for GADGET (SnapFormat = 3)
#
# Within each block, the particles will be ordered according to their particle type,
# i.e. gas particles will come first (type 0), then DM (type 1) particles, followed
# by disk (type 2) particles, and so on:
#
# gas   => 0
# halo  => 1
# disk  => 2
# bulge => 3
#####################################################################################

begin
	# Positions
    pos = convert(
		Array{Float32,2},
		vcat(
            rawIC_type0[:, 1:3],
			rawIC_type1[:, 1:3],
			rawIC_type2[:, 1:3],
			rawIC_type3[:, 1:3],
		)' .* u"kpc" ./ IU.x_cosmo,
	)

	# Velocities
	vel = convert(
		Array{Float32,2},
		vcat(
			rawIC_type0[:, 4:6],
			rawIC_type1[:, 4:6],
			rawIC_type2[:, 4:6],
			rawIC_type3[:, 4:6],
		)' .* u"km*s^-1" ./ IU.v_cosmo,
	)

	# Internal energy (for gas particles only)
	u_factor = uconvert(u"erg*g^-1", 1.0*u"km*s^-1"^2)
	u = convert(
		Array{Float32,1},
		rawIC_type0[:, 8] .* u_factor ./ IU.U_cgs,
	)

	# IDs
	id = convert(Array{UInt32,1}, 1:(s0 + s1 + s2 + s3))

	h5open(out_file * ".hdf5", "w") do fh5
		g0 = create_group(fh5, "PartType0")
		g1 = create_group(fh5, "PartType1")
		g2 = create_group(fh5, "PartType2")
		g3 = create_group(fh5, "PartType3")
		head = create_group(fh5, "Header")
	
		g0["Coordinates"] = collect(pos[:, 1:s0]')
		g1["Coordinates"] = collect(pos[:, (s0 + 1):(s0 + s1)]')
		g2["Coordinates"] = collect(pos[:, (s0 + s1 + 1):(s0 + s1 + s2)]')
		g3["Coordinates"] = collect(pos[:, (s0 + s1 + s2 + 1):(s0 + s1 + s2 + s3)]')
	
		g0["Velocities"] = collect(vel[:, 1:s0]')
		g1["Velocities"] = collect(vel[:, (s0 + 1):(s0 + s1)]')
		g2["Velocities"] = collect(vel[:, (s0 + s1 + 1):(s0 + s1 + s2)]')
		g3["Velocities"] = collect(vel[:, (s0 + s1 + s2  + 1):(s0 + s1 + s2 + s3)]')
	
		g0["ParticleIDs"] = id[1:s0]
		g1["ParticleIDs"] = id[(s0 + 1):(s0 + s1)]
		g2["ParticleIDs"] = id[(s0 + s1 + 1):(s0 + s1 + s2)]
		g3["ParticleIDs"] = id[(s0 + s1 + s2 + 1):(s0 + s1 + s2 + s3)]
	
		write_attribute(head, "NumPart_ThisFile", header.npart)
		write_attribute(head, "NumPart_Total", header.nall)
		write_attribute(head, "MassTable", header.massarr)
		write_attribute(head, "Time", header.time)
		write_attribute(head, "Redshift", header.z)
		write_attribute(head, "BoxSize", header.boxsize)
		write_attribute(head, "NumFilesPerSnapshot", header.num_files)
		write_attribute(head, "NumPart_Total_HighWord", header.npartTotalHighWord)
		write_attribute(head, "Omega0", header.omega_0)
	    write_attribute(head, "OmegaLambda", header.omega_l)
	    write_attribute(head, "HubbleParam", header.h0)
		write_attribute(head, "Flag_Sfr", header.flag_sfr)
		write_attribute(head, "Flag_Cooling", header.flag_cooling)
		write_attribute(head, "Flag_StellarAge", header.flag_stellarage)
		write_attribute(head, "Flag_Feedback", header.flag_feedback)
		write_attribute(head, "Flag_DoublePrecision", header.flag_doubleprecision)
		write_attribute(head, "Flag_Metals", header.flag_metals)
	end
end;

# ╔═╡ 39f2147d-ac35-4b14-ade6-ee07d17f84c6
#####################################################################################
# Consistency test
#####################################################################################

h5open(out_file * ".hdf5", "r") do fh5
	hdf5_pos = fh5["PartType0/Coordinates"][150:160, 2] * IU.x_cosmo
	hdf5_vel = fh5["PartType0/Velocities"][150:160, 2] * IU.v_cosmo

	@assert all(
		isapprox.(rawIC_type0[150:160, 2], ustrip(hdf5_pos), rtol = 10^-5)
	) && all(
		isapprox.(rawIC_type0[150:160, 5], ustrip(hdf5_vel), rtol = 10^-5)
	)
end;

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
HDF5 = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
UnitfulAstro = "6112ee07-acf9-5e0f-b108-d242c714bf9f"

[compat]
DataFrames = "~1.7.0"
DelimitedFiles = "~1.9.1"
HDF5 = "~0.17.2"
Unitful = "~1.22.0"
UnitfulAstro = "~1.2.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.4"
manifest_format = "2.0"
project_hash = "06464aa013b14e666607858ea8d75a76b4bd591f"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "fb61b4812c49343d7ef0b533ba982c46021938a6"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.7.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.HDF5]]
deps = ["Compat", "HDF5_jll", "Libdl", "MPIPreferences", "Mmap", "Preferences", "Printf", "Random", "Requires", "UUIDs"]
git-tree-sha1 = "e856eef26cf5bf2b0f95f8f4fc37553c72c8641c"
uuid = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
version = "0.17.2"

    [deps.HDF5.extensions]
    MPIExt = "MPI"

    [deps.HDF5.weakdeps]
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"

[[deps.HDF5_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "LibCURL_jll", "Libdl", "MPICH_jll", "MPIPreferences", "MPItrampoline_jll", "MicrosoftMPI_jll", "OpenMPI_jll", "OpenSSL_jll", "TOML", "Zlib_jll", "libaec_jll"]
git-tree-sha1 = "e94f84da9af7ce9c6be049e9067e511e17ff89ec"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.14.6+0"

[[deps.Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f93a9ce66cd89c9ba7a4695a47fd93b4c6bc59fa"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.12.0+0"

[[deps.InlineStrings]]
git-tree-sha1 = "6a9fde685a7ac1eb3495f8e812c5a7c3711c2d5e"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.3"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MPICH_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "3aa3210044138a1749dbd350a9ba8680869eb503"
uuid = "7cb0a576-ebde-5e09-9194-50597f1243b4"
version = "4.3.0+1"

[[deps.MPIPreferences]]
deps = ["Libdl", "Preferences"]
git-tree-sha1 = "c105fe467859e7f6e9a852cb15cb4301126fac07"
uuid = "3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"
version = "0.1.11"

[[deps.MPItrampoline_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "ff91ca13c7c472cef700f301c8d752bc2aaff1a8"
uuid = "f1f71cc9-e9ae-5b93-9b94-4fe0e1ad3748"
version = "5.5.3+0"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.MicrosoftMPI_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bc95bf4149bf535c09602e3acdf950d9b4376227"
uuid = "9237b28f-5490-5468-be7b-bb81f5f5e6cf"
version = "10.1.4+3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenMPI_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML", "Zlib_jll"]
git-tree-sha1 = "da913f03f17b449951e0461da960229d4a3d1a8c"
uuid = "fe0851c0-eecd-5654-98d4-656369965a5c"
version = "5.0.7+1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a9697f1d06cc3eb3fb3ad49cc67f2cfabaac31ea"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.16+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "cc4054e898b852042d7b503313f7ad03de99c3dd"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1101cd475833706e4d0e7b122218257178f48f34"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "712fb0231ee6f9120e005ccd56297abbc053e7e0"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.8"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "725421ae8e530ec29bcbdddbe91ff8053421d023"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "c0667a8e676c53d390a09dc6870b3d8d6650e2bf"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.22.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulAngles]]
deps = ["Dates", "Unitful"]
git-tree-sha1 = "79875b1f2e4bf918f0702a5980816955066d9ae2"
uuid = "6fb2a4bd-7999-5318-a3b2-8ad61056cd98"
version = "0.7.2"

[[deps.UnitfulAstro]]
deps = ["Unitful", "UnitfulAngles"]
git-tree-sha1 = "da7577e6a726959b14f7451674d00b78d10ca30f"
uuid = "6112ee07-acf9-5e0f-b108-d242c714bf9f"
version = "1.2.1"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libaec_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f5733a5a9047722470b95a81e1b172383971105c"
uuid = "477f73a3-ac25-53e9-8cc3-50b2fa2566f0"
version = "1.1.3+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╠═568ed730-9319-11ec-09ca-b1e962cadf22
# ╠═5ca972c9-548a-445e-ba2e-e38151f7bedf
# ╠═d3602da9-8263-4638-978e-bb7da9998885
# ╠═b20de40e-21fe-4d4e-b0dc-fcb0d4ead59e
# ╠═dd16a119-4f6e-4510-ac4c-17e4ab2fdf80
# ╠═a8db93bd-b83d-4236-ab15-a51604199ed6
# ╠═6647922f-bcc9-4438-a454-d8dba7a2d103
# ╠═a0f96fb1-8fda-410f-8861-563803f3798e
# ╠═f3cc08bd-b2c6-42a2-a757-5ef4e7c685e9
# ╠═108bb88d-0d89-4601-a2b2-3af50fcd1f3b
# ╠═39f2147d-ac35-4b14-ade6-ee07d17f84c6
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
