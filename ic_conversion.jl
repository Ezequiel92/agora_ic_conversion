### A Pluto.jl notebook ###
# v1.0.3

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
DataFrames = "~1.8.2"
HDF5 = "~0.17.4"
Unitful = "~1.29.0"
UnitfulAstro = "~1.2.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.7"
manifest_format = "2.0"
project_hash = "29a1fde9069679e533664342824b53c189727961"

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
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.1+2"

[[deps.Crayons]]
git-tree-sha1 = "54b76cbb40d9a0f5368c880725b2f141da77c94f"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.2.0"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "5fab31e2e01e70ad66e3e24c968c264d1cf166d6"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.8.2"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "b0bc6d2cad1fed8b7fd59a1551a991cb3d2809e6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.6"

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
version = "1.7.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.HDF5]]
deps = ["Compat", "HDF5_jll", "Libdl", "MPIPreferences", "Mmap", "Preferences", "Printf", "Random", "Requires", "UUIDs"]
git-tree-sha1 = "26e37af34e8ccb7a8358bd8de1619c7bcb526738"
uuid = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
version = "0.17.4"

    [deps.HDF5.extensions]
    MPIExt = "MPI"

    [deps.HDF5.weakdeps]
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"

[[deps.HDF5_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LibCURL_jll", "Libdl", "MPIABI_jll", "MPICH_jll", "MPIPreferences", "MPItrampoline_jll", "MicrosoftMPI_jll", "OpenMPI_jll", "OpenSSL_jll", "TOML", "Zlib_jll", "aws_c_s3_jll", "dlfcn_win32_jll", "libaec_jll", "mpif_jll"]
git-tree-sha1 = "194d676302b9b6aa53ea1f98ae8607d5caa8de4f"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "2.2.2+0"

[[deps.Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XML2_jll", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "8015ef94425626913eaa6933eb32e7c9bbe08030"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.14.1+0"

[[deps.InlineStrings]]
git-tree-sha1 = "06b65886c7577a3784d616e29f1302c2e36e389d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.6"

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
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f88f3ccef05a6a72a0cf0ed417c8fd68530f4ab2"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MPIABI_jll]]
deps = ["Artifacts", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "c0ab44a826d1a8219715078f79c265a11292db86"
uuid = "b5ada748-db0f-5fc0-8972-9331c762740c"
version = "1.0.0+0"

[[deps.MPICH_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "07dbec8aab01696edc0151a401a6cdfe95b9b885"
uuid = "7cb0a576-ebde-5e09-9194-50597f1243b4"
version = "5.0.1+0"

[[deps.MPIPreferences]]
deps = ["Libdl", "Preferences"]
git-tree-sha1 = "8e98d5d80b87403c311fd51e8455d4546ba7a5f8"
uuid = "3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"
version = "0.1.12"

[[deps.MPItrampoline_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "675df097f8eeb28998b2cfe3b25655af73d5f7df"
uuid = "f1f71cc9-e9ae-5b93-9b94-4fe0e1ad3748"
version = "5.5.6+0"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

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
version = "2025.11.4"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenMPI_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML", "Zlib_jll"]
git-tree-sha1 = "fb9be749680dd1283049ee17d96bd0ec611bc50f"
uuid = "fe0851c0-eecd-5654-98d4-656369965a5c"
version = "5.0.12+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05f45c2e0de6259db764adbfd2f1dc6d3f8de13c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "2.0.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "edbeefc7a4889f528644251bdb5fc9ab5348bc2c"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "5005266de4bfe50e53ff44a5cb5c540b6e47a254"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.6.0"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "REPL", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1b8aa19f229b1cea7fc93874a52e49db6a854450"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "3.4.8"

    [deps.PrettyTables.extensions]
    PrettyTablesExcelExt = "XLSX"
    PrettyTablesTypstryExt = "Typstry"

    [deps.PrettyTables.weakdeps]
    Typstry = "f0ed7684-a786-439e-b1e3-3b82803b501e"
    XLSX = "fdbf4ff8-1666-58a4-91e7-1b58723a45e0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
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
git-tree-sha1 = "084c47c7c5ce5cfecefa0a98dff69eb3646b5a80"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.10"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "13cd91cc9be159e3f4d95b857fa2aa383b53772a"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "e2b53ce13a53367e96601081e33d34746b571bad"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.5"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "773065c6e0e903924a9d838259be74338422aef2"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.5.0"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

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
git-tree-sha1 = "a94d9bdda1b7bed0046cea645639ab3f62196fac"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.14.0"

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
git-tree-sha1 = "1f0f9f401753701a7e4113b5056ca38d33875b55"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.29.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    LatexifyExt = ["Latexify", "LaTeXStrings"]
    NaNMathExt = "NaNMath"
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.UnitfulAngles]]
deps = ["Dates", "Unitful"]
git-tree-sha1 = "79875b1f2e4bf918f0702a5980816955066d9ae2"
uuid = "6fb2a4bd-7999-5318-a3b2-8ad61056cd98"
version = "0.7.2"

[[deps.UnitfulAstro]]
deps = ["Unitful", "UnitfulAngles"]
git-tree-sha1 = "fbe44a0ade62ae5ed0240ad314dfdd5482b90b40"
uuid = "6112ee07-acf9-5e0f-b108-d242c714bf9f"
version = "1.2.2"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "b826aedca5fce2e3c91904e338c279e5c1a2da74"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.15.4+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "58972370b81423fc546c56a60ed1a009450177c3"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.19.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.aws_c_auth_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "aws_c_cal_jll", "aws_c_http_jll", "aws_c_sdkutils_jll"]
git-tree-sha1 = "8cab83c96af80a1be968251ce1a0548a7545484d"
uuid = "2b3700d1-4306-52e2-a478-c162f0c514be"
version = "0.9.6+0"

[[deps.aws_c_cal_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "aws_c_common_jll"]
git-tree-sha1 = "22c0f42f4a1f0dc5dcfa8fd267c4ac407c455e7a"
uuid = "70f11efc-bab2-57f1-b0f3-22aad4e67c4b"
version = "0.9.13+0"

[[deps.aws_c_common_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a759cb9bf456ad792cc7898a81ae333cce9ef02a"
uuid = "73048d1d-b8c4-5092-a58d-866c5e8d1e50"
version = "0.12.6+0"

[[deps.aws_c_compression_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "aws_c_common_jll"]
git-tree-sha1 = "7910c72f45f44afd297c39fe43b99c56d5ed22ec"
uuid = "73a04cd5-f3d7-5bac-9290-e8adb709f224"
version = "0.3.2+0"

[[deps.aws_c_http_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "aws_c_compression_jll", "aws_c_io_jll"]
git-tree-sha1 = "3fb8685778068de502c72fec5dd8075e037cee15"
uuid = "3254fc65-9028-534d-aa9d-d76d128babc6"
version = "0.10.15+0"

[[deps.aws_c_io_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "aws_c_cal_jll", "aws_c_common_jll", "s2n_tls_jll"]
git-tree-sha1 = "7e481d474b2087ee8bbf55b81bf9119f21e396d9"
uuid = "13c41daa-f319-5298-b5eb-5754e0170d52"
version = "0.26.3+0"

[[deps.aws_c_s3_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "aws_c_auth_jll", "aws_c_common_jll", "aws_c_http_jll", "aws_checksums_jll", "s2n_tls_jll"]
git-tree-sha1 = "3e9917ab25114feba657e71be41cad068b9f6595"
uuid = "bd1f34fb-993f-5903-a121-aaf302eed6d4"
version = "0.11.5+0"

[[deps.aws_c_sdkutils_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "aws_c_common_jll"]
git-tree-sha1 = "c43dfba2c1ab9ea9f02f2c80e86fa16f6460244e"
uuid = "1282aa60-004d-510b-9f52-12498d409daa"
version = "0.2.4+1"

[[deps.aws_checksums_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "aws_c_common_jll"]
git-tree-sha1 = "2570c8e23f4771a087b12a47edcaaa670ac05a01"
uuid = "b2a88e68-78e7-5e94-8c20-c02986ec140e"
version = "0.2.10+0"

[[deps.dlfcn_win32_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e141d67ffe550eadfb5af1bdbdaf138031e4805f"
uuid = "c4b69c83-5512-53e3-94e6-de98773c479f"
version = "1.4.2+0"

[[deps.libaec_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "60f4792734488db6f42e2c7699f1d4594780bd03"
uuid = "477f73a3-ac25-53e9-8cc3-50b2fa2566f0"
version = "1.1.7+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.mpif_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIABI_jll", "MPICH_jll", "MPIPreferences", "MPItrampoline_jll", "MicrosoftMPI_jll", "OpenMPI_jll", "TOML"]
git-tree-sha1 = "a06fcd368cfe6fe2c0eb7b63320d4d27ddcd010d"
uuid = "9aeb927a-4695-514f-a259-621a69f20ec0"
version = "1.0.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"

[[deps.s2n_tls_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8ee2baf330dde17c6b711cee37dc18053433ddcf"
uuid = "cddc5d3d-934d-5d3a-9747-62fc12ea3f48"
version = "1.7.10+0"
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
