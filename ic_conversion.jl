### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ 568ed730-9319-11ec-09ca-b1e962cadf22
using DelimitedFiles, GadgetIO, GadgetUnits, Unitful, UnitfulAstro, HDF5, DataFrames

# ╔═╡ b20de40e-21fe-4d4e-b0dc-fcb0d4ead59e
#
# As an example we use the low resolution ICs from the agora project
# https://sites.google.com/site/santacruzcomparisonproject/data
#
begin
    const out_file = "./output/ic_low"
    const SIM_COSMO = false
end;

# ╔═╡ 6647922f-bcc9-4438-a454-d8dba7a2d103
#
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
#
begin
    rawIC_type0 = readdlm("./agora_data/LOW/gas.dat")
    rawIC_type1 = readdlm("./agora_data/LOW/halo.dat")
    rawIC_type2 = readdlm("./agora_data/LOW/disk.dat")
    rawIC_type3 = readdlm("./agora_data/LOW/bulge.dat")

    s0 = size(rawIC_type0, 1)
    s1 = size(rawIC_type1, 1)
    s2 = size(rawIC_type2, 1)
    s3 = size(rawIC_type3, 1)
end;

# ╔═╡ f3cc08bd-b2c6-42a2-a757-5ef4e7c685e9
#
# Header
# Fields: https://ludwigboess.github.io/GadgetIO.jl/stable/file_infos/
#
header = SnapshotHeader(
    Int32[s0, s1, s2, s3, 0, 0],
    [
        rawIC_type0[1, 7], 
        rawIC_type1[1, 7], 
        rawIC_type2[1, 7], 
        rawIC_type3[1, 7],
        0.0,
        0.0,
    ] .* 10^9*u"Msun" ./ GadgetPhysicalUnits(a_scale = 1.0, hpar = 1.0).m_msun,
    0.0,
    0.0,
	convert(Int32, 1),
	convert(Int32, 1),
	UInt32[s0, s1, s2, s3, 0, 0],
	convert(Int32, 1),
	convert(Int32, 1),
	0.0,
	0.0,
	0.0,
	1.0,
	convert(Int32, 1),
	convert(Int32, 1),
	UInt32[0, 0, 0, 0, 0, 0],
	convert(Int32, 0),
	convert(Int32, 0),
	convert(Int32, 0),
	0.0f0,
	Int32[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
);

# ╔═╡ a0f96fb1-8fda-410f-8861-563803f3798e
#
# Unit struct
#
# Fields: https://ludwigboess.github.io/GadgetUnits.jl/stable/conversion_structs/
#
GU = GadgetPhysicalUnits(a_scale = SIM_COSMO ? header.time : 1.0, hpar = header.h0);

# ╔═╡ 108bb88d-0d89-4601-a2b2-3af50fcd1f3b
# 
# Write the ICs in the GADGET binary format (SnapFormat = 2)
#
# Within each block, the particles will be ordered according to their particle type,  
# i.e. gas particles will come first (type 0), then DM (type 1) particles, followed 
# by disk (type 2) particles, and so on:
#
# gas   => 0
# halo  => 1
# disk  => 2
# bulge => 3
#
begin
	# Positions
    pos = convert(
		Array{Float32,2}, 
		vcat(
            rawIC_type0[:, 1:3],
			rawIC_type1[:, 1:3],
			rawIC_type2[:, 1:3],
			rawIC_type3[:, 1:3],
		)' .* u"kpc" ./ GU.x_physical,
	)

	# Velocities
	vel = convert(
		Array{Float32,2}, 
		vcat(
			rawIC_type0[:, 4:6],
			rawIC_type1[:, 4:6],
			rawIC_type2[:, 4:6],
			rawIC_type3[:, 4:6],
		)' .* u"km/s" ./ GU.v_physical,
	)

	# Internal energy (for gas particles only)
	u_factor = uconvert(u"erg/Msun", 1.0*u"km/s"^2)
	u = convert(
		Array{Float32,1}, 
		rawIC_type0[:, 8] .* u_factor ./ (GU.E_cgs / GU.m_msun),
	)

	# IDs
	id = convert(Array{UInt32,1}, 1:(s0+s1+s2+s3))

	# Write output IC file
	open(out_file, "w") do f
		write_header(f, header)
		write_block(f, pos, "POS")
		write_block(f, vel, "VEL")
		write_block(f, id, "ID")
		write_block(f, u, "U")
	end
end;

# ╔═╡ 6ea6a52d-1ee4-4069-a5ca-4b0021a7b69c
#
# Write the ICs in HDF5 format for GADGET (SnapFormat = 3)
#
h5open(out_file * ".hdf5", "w") do fh5
	g0 = create_group(fh5, "PartType0")
	g1 = create_group(fh5, "PartType1")
	g2 = create_group(fh5, "PartType2")
	g3 = create_group(fh5, "PartType3")
	head = create_group(fh5, "Header")
		
	g0["Coordinates"] = collect(pos[:, 1:s0]')
	g1["Coordinates"] = collect(pos[:, (s0+1):(s0+s1)]')
	g2["Coordinates"] = collect(pos[:, (s0+s1+1):(s0+s1+s2)]')
	g3["Coordinates"] = collect(pos[:, (s0+s1+s2+1):(s0+s1+s2+s3)]')
		
	g0["Velocities"] = collect(vel[:, 1:s0]')
	g1["Velocities"] = collect(vel[:, (s0+1):(s0+s1)]')
	g2["Velocities"] = collect(vel[:, (s0+s1+1):(s0+s1+s2)]')
	g3["Velocities"] = collect(vel[:, (s0+s1+s2+1):(s0+s1+s2+s3)]')

	g0["ParticleIDs"] = id[1:s0]
	g1["ParticleIDs"] = id[(s0+1):(s0+s1)]
	g2["ParticleIDs"] = id[(s0+s1+1):(s0+s1+s2)]
	g3["ParticleIDs"] = id[(s0+s1+s2+1):(s0+s1+s2+s3)]
		
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
end;

# ╔═╡ 39f2147d-ac35-4b14-ade6-ee07d17f84c6
#
# Consistency test
#
begin
	
	info_pos = InfoLine("POS", Float32, 3, Int32[1, 1, 1, 1, 0, 0])
	input_pos = read_block(out_file, "POS", info = info_pos, parttype = 0)
	
	info_vel = InfoLine("VEL", Float32, 3, Int32[1, 1, 1, 1, 0, 0])
	input_vel = read_block(out_file, "VEL", info = info_vel, parttype = 0)
	
	@assert all(isapprox.(
		rawIC_type0[150:160, 2], 
		ustrip(input_pos[2, 150:160] * GU.x_physical), 
		rtol = 10^-5,
	)) && all(isapprox.(
		rawIC_type0[150:160, 5], 
		ustrip(input_vel[2, 150:160] * GU.v_physical), 
		rtol = 10^-5,
	))

	h5open(out_file * ".hdf5", "r") do fh5
		hdf5_pos = fh5["PartType0/Coordinates"][150:160, 2] * GU.x_physical
		hdf5_vel = fh5["PartType0/Velocities"][150:160, 2] * GU.v_physical

		@assert all(
			isapprox.(rawIC_type0[150:160, 2], ustrip(hdf5_pos), rtol = 10^-5)
		) && all(
			isapprox.(rawIC_type0[150:160, 5], ustrip(hdf5_vel), rtol = 10^-5)
		)
	end

end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
GadgetIO = "826b50da-1eb7-48f3-bd4b-d2350582c309"
GadgetUnits = "08630afb-492b-4c1a-9729-2a116101b53a"
HDF5 = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
UnitfulAstro = "6112ee07-acf9-5e0f-b108-d242c714bf9f"

[compat]
DataFrames = "~1.4.4"
GadgetIO = "~0.6.3"
GadgetUnits = "~0.2.3"
HDF5 = "~0.16.12"
Unitful = "~1.12.2"
UnitfulAstro = "~1.2.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "3e1d419c605125478b0b913ce649db8a9f3a21b0"

[[deps.ANSIColoredPrinters]]
git-tree-sha1 = "574baf8110975760d391c710b6341da1afa48d8c"
uuid = "a4c015fc-c6ff-483c-b24f-f7ea428134e9"
version = "0.0.1"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.CommonSolve]]
git-tree-sha1 = "9441451ee712d1aec22edad62db1a9af3dc8d852"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.3"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "00a2cccc7f098ff3b66806862d275ca3db9e6e5a"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.5.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.Cosmology]]
deps = ["QuadGK", "Unitful", "UnitfulAstro"]
git-tree-sha1 = "27996feb4391e8373417471feaa8638da773e589"
uuid = "76746363-e552-5dba-9a5a-cef6fa9cc5ab"
version = "1.0.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SnoopPrecompile", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d4f69885afa5e6149d0cab3818491565cf41446d"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.4.4"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Documenter]]
deps = ["ANSIColoredPrinters", "Base64", "Dates", "DocStringExtensions", "IOCapture", "InteractiveUtils", "JSON", "LibGit2", "Logging", "Markdown", "REPL", "Test", "Unicode"]
git-tree-sha1 = "6030186b00a38e9d0434518627426570aac2ef95"
uuid = "e30172f5-a6a5-5a46-863b-614d45cd2de4"
version = "0.27.23"

[[deps.DocumenterTools]]
deps = ["AbstractTrees", "Base64", "DocStringExtensions", "Documenter", "FileWatching", "Gumbo", "LibGit2", "OpenSSH_jll", "Sass"]
git-tree-sha1 = "d80421b418f6cc6321acd131011ea68e890f46f3"
uuid = "35a29f4d-8980-5a13-9543-d66fff28ecb8"
version = "0.1.16"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GadgetIO]]
deps = ["Dates", "DelimitedFiles", "Documenter", "DocumenterTools", "Downloads", "LinearAlgebra", "Printf", "ProgressMeter", "Test"]
git-tree-sha1 = "9ee97a68166313d5c7b9b3fe8eea02e17d3fa321"
uuid = "826b50da-1eb7-48f3-bd4b-d2350582c309"
version = "0.6.3"

[[deps.GadgetUnits]]
deps = ["Cosmology", "Documenter", "DocumenterTools", "Downloads", "GadgetIO", "Roots", "Test", "Unitful", "UnitfulAstro"]
git-tree-sha1 = "1a8885dfbe08bc3b8d27d6684c8acf8216f40fb3"
uuid = "08630afb-492b-4c1a-9729-2a116101b53a"
version = "0.2.4"

[[deps.Gumbo]]
deps = ["AbstractTrees", "Gumbo_jll", "Libdl"]
git-tree-sha1 = "a1a138dfbf9df5bace489c7a9d5196d6afdfa140"
uuid = "708ec375-b3d6-5a57-a7ce-8257bf98657a"
version = "0.8.2"

[[deps.Gumbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "29070dee9df18d9565276d68a596854b1764aa38"
uuid = "528830af-5a63-567c-a44a-034ed33b8444"
version = "0.10.2+0"

[[deps.HDF5]]
deps = ["Compat", "HDF5_jll", "Libdl", "Mmap", "Random", "Requires", "UUIDs"]
git-tree-sha1 = "b5df7c3cab3a00c33c2e09c6bd23982a75e2fbb2"
uuid = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
version = "0.16.13"

[[deps.HDF5_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "OpenSSL_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "4cc2bb72df6ff40b055295fdef6d92955f9dede8"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.12.2+2"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InvertedIndices]]
git-tree-sha1 = "82aec7a3dd64f4d9584659dc0b62ef7db2ef3e19"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.2.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenSSH_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "OpenSSL_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1b2f042897343a9dfdcc9366e4ecbd3d00780c49"
uuid = "9bd350c2-7e96-507f-8002-3f2e150b4e1b"
version = "8.9.0+1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6e9dba33f9f2c44e08a020b0caf6903be540004"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.19+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "8175fc2b118a3755113c8e68084dc1a9e63c61ee"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "LaTeXStrings", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "96f6db03ab535bdb901300f88335257b0018689d"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "de191bc385072cc6c7ed3ffdc1caeed3f22c74d4"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.7.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Roots]]
deps = ["ChainRulesCore", "CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "a3db467ce768343235032a1ca0830fc64158dadf"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.0.8"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Sass]]
deps = ["libsass_jll"]
git-tree-sha1 = "aa841c3738cec78b5dbccd56dda332710f35f6a5"
uuid = "322a6be2-4ae8-5d68-aaf1-3e960788d1d9"
version = "0.2.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d670a70dd3cdbe1c1186f2f17c9a68a7ec24838c"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.12.2"

[[deps.UnitfulAngles]]
deps = ["Dates", "Unitful"]
git-tree-sha1 = "d6cfdb6ddeb388af1aea38d2b9905fa014d92d98"
uuid = "6fb2a4bd-7999-5318-a3b2-8ad61056cd98"
version = "0.6.2"

[[deps.UnitfulAstro]]
deps = ["Unitful", "UnitfulAngles"]
git-tree-sha1 = "05adf5e3a3bd1038dd50ff6760cddd42380a7260"
uuid = "6112ee07-acf9-5e0f-b108-d242c714bf9f"
version = "1.2.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libsass_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "941afb93587dcec07f89e511057f5efc0bec6f0d"
uuid = "47bcb7c8-5119-555a-9eeb-0afcc36cd728"
version = "3.6.4+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═568ed730-9319-11ec-09ca-b1e962cadf22
# ╠═b20de40e-21fe-4d4e-b0dc-fcb0d4ead59e
# ╠═6647922f-bcc9-4438-a454-d8dba7a2d103
# ╠═f3cc08bd-b2c6-42a2-a757-5ef4e7c685e9
# ╠═a0f96fb1-8fda-410f-8861-563803f3798e
# ╠═108bb88d-0d89-4601-a2b2-3af50fcd1f3b
# ╠═6ea6a52d-1ee4-4069-a5ca-4b0021a7b69c
# ╠═39f2147d-ac35-4b14-ade6-ee07d17f84c6
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
