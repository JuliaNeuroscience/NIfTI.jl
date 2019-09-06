#!/bin/bash
# -*- mode: julia -*-
#=
exec julia --color=yes --startup-file=no -e "include(popfirst!(ARGS))" "${BASH_SOURCE[0]}" "$@"
=#

# dicom2nifti.jl
# Simple DICOM to NIfTI converter

# This script walks through the current working directory looking for DICOM
# files and sorts them into volumes. It then creates directories corresponding
# to each protocol, and saves nii files named by series number. It has
# (barely) been tested with DICOM files from a Siemens scanner. It may not
# work with other scanners. It may give incorrect information. Test thoroughly
# and use at your own risk.

using DICOM, NIfTI, ArgParse, LinearAlgebra

function parse_commandline()
    s = ArgParseSettings()
    s.description = "DICOM to NIfTI converter"

    @add_arg_table s begin
        "--extension"
            help = "assume files with given extension are DICOM files. Set to blank for no extension."
            default = ".dcm"
        "--sphinx"
            help = "correct for sphinx position (for monkeys)"
            action = :store_true
        "--axes"
            help = """orientation of tensor in the resulting volume. This is a three-character label
            that determines how the data is stored, and how it is displayed in applications that do
            not make use of the orientation matrix in the NIfTI header. Valid characters are:
            L (= left), R (= right), A (= anterior), P (= posterior), I (= inferior), and
            S (= superior). If not specified, the orientation in the DICOM file is preserved."""
            default = nothing
        "dicomdir"
            help = "directory of DICOM files"
            required = true
        "niftidir"
            help = "directory in which to put NIfTI files"
            required = true
    end

    parse_args(s)
end

function main()
    cmd = parse_commandline()

    if cmd["axes"] == nothing
        axes = nothing
    else
        axes = uppercase(cmd["axes"])
        if length(axes) != 3 || (!('R' in axes) && !('L' in axes)) ||
                (!('A' in axes) && !('P' in axes)) ||
                (!('S' in axes) && !('I' in axes))
            error("Invalid axes provided")
        end
    end

    transform = cmd["sphinx"] ? [1 0 0; 0 0 1; 0 -1 0] : [-1 0 0; 0 -1 0; 0 0 1]

    # Load all DICOMs into an array, indexed by series number
    dicoms = Dict{Int, Any}()
    for (root, dirs, files) in walkdir(cmd["dicomdir"])
        for file in files
            fpath = joinpath(root, file)
            if !endswith(fpath, cmd["extension"])
                return
            end

            local d
            println(fpath)
            d = DICOM.dcm_parse(fpath)

            if !(tag"Series Number" in keys(d))
                println("No series number found for $fpath; ignoring")
                return
            end
            series_number = d[tag"Series Number"]

            dicom_arr = get(dicoms, series_number, nothing)
            if dicom_arr == nothing
                dicoms[series_number] = dicom_arr = []
            end
            push!(dicom_arr, d)
        end
    end

    # Sort DICOMs
    for (series_number, slices) in dicoms
        d = slices[1]

        protocol_name = d[tag"Protocol Name"]
        orientation = d[tag"Image Orientation (Patient)"]
        pixel_spacing = d[tag"Pixel Spacing"]
        slice_spacing = slice_thickness = d[tag"Slice Thickness"]
        if tag"Spacing Between Slices" in keys(d)
            slice_spacing = d[tag"Spacing Between Slices"]
        end
        time_step = false
        if tag"Repetition Time" in keys(d)
            time_step = d[tag"Repetition Time"] # This is the TR, but may not be the time step
        end
        phase_encoding_direction = d[tag"In-plane Phase Encoding Direction"]
        pixel_rows = d[tag"Rows"]
        pixel_cols = d[tag"Columns"]

        # Convert orientation to RAS
        orientation = transform*reshape(orientation, (3, 2))

        # Determine Z as cross product of X and Y
        orientation = hcat(orientation, cross(orientation[:, 1], orientation[:, 2]))

        # Scale by pixel size
        orientation = orientation*[pixel_spacing[1] 0 0; 0 pixel_spacing[2] 0; 0 0 slice_spacing]

        # Find Z coordinates of each slice
        image_positions = transform * hcat([d[tag"Image Position (Patient)"] for d in slices]...)
        z_coords = [dot(image_positions[:, i], orientation[:, 3])
            for i = 1:size(image_positions, 2)]

        # Sort Z coordinates
        p = sortperm(z_coords)

        # Concatenate slices along Z axis according to increasing Z coordinate
        raw = cat([slices[i][tag"Pixel Data"] for i in p]..., dims = 3)

        orientation = hcat(orientation, image_positions[:, p[1]])

        dim_info = phase_encoding_direction == "ROW" ? (1, 2, 3) :
                phase_encoding_direction == "COL" ? (2, 1, 3) :
                (0, 0, 0)

        voxel_size = tuple(pixel_spacing..., slice_spacing)

        # Permute according to axes specified on command line
        if axes != nothing
            # Determine permutation of current volume to RAS
            ras = Array{Int,1}(undef, 3)
            sign = Array{Bool,1}(undef, 3)
            rg = [1:3;]
            abs_orientation = abs.(orientation)
            for i = 1:3
                idx = findmax(abs_orientation[i, rg])[2]
                ras[i] = splice!(rg, idx)
                sign[i] = orientation[i, ras[i]] >= 0
            end

            # Determine permutation of RAS to user-preferred orientation
            specified_axes = Int[(c in "RL" ? 1 :
                                  c in "AP" ? 2 : 3)
                                 for c in axes]
            specified_sign = Bool[c in "RAS" for c in axes]
            output_axes = ras[specified_axes]

            # Permute
            raw = permutedims(raw, output_axes)
            orientation[:, 1:3] = orientation[:, output_axes]
            voxel_size = voxel_size[output_axes]
            dim_info = dim_info[output_axes]

            # Flip signs and dimensions
            flip_sign = findall(sign[specified_axes] .!= specified_sign)
            orientation[:, flip_sign] = -orientation[:, flip_sign]
            for dim in flip_sign
                raw = reverse(raw, dims=dim)
                orientation[dim, 4] += (sign[specified_axes[dim]] ? 1 : -1) *
                    (size(raw, dim) - 1)*voxel_size[dim]
            end
        end

        # Create a directory for each protocol
        protocol_dir = joinpath(cmd["niftidir"], replace(protocol_name, "/" => ""))

        if !isdir(cmd["niftidir"])
            mkdir(cmd["niftidir"])
        end
        if !isdir(protocol_dir)
            mkdir(protocol_dir)
        end

        # Write NIfTI volumes
        # convert to Float32 as mricron and fsleyes apparently do not support UInt16
        ni = NIVolume(Array{Float32, 3}(raw); voxel_size=voxel_size,
            orientation=convert(Array{Float32, 2}, orientation), dim_info=dim_info,
            time_step=time_step != false ? time_step : 0f0)
        niwrite(joinpath(protocol_dir, "$(series_number).nii"), ni)
    end
end

main()
