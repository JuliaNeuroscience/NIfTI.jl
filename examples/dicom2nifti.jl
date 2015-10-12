# dicom2nifti.jl
# Simple DICOM to NIfTI converter

# This script walks through the current working directory looking for DICOM
# files and sorts them into volumes. It then creates directories corresponding
# to each protocol, and saves nii files named by series number. It has
# (barely) been tested with DICOM files from a Siemens scanner. It may not
# work with other scanners. It may give incorrect information. Test thoroughly
# and use at your own risk.

using DICOM, NIfTI, ArgParse

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

function walk(fn::Function, path::String)
    for fname in readdir(path)
        if fname == "." || fname == ".."
            continue
        end
        fpath = joinpath(path, fname)
        if isdir(lstat(fpath))
            walk(fn, fpath)
        else
            fn(fpath)
        end
    end
end

function findtag(d, tag, from, what)
    tag = lookup(d, tag)
    if tag == false
        error("$what missing from $from")
    end
    tag.data
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
    walk(cmd["dicomdir"]) do fpath
        if !endswith(fpath, cmd["extension"])
            return
        end

        local d
        f = open(fpath, "r")
        d = DICOM.dcm_parse(f)

        series_number_index = findfirst((x) -> x.tag == (0x0020, 0x0011), d)
        if series_number_index == 0
            println("No series number found for $fpath; ignoring")
            return
        end
        series_number = d[series_number_index].data[1]

        dicom_arr = get(dicoms, series_number, nothing)
        if dicom_arr == nothing
            dicoms[series_number] = dicom_arr = {}
        end
        push!(dicom_arr, d)
    end

    # Sort DICOMs
    for (series_number, slices) in dicoms
        d = slices[1]

        protocol_name = findtag(d, (0x0018, 0x1030), series_number, "Protocol Name")[1]
        orientation = findtag(d, (0x0020, 0x0037), series_number,
            "Image Orientation (Patient)")
        pixel_spacing = findtag(d, (0x0028, 0x0030), series_number, "Pixel Spacing")
        slice_thickness = findtag(d, (0x0018, 0x0050), series_number, "Slice Thickness")[1]
        time_step = lookup(d, (0x0018, 0x0080)) # This is the TR, but may not be the time step
        phase_encoding_direction = lookup(d, (0x0018, 0x1312))

        # Convert orentation to RAS
        orientation = transform*reshape(orientation, (3, 2))

        # Determine Z as cross product of X and Y
        orientation = hcat(orientation, cross(orientation[:, 1], orientation[:, 2]))

        # Scale by pixel size
        orientation = orientation*[pixel_spacing[1] 0 0; 0 pixel_spacing[2] 0; 0 0 slice_thickness]

        # Find Z coordinates of each slice
        image_positions = transform*hcat(
            [float32(findtag(d, (0x0020, 0x0032), series_number, "Image Position (Patient)"))
            for d in slices]...)
        z_coords = [dot(image_positions[:, i], orientation[:, 3])
            for i = 1:size(image_positions, 2)]

        # Sort Z coordinates
        p = sortperm(z_coords)

        # Concatenate slices along Z axis according to increasing Z coordinate
        raw = cat(3,
            [findtag(slices[i], (0x7FE0, 0x0010), series_number, "Pixel Data")[1] for i in p]...)

        orientation = float32(hcat(orientation, image_positions[:, p[1]]))

        dim_info = phase_encoding_direction.data[1] == "ROW" ? (1, 2, 3) :
                phase_encoding_direction.data[1] == "COL" ? (2, 1, 3) :
                (0, 0, 0)

        voxel_size = tuple(pixel_spacing..., slice_thickness)

        # Permute according to axes specified on command line
        if axes != nothing
            # Determine permutation of current volume to RAS
            ras = Array(Int, 3)
            sign = Array(Bool, 3)
            rg = [1:3;]
            abs_orientation = abs(orientation)
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
            flip_sign = find(sign[specified_axes] .!= specified_sign)
            orientation[:, flip_sign] = -orientation[:, flip_sign]
            for dim in flip_sign
                raw = flipdim(raw, dim)
                orientation[dim, 4] += (sign[specified_axes[dim]] ? 1 : -1) * 
                    (size(raw, dim) - 1)*voxel_size[dim]
            end
        end

        # Create a directory for each protocol
        protocol_dir = joinpath(cmd["niftidir"], replace(protocol_name, "/", ""))
        if !isdir(protocol_dir)
            mkdir(protocol_dir)
        end

        # Write NIfTI volumes
        ni = NIVolume(raw; voxel_size=voxel_size,
            orientation=orientation, dim_info=dim_info,
            time_step=time_step != false && !isempty(time_step.data) ? time_step.data[1] : 0f0)
        niwrite(joinpath(protocol_dir, "$(series_number).nii"), ni)
    end
end

main()