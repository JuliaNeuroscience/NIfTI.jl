
module NIfTI

using CodecZlib, Mmap, MappedArrays, TranscodingStreams

import Base.getindex, Base.size, Base.ndims, Base.length, Base.write, Base64
export NIVolume, niread, niwrite, voxel_size, time_step, vox, getaffine, setaffine

include("parsers.jl")
include("header.jl")
include("extensions.jl")
include("volume.jl")

function string_tuple(x::String, n::Int)
    a = codeunits(x)
    padding = zeros(UInt8, n-length(a))
    (a..., padding...)
end
string_tuple(x::AbstractString) = string_tuple(bytestring(x))

struct NIVolume{T<:Number,N,R} <: AbstractArray{T,N}
    header::NIfTI1Header
    extensions::Vector{NIfTIExtension}
    raw::R
end

NIVolume(header::NIfTI1Header, extensions::Vector{NIfTIExtension}, raw::AbstractArray{T,N}) where {T<:Number,N} =
    NIVolume{typeof(one(T)*1f0+1f0),N,typeof(raw)}(header, extensions, raw)
NIVolume(header::NIfTI1Header, raw::AbstractArray{T,N}) where {T<:Number,N} =
    NIVolume{typeof(one(T)*1f0+1f0),N,typeof(raw)}(header, NIfTIExtension[], raw)
NIVolume(header::NIfTI1Header, extensions::Vector{NIfTIExtension}, raw::AbstractArray{Bool,N}) where {N} =
    NIVolume{Bool,N,typeof(raw)}(header, extensions, raw)
NIVolume(header::NIfTI1Header, raw::AbstractArray{Bool,N}) where {N} =
    NIVolume{Bool,N,typeof(raw)}(header, NIfTIExtension[], raw)



# Always in mm
function voxel_size(hdr)
    [hdr.pixdim[i] * SPATIAL_UNIT_MULTIPLIERS[hdr.xyzt_units & Int8(3)] for i = 2:min(hdr.dim[1], 3)+1]
end

# Always in ms
time_step(hdr) = hdr.pixdim[5] * TIME_UNIT_MULTIPLIERS[hdr.xyzt_units >> 3]

function to_dim_info(dim_info::Tuple{Integer,Integer,Integer})
    if dim_info[1] > 3 || dim_info[1] < 0
        error("Invalid frequency dimension $(dim_info[1])")
    elseif dim_info[2] > 3 || dim_info[2] < 0
        error("Invalid phase dimension $(dim_info[2])")
    elseif dim_info[3] > 3 || dim_info[3] < 0
        error("Invalid slice dimension $(dim_info[3])")
    end

    return Int8(dim_info[1] | (dim_info[2] << 2) | (dim_info[3] << 4))
end

# Returns or sets dim_info as a tuple whose values are the frequency, phase, and slice dimensions
function dim_info(header::NIfTI1Header)
    return (
        header.dim_info & int8(3),
        (header.dim_info >> 2) & int8(3),
        (header.dim_info >> 4) & int8(3)
    )
end
function dim_info(header::NIfTI1Header, dim_info::Tuple{T, T, T}) where {T<:Integer}
    header.dim_info = to_dim_info(dim_info)
end

# Gets the size of a type in bits
nibitpix(t::Type) = Int16(sizeof(t)*8)
nibitpix(::Type{Bool}) = Int16(1)

# Convert a NIfTI header to a 4x4 affine transformation matrix
function getaffine(h::NIfTI1Header)
    pixdim = convert(Tuple{Vararg{Float64}}, h.pixdim)
    # See documentation at
    # http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/qsform.html
    if h.sform_code > 0
        # Method 3
        return Float64[
            collect(h.srow_x)'
            collect(h.srow_y)'
            collect(h.srow_z)'
            0 0 0 1
        ]
    elseif h.qform_code > 0
        # Method 2
        # http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/quatern.html
        b::Float64 = h.quatern_b
        c::Float64 = h.quatern_c
        d::Float64 = h.quatern_d
        a = sqrt(1 - b*b - c*c - d*d)
        A = [
            a*a+b*b-c*c-d*d   2*b*c-2*a*d       2*b*d+2*a*c
            2*b*c+2*a*d       a*a+c*c-b*b-d*d   2*c*d-2*a*b
            2*b*d-2*a*c       2*c*d+2*a*b       a*a+d*d-c*c-b*b
        ]
        B = [
            pixdim[2]
            pixdim[3]
            pixdim[1]*pixdim[4]
        ]
        return vcat(hcat(A*B, [
            h.qoffset_x
            h.qoffset_y
            h.qoffset_z
        ]), [0 0 0 1])
    elseif h.qform_code == 0
        # Method 1
        return Float64[
             pixdim[2] 0         0          0 0
                       0 pixdim[3]          0 0
                       0         0  bitpix[4] 0
                       0         0          0 1
        ]
    end
end

# Set affine matrix of NIfTI header
function setaffine(h::NIfTI1Header, affine::Array{T,2}) where {T}
    size(affine, 1) == size(affine, 2) == 4 ||
        error("affine matrix must be 4x4")
    affine[4, 1] == affine[4, 2] == affine[4, 3] == 0 && affine[4, 4] == 1 ||
        error("last row of affine matrix must be [0 0 0 1]")
    h.qform_code = one(Int16)
    h.sform_code = one(Int16)
    h.pixdim = (zero(Float32), h.pixdim[2:end]...,)
    h.quatern_b = zero(Float32)
    h.quatern_c = zero(Float32)
    h.quatern_d = zero(Float32)
    h.qoffset_x = zero(Float32)
    h.qoffset_y = zero(Float32)
    h.qoffset_z = zero(Float32)
    h.srow_x = convert(Tuple{Vararg{Float32}}, (affine[1, :]...,))
    h.srow_y = convert(Tuple{Vararg{Float32}}, (affine[2, :]...,))
    h.srow_z = convert(Tuple{Vararg{Float32}}, (affine[3, :]...,))
    h
end

# Constructor
function NIVolume(
    # Optional MRI volume; if not given, an empty volume is used
    raw::AbstractArray{T}=Int16[],
    extensions::Union{Vector{NIfTIExtension},Nothing}=nothing;

    # Fields specified as UNUSED in NIfTI1 spec
    data_type::AbstractString="", db_name::AbstractString="", extents::Integer=Int32(0),
    session_error::Integer=Int16(0), regular::Integer=Int8(0), glmax::Integer=Int32(0),
    glmin::Integer=Int16(0),

    # The frequency encoding, phase encoding and slice dimensions.
    dim_info::NTuple{3, Integer}=(0, 0, 0),
    # Describes data contained in the file; for valid values, see
    # http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/group__NIfTI1__INTENT__CODES.html
    intent_p1::Real=0f0, intent_p2::Real=0f0, intent_p3::Real=0f0,
    intent_code::Integer=Int16(0), intent_name::AbstractString="",
    # Information about which slices were acquired, in case the volume has been padded
    slice_start::Integer=Int16(0), slice_end::Integer=Int16(0), slice_code::Int8=Int8(0),
    # The size of each voxel and the time step. These are formulated in mm unless xyzt_units is
    # explicitly specified
    voxel_size::NTuple{3, Real}=(1f0, 1f0, 1f0), time_step::Real=0f0, xyzt_units::Int8=Int8(18),
    # Slope and intercept by which volume shoudl be scaled
    scl_slope::Real=1f0, scl_inter::Real=0f0,
    # These describe how data should be scaled when displayed on the screen. They are probably
    # rarely used
    cal_max::Real=0f0, cal_min::Real=0f0,
    # The amount of time taken to acquire a slice
    slice_duration::Real=0f0,
    # Indicates a non-zero start point for time axis
    toffset::Real=0f0,

    # "Any text you like"
    descrip::AbstractString="",
    # Name of auxiliary file
    aux_file::AbstractString="",

    # Parameters for Method 2. See the NIfTI spec
    qfac::Float32=0f0, quatern_b::Real=0f0, quatern_c::Real=0f0, quatern_d::Real=0f0,
    qoffset_x::Real=0f0, qoffset_y::Real=0f0, qoffset_z::Real=0f0,
    # Orientation matrix for Method 3
    orientation::Union{Matrix{Float32},Nothing}=nothing
) where {T <: Number}

    local t
    if isempty(raw)
        raw = Int16[]
        t = Int16
    else
        t = T
    end

    if extensions == nothing
        extensions = NIfTIExtension[]
    end

    method2 = qfac != 0 || quatern_b != 0 || quatern_c != 0 || quatern_d != 0 || qoffset_x != 0 ||
        qoffset_y != 0 || qoffset_z != 0
    method3 = orientation != nothing

    if method2 && method3
        error("Orientation parameters for Method 2 and Method 3 are mutually exclusive")
    end

    if method3
        if size(orientation) != (3, 4)
            error("Orientation matrix must be of dimensions (3, 4)")
        end
    else
        orientation = zeros(Float32, 3, 4)
    end

    if slice_start == 0 && slice_end == 0 && dim_info[3] != 0
        slice_start = 0
        slice_end = size(raw, dim_info[3]) - 1
    end

    if isempty(extensions)
        vox_offset = 352
    else
        vox_offset = mapreduce(esize, +, extensions) + 352
    end

    NIVolume(
        NIfTI1Header(
            string_tuple(data_type, 10),
            string_tuple(db_name, 18),
            extents,
            session_error,
            regular,
            to_dim_info(dim_info),
            to_dim_i16(size(raw)),
            intent_p1,
            intent_p2,
            intent_p3,
            intent_code,
            eltype_to_i16(t),
            nibitpix(t),
            slice_start,
            (qfac, voxel_size..., time_step, 0, 0, 0),
            vox_offset,
            scl_slope, scl_inter,
            slice_end,
            slice_code,
            xyzt_units,
            cal_max,
            cal_min,
            slice_duration,
            toffset,
            glmax, glmin,
            string_tuple(descrip, 80),
            string_tuple(aux_file, 24),
            (method2 || method3),
            method3,
            quatern_b,
            quatern_c,
            quatern_d,
            qoffset_x,
            qoffset_y,
            qoffset_z,
            (orientation[1, :]...,),
            (orientation[2, :]...,),
            (orientation[3, :]...,),
            string_tuple(intent_name, 16), NP1_MAGIC),
        extensions,
        raw
    )
end

# Validates the header of a volume and updates it to match the volume's contents
function niupdate(vol::NIVolume{T}) where {T}
    vol.header.dim = nidim(vol.raw)
    vol.header.datatype = eltype_to_int16(T)
    vol.header.bitpix = nibitpix(T)
    vol.header.vox_offset = isempty(vol.extensions) ? Int32(352) :
        Int32(mapreduce(esize, +, vol.extensions) + SIZEOF_HDR1)
    vol
end

# Avoid method ambiguity
function Base.write(io::Base64.Base64EncodePipe, vol::NIVolume{UInt8,1})
    return invoke(write, (IO, NIVolume{UInt8,1}), io, vol)
end

# Convenience function to write a NIfTI file given a path
function niwrite(path::AbstractString, vol::NIVolume)
    if split(path,".")[end] == "gz"
        io = GzipCompressorStream(open(path, "w"))
        niwrite(io, vol)
        close(io)
    else
        io = open(path, "w")
        niwrite(io, vol)
        close(io)
    end
end

function niwrite(io, vol::NIVolume)
    write_header(io, vol.header)
    write_extension(io, vol.extensions)
    write_volume(io, vol.raw)
end

function niread(file::AbstractString; mmap::Bool=false, mode::AbstractString="r")
    io, hdr, swapped = niopen(file, mode)
    ex = read_extensions(io, hdr.vox_offset - (sizeof_hdr(hdr) + 4))

    if is_volume_separate(hdr.magic)
        vol = read_volume(hdr_to_img(file), mode, hdr, mmap)
    else
        vol = read_volume(io, hdr, mmap)
    end

    if swapped && sizeof(eltype(vol)) > 1
        return NIVolume(hdr, ex, mappedarray(ntoh, hton, vol))
    else
        return NIVolume(hdr, ex, vol)
    end
end

# Allow file to be indexed like an array, but with indices yielding scaled data
@inline function getindex(f::NIVolume{T}, idx::Vararg{Int}) where {T}
    return getindex(f.raw, idx...,) * f.header.scl_slope + f.header.scl_inter
end

add1(x::Union{AbstractArray{T},T}) where {T<:Integer} = x + 1
add1(::Colon) = Colon()
@inline vox(f::NIVolume, args...,) = getindex(f, map(add1, args)...,)
size(f::NIVolume) = size(f.raw)
size(f::NIVolume, d) = size(f.raw, d)
ndims(f::NIVolume) = ndims(f.raw)
length(f::NIVolume) = length(f.raw)
lastindex(f::NIVolume) = lastindex(f.raw)

end

