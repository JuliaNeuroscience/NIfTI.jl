module NIfTI

using CodecZlib, Mmap, MappedArrays, TranscodingStreams

import Base.getindex, Base.size, Base.ndims, Base.length, Base.write, Base64
export NIVolume, niread, niwrite, voxel_size, time_step, vox, getaffine, setaffine

include("parsers.jl")
include("extensions.jl")
include("volume.jl")
include("headers.jl")

"""
    NIVolume{T<:Number,N,R} <: AbstractArray{T,N}
An `N`-dimensional NIfTI volume, with raw data of type 
`R`. Note that if `R <: Number`, it will be converted to `Float32`. Additionally, the header is automatically
updated to be consistent with the raw volume. 

# Members
- `header`: a `NIfTIHeader`
- `extensions`: a Vector of `NIfTIExtension`s 
- `raw`: Raw data of type `R`
"""
struct NIVolume{T<:Number,N,R} <: AbstractArray{T,N}
    header::NIfTIHeader
    extensions::Vector{NIfTIExtension}
    raw::R
end
NIVolume(header::NIfTIHeader, extensions::Vector{NIfTIExtension}, raw::R) where {R} =
    niupdate(new(header, extensions, raw))

NIVolume(header::NIfTIHeader, extensions::Vector{NIfTIExtension}, raw::AbstractArray{T,N}) where {T<:Number,N} =
    NIVolume{typeof(one(T) * 1.0f0 + 1.0f0),N,typeof(raw)}(header, extensions, raw)
NIVolume(header::NIfTIHeader, raw::AbstractArray{T,N}) where {T<:Number,N} =
    NIVolume{typeof(one(T) * 1.0f0 + 1.0f0),N,typeof(raw)}(header, NIfTIExtension[], raw)

NIVolume(header::NIfTIHeader, extensions::Vector{NIfTIExtension}, raw::AbstractArray{Bool,N}) where {N} =
    NIVolume{Bool,N,typeof(raw)}(header, extensions, raw)
NIVolume(header::NIfTIHeader, raw::AbstractArray{Bool,N}) where {N} =
    NIVolume{Bool,N,typeof(raw)}(header, NIfTIExtension[], raw)

header(x::NIVolume) = getfield(x, :header)

include("coordinates.jl")

# Always in ms
"""
    time_step(header::NIfTIHeader)

Get the TR **in ms** from a `NIfTIHeader`.
"""
time_step(header::NIfTIHeader) =
    header.pixdim[5] * TIME_UNIT_MULTIPLIERS[header.xyzt_units>>3]

# Gets the size of a type in bits
nibitpix(t::Type) = Int16(sizeof(t) * 8)
nibitpix(::Type{Bool}) = Int16(1)

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
    dim_info::NTuple{3,Integer}=(0, 0, 0),
    # Describes data contained in the file; for valid values, see
    # http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/group__NIfTI1__INTENT__CODES.html
    intent_p1::Real=0.0f0, intent_p2::Real=0.0f0, intent_p3::Real=0.0f0,
    intent_code::Integer=Int16(0), intent_name::AbstractString="",
    # Information about which slices were acquired, in case the volume has been padded
    slice_start::Integer=Int16(0), slice_end::Integer=Int16(0), slice_code=UInt8(0),
    # The size of each voxel and the time step. These are formulated in mm unless xyzt_units is
    # explicitly specified
    voxel_size::NTuple{3,Real}=(1.0f0, 1.0f0, 1.0f0), time_step::Real=0.0f0, xyzt_units=UInt8(18),
    # Slope and intercept by which volume shoudl be scaled
    scl_slope::Real=1.0f0, scl_inter::Real=0.0f0,
    # These describe how data should be scaled when displayed on the screen. They are probably
    # rarely used
    cal_max::Real=0.0f0, cal_min::Real=0.0f0,
    # The amount of time taken to acquire a slice
    slice_duration::Real=0.0f0,
    # Indicates a non-zero start point for time axis
    toffset::Real=0.0f0,

    # "Any text you like"
    descrip::AbstractString="",
    # Name of auxiliary file
    aux_file::AbstractString="",

    # Parameters for Method 2. See the NIfTI spec
    qfac::Float32=0.0f0, quatern_b::Real=0.0f0, quatern_c::Real=0.0f0, quatern_d::Real=0.0f0,
    qoffset_x::Real=0.0f0, qoffset_y::Real=0.0f0, qoffset_z::Real=0.0f0,
    # Orientation matrix for Method 3
    orientation::Union{Matrix{Float32},Nothing}=nothing) where {T<:Number}
    local t
    if isempty(raw)
        raw = Int16[]
        t = Int16
    else
        t = T
    end

    if extensions === nothing
        extensions = NIfTIExtension[]
    end

    method2 = qfac != 0 || quatern_b != 0 || quatern_c != 0 || quatern_d != 0 || qoffset_x != 0 ||
              qoffset_y != 0 || qoffset_z != 0
    method3 = orientation !== nothing

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

    NIVolume(NIfTI1Header(SIZEOF_HDR1, string_tuple(data_type, 10), string_tuple(db_name, 18), extents, session_error,
            regular, to_dim_info(dim_info), to_dim_i16(size(raw)), intent_p1, intent_p2,
            intent_p3, intent_code, eltype_to_int16(t), nibitpix(t),
            slice_start, (qfac, voxel_size..., time_step, 0, 0, 0), 352,
            scl_slope,
            scl_inter,
            slice_end,
            UInt8(slice_code),
            UInt8(xyzt_units),
            cal_max, cal_min, slice_duration,
            toffset, glmax, glmin, string_tuple(descrip, 80), string_tuple(aux_file, 24), (method2 || method3),
            method3, quatern_b, quatern_c, quatern_d,
            qoffset_x, qoffset_y, qoffset_z, (orientation[1, :]...,),
            (orientation[2, :]...,), (orientation[3, :]...,), string_tuple(intent_name, 16), NP1_MAGIC), extensions, raw)
end

# Validates the header of a volume and updates it to match the volume's contents
function niupdate(vol::NIVolume{T,N,R}) where {T,N,R}
    vol.header.dim = to_dim_i16(size(vol.raw))
    vol.header.datatype = eltype_to_int16(eltype(R))
    vol.header.bitpix = nibitpix(T)
    vol.header.vox_offset = isempty(vol.extensions) ? vol.header.sizeof_hdr + 4 :
                            Int32(mapreduce(esize, +, vol.extensions) + vol.header.sizeof_hdr)
    vol
end

# Avoid method ambiguity
write(io::Base64.Base64EncodePipe, vol::NIVolume{UInt8,1}) = invoke(write, (IO, NIVolume{UInt8,1}), io, vol)

# Write a NIfTI file
function write(io::IO, vol::NIVolume)
    write(io, niupdate(vol).header)
    if isempty(vol.extensions)
        write(io, Int32(0))
    else
        for ex in vol.extensions
            sz = esize(ex)
            write(io, Int32(sz))
            write(io, Int32(ex.ecode))
            write(io, ex.edata)
            write(io, zeros(UInt8, sz - length(ex.edata)))
        end
    end
    if eltype(vol.raw) == Bool
        write(io, BitArray(vol.raw))
    else
        write(io, vol.raw)
    end
end

"""
    niwrite(path::AbstractString, vol::NIVolume)   

Write a NIVolume to a file specified by `path`.
"""
function niwrite(path::AbstractString, vol::NIVolume)
    if split(path, ".")[end] == "gz"
        io = open(path, "w")
        stream = GzipCompressorStream(io)
        write(stream, vol)
        close(stream)
        close(io)
    else
        io = open(path, "w")
        write(io, vol)
        close(io)
    end
end

# Read header from a NIfTI file
function read_header(io::IO)
    header_size = peek(io, Int32)
    if header_size == SIZEOF_HDR1 || bswap(header_size) == SIZEOF_HDR1
        header, swapped = read(io, NIfTI1Header)
        #error("This is not a NIfTI-1 file")
        return header, swapped
    end

    if header_size == SIZEOF_HDR2 || bswap(header_size) == SIZEOF_HDR2
        header, swapped = read(io, NIfTI2Header)
        #error("This is not a NIfTI-1 file")
        return header, swapped
    end
    error("This is not a NIfTI file")
end

# Look for a gzip header in an IOStream
function isgz(io::IO)
    try
        ret = read(io, UInt8) == 0x1F && read(io, UInt8) == 0x8B
        seekstart(io)
        ret
    catch err
        if isa(err, EOFError)
            @debug "reading the file resulted in an EOF error and \nthe end of the file was read.\nNo more data was available to read from the filestream.\nIt is likely that the file was corrupted or is empty (0 bytes)."
            rethrow(err)
        end
    end
end

"""
    niread(file; mmap=false, mode="r")

Read a NIfTI file to a NIVolume. Set `mmap=true` to memory map the volume.
"""
function niread(file::AbstractString; mmap::Bool=false, mode::AbstractString="r")
    io = niopen(file, mode)
    hdr, swapped = read_header(io)
    ex = read_extensions(io, hdr.vox_offset - hdr.sizeof_hdr - 4)

    if hdr.magic === NP1_MAGIC || hdr.magic == NP2_MAGIC
        vol = read_volume(io, to_eltype(hdr.datatype), to_dimensions(hdr.dim), mmap)
    else
        vol = read_volume(niopen(hdr_to_img(file), mode), to_eltype(hdr.datatype), to_dimensions(hdr.dim), mmap)
    end

    if swapped && sizeof(eltype(vol)) > 1
        return NIVolume(hdr, ex, mappedarray(ntoh, hton, vol))
    else
        return NIVolume(hdr, ex, vol)
    end
end

# Allow file to be indexed like an array, but with indices yielding scaled data
@inline getindex(f::NIVolume{T}, idx::Vararg{Int}) where {T} =
    f.header.scl_slope == zero(typeof(f.header.scl_slope)) ?
    getindex(f.raw, idx...,) :
    getindex(f.raw, idx...,) * f.header.scl_slope + f.header.scl_inter

add1(x::Union{AbstractArray{T},T}) where {T<:Integer} = x + 1
add1(::Colon) = Colon()

"""
    vox(f::NIVolume, args...,)

Get the value of a voxel from volume `f`, scaled by slope and intercept given in header, with 0-based indexing.
`args` are the voxel indices and the length should be the number of dimensions in `f`.
"""
@inline vox(f::NIVolume, args...,) = getindex(f, map(add1, args)...,)

size(f::NIVolume) = size(f.raw)
size(f::NIVolume, d) = size(f.raw, d)
ndims(f::NIVolume) = ndims(f.raw)
length(f::NIVolume) = length(f.raw)
lastindex(f::NIVolume) = lastindex(f.raw)

end
