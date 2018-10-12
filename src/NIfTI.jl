# NIfTI.jl
# Methods for reading NIfTI MRI files in Julia

module NIfTI

using GZip, Mmap

import Base.getindex, Base.size, Base.ndims, Base.length, Base.write, Base64
export NIVolume, niread, niwrite, voxel_size, time_step, vox, getaffine, setaffine

function define_packed(ty::DataType)
    packed_offsets = cumsum([sizeof(x) for x in ty.types])
    sz = pop!(packed_offsets)
    pushfirst!(packed_offsets, 0)

    @eval begin
        function Base.read(io::IO, ::Type{$ty})
            bytes = read!(io, Array{UInt8}(undef, $sz...))
            hdr = $(Expr(:new, ty, [:(unsafe_load(convert(Ptr{$(ty.types[i])}, pointer(bytes)+$(packed_offsets[i])))) for i = 1:length(packed_offsets)]...,))
            if hdr.sizeof_hdr == ntoh(Int32(348))
                return byteswap(hdr), true
            end
            hdr, false
        end
        function Base.write(io::IO, x::$ty)
            bytes = UInt8[]
            for name in fieldnames($ty)
                append!(bytes, reinterpret(UInt8, [getfield(x,name)]))
            end
            write(io, bytes)
            $sz
        end
    end
    nothing
end

mutable struct NIfTI1Header
    sizeof_hdr::Int32

    data_type::NTuple{10,UInt8}
    db_name::NTuple{18,UInt8}
    extents::Int32
    session_error::Int16
    regular::Int8

    dim_info::Int8
    dim::NTuple{8,Int16}
    intent_p1::Float32
    intent_p2::Float32
    intent_p3::Float32
    intent_code::Int16
    datatype::Int16
    bitpix::Int16
    slice_start::Int16
    pixdim::NTuple{8,Float32}
    vox_offset::Float32
    scl_slope::Float32
    scl_inter::Float32
    slice_end::Int16
    slice_code::Int8
    xyzt_units::Int8
    cal_max::Float32
    cal_min::Float32
    slice_duration::Float32
    toffset::Float32

    glmax::Int32
    glmin::Int32

    descrip::NTuple{80,UInt8}
    aux_file::NTuple{24,UInt8}

    qform_code::Int16
    sform_code::Int16
    quatern_b::Float32
    quatern_c::Float32
    quatern_d::Float32
    qoffset_x::Float32
    qoffset_y::Float32
    qoffset_z::Float32

    srow_x::NTuple{4,Float32}
    srow_y::NTuple{4,Float32}
    srow_z::NTuple{4,Float32}

    intent_name::NTuple{16,UInt8}

    magic::NTuple{4,UInt8}
end
define_packed(NIfTI1Header)

# byteswapping

function byteswap(hdr::NIfTI1Header)
    for fn in fieldnames(typeof(hdr))
        val = getfield(hdr, fn)
        if isa(val, Number) && sizeof(val) > 1
            setfield!(hdr, fn, ntoh(val))
        elseif isa(val, NTuple) && sizeof(eltype(val)) > 1
            setfield!(hdr, fn, map(ntoh, val))
        end
    end
    hdr
end

const SIZEOF_HDR = Int32(348)

const NIfTI_DT_BITSTYPES = Dict{Int16,Type}([
    (Int16(2), UInt8),
    (Int16(4), Int16),
    (Int16(8), Int32),
    (Int16(16), Float32),
    (Int16(32), ComplexF32),
    (Int16(64), Float64),
    (Int16(256), Int8),
    (Int16(512), UInt16),
    (Int16(768), UInt32),
    (Int16(1024), Int64),
    (Int16(1280), UInt64),
    (Int16(1792), ComplexF64)
])
const NIfTI_DT_BITSTYPES_REVERSE = Dict{Type,Int16}()
for (k, v) in NIfTI_DT_BITSTYPES
    NIfTI_DT_BITSTYPES_REVERSE[v] = k
end

const NP1_MAGIC = (0x6e,0x2b,0x31,0x00)
const NI1_MAGIC = (0x6e,0x69,0x31,0x00)

function string_tuple(x::String, n::Int)
    a = codeunits(x)
    padding = zeros(UInt8, n-length(a))
    (a..., padding...)
end
string_tuple(x::AbstractString) = string_tuple(bytestring(x))

mutable struct NIfTI1Extension
    ecode::Int32
    edata::Vector{UInt8}
end

mutable struct NIVolume{T<:Number,N,R} <: AbstractArray{T,N}
    header::NIfTI1Header
    extensions::Vector{NIfTI1Extension}
    raw::R

end
NIVolume(header::NIfTI1Header, extensions::Vector{NIfTI1Extension}, raw::R) where {R}=
    niupdate(new(header, extensions, raw))

NIVolume(header::NIfTI1Header, extensions::Vector{NIfTI1Extension}, raw::AbstractArray{T,N}) where {T<:Number,N} =
    NIVolume{typeof(one(T)*1f0+1f0),N,typeof(raw)}(header, extensions, raw)
NIVolume(header::NIfTI1Header, raw::AbstractArray{T,N}) where {T<:Number,N} =
    NIVolume{typeof(one(T)*1f0+1f0),N,typeof(raw)}(header, NIfTI1Extension[], raw)

# Conversion factors to mm/ms
# http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/xyzt_units.html
const SPATIAL_UNIT_MULTIPLIERS = [
    1000,   # 1 => NIfTI_UNITS_METER
    1,      # 2 => NIfTI_UNITS_MM
    0.001   # 3 => NIfTI_UNITS_MICRON
]
const TIME_UNIT_MULTIPLIERS = [
    1000,   # NIfTI_UNITS_SEC
    1,      # NIfTI_UNITS_MSEC
    0.001,  # NIfTI_UNITS_USEC
    1,      # NIfTI_UNITS_HZ
    1,      # NIfTI_UNITS_PPM
    1       # NIfTI_UNITS_RADS
]

# Always in mm
voxel_size(header::NIfTI1Header) =
    [header.pixdim[i] * SPATIAL_UNIT_MULTIPLIERS[header.xyzt_units & Int8(3)] for i = 2:min(header.dim[1], 3)+1]

# Always in ms
time_step(header::NIfTI1Header) =
    header.pixdim[5] * TIME_UNIT_MULTIPLIERS[header.xyzt_units >> 3]

function to_dim_info(dim_info::Tuple{Integer,Integer,Integer})
    if dim_info[1] > 3 || dim_info[1] < 0
        error("Invalid frequency dimension $(dim_info[1])")
    elseif dim_info[2] > 3 || dim_info[2] < 0
        error("Invalid phase dimension $(dim_info[2])")
    elseif dim_info[3] > 3 || dim_info[3] < 0
        error("Invalid slice dimension $(dim_info[3])")
    end

    Int8(dim_info[1] | (dim_info[2] << 2) | (dim_info[3] << 4))
end

# Returns or sets dim_info as a tuple whose values are the frequency, phase, and slice dimensions
dim_info(header::NIfTI1Header) = (header.dim_info & int8(3), (header.dim_info >> 2) & int8(3),
    (header.dim_info >> 4) & int8(3))
dim_info(header::NIfTI1Header, dim_info::Tuple{T, T, T}) where {T<:Integer} =
    header.dim_info = to_dim_info(dim_info)

# Gets dim to be used in header
function nidim(x::AbstractArray)
    dim = ones(Int16, 8)
    dim[1] = ndims(x)
    dim[2:dim[1]+1] = [size(x)...]
    (dim...,)
end

# Gets datatype to be used in header
function nidatatype(t::Type)
    t = get(NIfTI_DT_BITSTYPES_REVERSE, t, nothing)
    if t == nothing
        error("Unsupported data type $T")
    end
    t
end

# Gets the size of a type in bits
nibitpix(t::Type) = Int16(sizeof(t)*8)

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
    extensions::Union{Vector{NIfTI1Extension},Nothing}=nothing;

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
    orientation::Union{Matrix{Float32},Nothing}=nothing) where {T <: Number}
    local t
    if isempty(raw)
        raw = Int16[]
        t = Int16
    else
        t = T
    end

    if extensions == nothing
        extensions = NIfTI1Extension[]
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

    NIVolume(NIfTI1Header(SIZEOF_HDR, string_tuple(data_type, 10), string_tuple(db_name, 18), extents, session_error,
        regular, to_dim_info(dim_info), nidim(raw), intent_p1, intent_p2,
        intent_p3, intent_code, nidatatype(t), nibitpix(t),
        slice_start, (qfac, voxel_size..., time_step, 0, 0, 0), 352,
        scl_slope, scl_inter, slice_end, slice_code,
        xyzt_units, cal_max, cal_min, slice_duration,
        toffset, glmax, glmin, string_tuple(descrip, 80), string_tuple(aux_file, 24), (method2 || method3),
        method3, quatern_b, quatern_c, quatern_d,
        qoffset_x, qoffset_y, qoffset_z, (orientation[1, :]...,),
        (orientation[2, :]...,), (orientation[3, :]...,), string_tuple(intent_name, 16), NP1_MAGIC), extensions, raw)
end

# Calculates the size of a NIfTI extension
esize(ex::NIfTI1Extension) = 8 + ceil(Int, length(ex.edata)/16)*16

# Validates the header of a volume and updates it to match the volume's contents
function niupdate(vol::NIVolume{T}) where {T}
    vol.header.sizeof_hdr = SIZEOF_HDR
    vol.header.dim = nidim(vol.raw)
    vol.header.datatype = nidatatype(T)
    vol.header.bitpix = nibitpix(T)
    vol.header.vox_offset = isempty(vol.extensions) ? Int32(352) :
        Int32(mapreduce(esize, +, vol.extensions) + SIZEOF_HDR)
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
    write(io, vol.raw)
end

# Convenience function to write a NIfTI file given a path
function niwrite(path::AbstractString, vol::NIVolume)
    if split(path,".")[end] == "gz"
        iogz = gzopen(path, "w9")
        write(iogz, vol)
        close(iogz)
    else
        io = open(path, "w")
        write(io, vol)
        close(io)
    end
end

# Read header from a NIfTI file
function read_header(io::IO)
    header, swapped = read(io, NIfTI1Header)
    if header.sizeof_hdr != SIZEOF_HDR
        error("This is not a NIfTI-1 file")
    end
    header, swapped
end

# Read extension fields following NIfTI header
function read_extensions(io::IO, header::NIfTI1Header)
    if eof(io)
        return NIfTI1Extension[]
    end

    extension = read!(io, Array{UInt8}(undef, 4))
    if extension[1] != 1
        return NIfTI1Extension[]
    end

    extensions = NIfTI1Extension[]
    should_bswap = header.dim[1] > 7
    while header.magic == NP1_MAGIC ? position(io) < header.vox_offset : !eof(io)
        esize = read(io, Int32)
        ecode = read(io, Int32)

        if should_bswap
            esize = bswap(esize)
            ecode = bswap(ecode)
        end

        push!(extensions, NIfTI1Extension(ecode, read!(io, Array{UInt8}(undef, esize-8))))
    end
    extensions
end

# Look for a gzip header in an IOStream
function isgz(io::IO)
    ret = read(io, UInt8) == 0x1F && read(io, UInt8) == 0x8B
    seek(io, 0)
    ret
end

# Read a NIfTI file. The optional mmap argument determines whether the contents are read in full
# (if false) or mmaped from the disk (if true).
"""

"""
function niread(file::AbstractString; mmap::Bool=false)
    file_io = open(file, "r")
    header_gzipped = isgz(file_io)
    header_io = header_gzipped ? gzdopen(file_io) : file_io
    header, swapped = read_header(header_io)
    extensions = read_extensions(header_io, header)
    dims = convert(Tuple{Vararg{Int}}, header.dim[2:header.dim[1]+1])

    if !haskey(NIfTI_DT_BITSTYPES, header.datatype)
        error("data type $(header.datatype) not yet supported")
    end
    dtype = NIfTI_DT_BITSTYPES[header.datatype]

    local volume
    if header.magic == NP1_MAGIC
        if mmap
            if header_gzipped
                close(header_io)
                close(file_io)
                error("cannot mmap a gzipped NIfTI file")
            else
                volume = Mmap.mmap(header_io, Array{dtype,length(dims)}, dims, Int(header.vox_offset))
            end
        else
            seek(header_io, Int(header.vox_offset))
            volume = read!(header_io, Array{dtype}(undef, dims))
            if !eof(header_io)
                warn("file size does not match length of data; some data may be ignored")
            end
            close(header_io)
            !header_gzipped || close(file_io)
        end
    elseif header.magic == NI1_MAGIC
        close(header_io)
        !header_gzipped || close(file_io)

        volume_name = replace(file, r"\.\w+(\.gz)?$" => "")*".img"
        if !isfile(volume_name)
            if isfile(volume_name*".gz")
                volume_name *= ".gz"
            else
                error("NIfTI file is dual file storage, but $volume_name does not exist")
            end
        end

        volume_io = open(volume_name, "r")
        volume_gzipped = isgz(volume_io)
        if mmap
            if volume_gzipped
                close(volume_io)
                error("cannot mmap a gzipped NIfTI file")
            else
                volume = Mmap.mmap(volume_io, Array{dtype,length(dims)}, dims)
            end
        else
            if volume_gzipped
                volume_gz_io = gzdopen(volume_io)
                volume = read!(volume_gz_io, Array{dtype}(undef, dims))
                close(volume_gz_io)
            else
                volume = read!(volume_io, Array{dtype}(undef, dims))
            end
            close(volume_io)
        end
    end
    if swapped && sizeof(eltype(volume)) > 1
        volume = mappedarray((ntoh, hton), volume)
    end

    return NIVolume(header, extensions, volume)
end

# Allow file to be indexed like an array, but with indices yielding scaled data
@inline getindex(f::NIVolume{T}, idx::Vararg{Int}) where {T} =
    getindex(f.raw, idx...,) * f.header.scl_slope + f.header.scl_inter

add1(x::Union{AbstractArray{T},T}) where {T<:Integer} = x + 1
add1(::Colon) = Colon()
@inline vox(f::NIVolume, args...,) = getindex(f, map(add1, args)...,)
size(f::NIVolume) = size(f.raw)
size(f::NIVolume, d) = size(f.raw, d)
ndims(f::NIVolume) = ndims(f.raw)
length(f::NIVolume) = length(f.raw)
lastindex(f::NIVolume) = lastindex(f.raw)

end
