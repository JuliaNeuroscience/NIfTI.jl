module NIfTI

using CodecZlib
using FileIO
using Mmap
using MappedArrays
using MetadataArrays
using TranscodingStreams

import Base.write, Base64

const SIZEOF_HDR1 = Int32(348)
const SIZEOF_HDR2 = Int32(540)

const NP1_MAGIC = (0x6e,0x2b,0x31,0x00)
const NI1_MAGIC = (0x6e,0x69,0x31,0x00)
const NP2_MAGIC = (0x6e,0x2b,0x32,0x00,0x0d,0x0a,0x1a,0x0a)
const NI2_MAGIC = (0x6e,0x69,0x32,0x00,0x0d,0x0a,0x1a,0x0a)

_bswap(x) = bswap(x)
_bswap(x::Tuple) = map(bswap, x)

struct NIfTI1Header
    data_type::NTuple{10,UInt8}  # unused
    db_name::NTuple{18,UInt8}  # unused
    extents::Int32  # unused
    session_error::Int16  # unused
    regular::Int8  # unused

    dim_info::Int8
    dim::NTuple{8,Int16}
    intent_parameters::NTuple{3, Float32}
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

    srow_x::NTuple{4, Float32}
    srow_y::NTuple{4, Float32}
    srow_z::NTuple{4, Float32}

    intent_name::NTuple{16, UInt8}

    magic::NTuple{4,UInt8}
end

struct NIfTI2Header
    magic::NTuple{8,UInt8}
    datatype::Int16
    bitpix::Int16
    dim::NTuple{8,Int64}
    intent_parameters::NTuple{3, Float64}
    pixdim::NTuple{8, Float64}
    vox_offset::Int64
    scl_slope::Float64
    scl_inter::Float64
    cal_max::Float64
    cal_min::Float64
    slice_duration::Float64
    toffset::Float64
    slice_start::Int64
    slice_end::Int64
    descrip::NTuple{80,UInt8}
    aux_file::NTuple{24,UInt8}
    qform_code::Int32
    sform_code::Int32
    quatern_b::Float64
    quatern_c::Float64
    quatern_d::Float64
    qoffset_x::Float64
    qoffset_y::Float64
    qoffset_z::Float64
    srow_x::NTuple{4, Float64}
    srow_y::NTuple{4, Float64}
    srow_z::NTuple{4, Float64}
    slice_code::Int32
    xyzt_units::Int32
    intent_code::Int32
    intent_name::NTuple{16, UInt8}
    dim_info::UInt8
    unused_str::NTuple{15, UInt8}
end

const NIfTIHeader = Union{NIfTI1Header, NIfTI2Header}

function define_read_header()
    packed_offsets1 = cumsum([sizeof(x) for x in NIfTI1Header.types])
    sz1 = pop!(packed_offsets1)
    pushfirst!(packed_offsets1, 0)

    packed_offsets2 = cumsum([sizeof(x) for x in NIfTI2Header.types])
    sz2 = pop!(packed_offsets2)
    pushfirst!(packed_offsets2, 0)

    @eval begin
        function read_header(io::IO, v::Int, swap::Bool)
            if v === 1
                bytes = read!(io, Array{UInt8}(undef, $sz1...))
                if swap
                    hdr = $(Expr(:new, NIfTI1Header, [:(_bswap(unsafe_load(convert(Ptr{$(NIfTI1Header.types[i])}, pointer(bytes)+$(packed_offsets1[i]))))) for i = 1:length(packed_offsets1)]...,))
                else
                    hdr = $(Expr(:new, NIfTI1Header, [:(unsafe_load(convert(Ptr{$(NIfTI1Header.types[i])}, pointer(bytes)+$(packed_offsets1[i])))) for i = 1:length(packed_offsets1)]...,))
                end
            else
                bytes = read!(io, Array{UInt8}(undef, $sz2...))
                if swap
                    hdr = $(Expr(:new, NIfTI2Header, [:(_bswap(unsafe_load(convert(Ptr{$(NIfTI2Header.types[i])}, pointer(bytes)+$(packed_offsets2[i]))))) for i = 1:length(packed_offsets2)]...,))
                else
                    hdr = $(Expr(:new, NIfTI2Header, [:(unsafe_load(convert(Ptr{$(NIfTI2Header.types[i])}, pointer(bytes)+$(packed_offsets2[i])))) for i = 1:length(packed_offsets2)]...,))
                end
            end
        end
    end
end
define_read_header()

include("intent.jl")
include("parsers.jl")
include("extensions.jl")
include("headers.jl")
include("coordinates.jl")

# Always in ms
"""
    time_step(header::NIfTI1Header)

Get the TR **in ms** from a `NIfTI1Header`.
"""
time_step(header::NIfTI1Header) =
    header.pixdim[5] * TIME_UNIT_MULTIPLIERS[header.xyzt_units >> 3]

# Gets the size of a type in bits
nibitpix(t::Type) = Int16(sizeof(t)*8)
nibitpix(::Type{Bool}) = Int16(1)

# Avoid method ambiguity
# function write(io::Base64.Base64EncodePipe, vol::NIVolume{UInt8,1})
#     invoke(write, (IO, NIVolume{UInt8,1}), io, vol)
# end

"""
    NIfTI.save(path::AbstractString, img::AbstractArray)

Write `img` to a NIfTI file at the specified `path`.
"""
function save(f::AbstractString, img; kwargs...)
    save(File{format"NII"}(f), img; kwargs...)
end
function save(f::File{format"NII"}, img; kwargs...)
    open(f, "w") do io
        save(io, img; kwargs...)
    end
end
function save(io::Stream, img; kwargs...)
    exts = extensions(img)
    hdr_version = 1
    if hdr_version === 1
        bytes = UInt8[]
        write(io, SIZEOF_HDR1)
        hdr = NIfTI1Header(img; kwargs...)
        for name in fieldnames(NIfTI1Header)
            append!(bytes, reinterpret(UInt8, [getfield(hdr, name)]))
        end
        write(io, bytes)
    else
        bytes = UInt8[]
        write(io, SIZEOF_HDR2)
        hdr = NIfTI2Header(img; kwargs...)
        for name in fieldnames(NIfTI1Header)
            append!(bytes, reinterpret(UInt8, [getfield(hdr, name)]))
        end
        write(io, bytes)
    end
    if isempty(exts)
        write(io, Int32(0))
    else
        for ex in exts
            sz = esize(ex)
            write(io, Int32(sz))
            write(io, Int32(ex.ecode))
            write(io, ex.edata)
            write(io, zeros(UInt8, sz - length(ex.edata)))
        end
    end
    T = to_eltype(hdr)
    if eltype(img) <: T
        if T <: Bool
            write(io, BitArray(img))
        else
            write(io, img)
        end
    else
        for v in img
            write(io, convert(T, v))
        end
    end
    nothing
end

"""
    NIfTI.load(file, mode::AbstractString="r"; mmap=false, scale=false)

Read a NIfTI file to an array. Set `mmap=true` to memory map the volume.
"""
function load(f::AbstractString, mode::AbstractString="r"; mmap::Bool=false, scale::Bool=false)
    open(f, mode) do io
        load(Stream{format"NII"}(splitext(f)[end] == ".gz" ? GzipDecompressorStream(io) : io, f); mmap=mmap, scale=scale)
    end
end
function load(f::File{format"NII"}, mode::AbstractString="r"; mmap::Bool=false, scale::Bool=false)
    open(f, mode) do io
        load(io; mmap=mmap, scale=scale)
    end
end
function load(s::Stream{format"NII"}; mmap::Bool=false, scale::Bool=false)
    io = stream(s)
    hdr_size = read(io, Int32)
    if hdr_size === SIZEOF_HDR1
        v = 1
        swapped = false
    elseif hdr_size === SIZEOF_HDR2
        v = 2
        swapped = false
    else
        hdr_size = bswap(hdr_size)
        if hdr_size === SIZEOF_HDR1
            v = 1
            swapped = true
        elseif hdr_size === SIZEOF_HDR
            v = 2
            swapped = true
        else
            error("This is not a NIfTI-1 file")
        end
    end
    hdr = read_header(io, v, swapped)
    isgzipped = io isa GzipDecompressorStream
    isgzipped && mmap && error("cannot mmap a gzipped NIfTI file")
    scale && mmap && throw(ArgumentError("Cannot scale mapped data"))
    offset = hdr.vox_offset - 352
    exts = Extension[]
    if !eof(io)
        b1 = read(io, UInt8)
        # GZIP doesn't skip so we read and throw away
        read(io, UInt8)
        read(io, UInt8)
        read(io, UInt8)
        if b1 !== UInt8(0)
            counter = 0
            while counter < (n - 1)
                esize = read(io, Int32)
                ec = read(io, Int32)
                data_i = Array{UInt8}(undef, esize - 8)
                read!(io, data_i)
                push!(exts, Extension(ec, data_i))
                counter += esize
            end
        end
    end

    T = to_eltype(hdr.datatype)

    dimensions = hdr.dim
    N = Int(getfield(dimensions, 1))
    sz = ntuple(i -> getfield(dimensions, i + 1), N)

    if (v === 1 && hdr.magic === NP1_MAGIC) || hdr.mage === NP2_MAGIC
        volio = io
    else
        volio = hdr_to_img(filename(io))
    end
    ArrayType = T === Bool ? BitArray{N} : Array{T,N}
    vol = mmap ? Mmap.mmap(io, ArrayType, sz) : read!(io, ArrayType(undef, sz))
    if swapped
        vol = mappedarray(bswap, vol)
    end

    if scale && hdr.scl_slope != 0
        slope = hdr.scl_slope
        inter = hdr.scl_inter
        if swapped && sizeof(T) > 1
            vol .= slope .* bswap.(vol) .+ inter
        else
            vol .= slope .* vol .+ inter
        end
        return MetadataArray(vol, (header=hdr, extension=exts))
    elseif swapped && sizeof(eltype(vol)) > 1
        return MetadataArray(mappedarray(ntoh, hton, vol), (header=hdr, extension=exts))
    else
        return MetadataArray(vol, (header=hdr, extension=exts))
    end
end

end
