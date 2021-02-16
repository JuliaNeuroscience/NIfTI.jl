
###
### read
###
#=
@generated read_header(io, ::Type{T}) where {T} = _read_header(T)
@generated read_header_swap(io, ::Type{T}) where {T} = _read_header_swap(T)

function _read_header(::Type{T}) where {T}
    if T <: NIfTI1Header
        out = Expr(:new, T, Int32(348))
    else
        out = Expr(:new, T, Int32(540))
    end
    for i in 2:nfields(T)
        T_i = fieldtype(T, i)
        if T_i <: Tuple
            t = Expr(:tuple)
            for p in T_i.parameters
                push!(t.args, :(read(io, $p)))
            end
            push!(out.args, t)
        else
            push!(out.args, :(read(io, $T_i)))
        end
    end
    return out
end

function _read_header_swap(::Type{T}) where {T}
    if T <: NIfTI1Header
        out = Expr(:new, T, Int32(348))
    else
        out = Expr(:new, T, Int32(540))
    end
    for i in 2:nfields(T)
        T_i = fieldtype(T, i)
        if T_i <: Tuple
            t = Expr(:tuple)
            for p in T_i.parameters
                push!(t.args, :(bswap(read(io, $p))))
            end
            push!(out.args, t)
        else
            push!(out.args, :(bswap(read(io, $T_i))))
        end
    end
    return out
end

###
### write
###
@generated write_header(io, x::X) where {X} = _generate_write(X)


function _generate_write(::Type{T}) where {T}
    if T <: NIfTI1Header
        out = Expr(:block, write(io, Int32(348)))
    else
        out = Expr(:block, write(io, Int32(540)))
    end
    for i in 2:nfields(T)
        T_i = fieldtype(T, i)
        n = fieldname(T, i)
        if T_i <: Tuple
            push!(out.args, write(io, Ref(x.$n)))
        else
            push!(out.args, write(io, x.$n))
        end
    end
    return out
end
=#

struct NIfTI1Header
    data_type::NTuple{10,UInt8}
    db_name::NTuple{18,UInt8}
    extents::Int32
    session_error::Int16
    regular::Int8

    dim_info::UInt8
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
    xyzt_units::UInt8
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

struct NIfTI2Header
    magic::NTuple{8,UInt8}
    datatype::Int16
    bitpix::Int16
    dim::NTuple{8,Int64}
    intent_p1::Float64
    intent_p2::Float64
    intent_p3::Float64
    pixdim::NTuple{8,Float64}
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
    srow_x::NTuple{4,Float64}
    srow_y::NTuple{4,Float64}
    srow_z::NTuple{4,Float64}
    slice_code::Int32
    xyzt_units::Int32
    intent_code::Int32
    intent_name::NTuple{16,UInt8}
    dim_info::UInt8
    unused_str::NTuple{15,UInt8}
end

sizeof_hdr(::NIfTI1Header) = Int32(348)

sizeof_hdr(::NIfTI2Header) = Int32(540)


function define_packed(T::DataType)
    offsets = [sizeof(x) for x in T.types]
    packed_offsets = cumsum(offsets)
    sz = pop!(packed_offsets)
    pushfirst!(packed_offsets, 0)

    @eval begin
        function read_header(io::IO, ::Type{$T})
            bytes = read!(io, Array{UInt8}(undef, ($sz,)))
            ptr = pointer(bytes)
            hdr = $(Expr(:new, T, [:(unsafe_load(convert(Ptr{$(T.types[i])}, ptr+$(packed_offsets[i])))) for i = 1:length(packed_offsets)]...,))
            return hdr
        end

        function read_header_swap(io::IO, ::Type{$T})
            bytes = read!(io, Array{UInt8}(undef, ($sz,)))
            ptr = pointer(bytes)
            hdr = $(Expr(:new, T, [:(_byteswap(unsafe_load(convert(Ptr{$(T.types[i])}, ptr+$(packed_offsets[i]))))) for i = 1:length(packed_offsets)]...,))
            return hdr
        end
        #=
        function Base.write(io::IO, x::$T)
            bytes = $(ifelse(T<:NIfTI1Header, UInt8[0x5c,0x01,0x00,0x00], UInt8[0x1c,0x02,0x00,0x00]))
            for name in fieldnames($T)
                append!(bytes, reinterpret(UInt8, [getfield(x,name)]))
            end
            write(io, bytes)
            $(sz + 4)
        end
        =#
    end
    nothing
end

function write_header(io, x::NIfTI2Header)
    write(io, Int32(540))
    _write_header(io, x)
end

function write_header(io, x::NIfTI1Header)
    write(io, Int32(348))
    _write_header(io, x)
end


@generated function _write_header(io, x::T) where {T}
    out = Expr(:block)
    for i in 1:nfields(T)
        T_i = fieldtype(T, i)
        n = fieldname(T, i)
        if T_i <: Tuple
            push!(out.args, :(write(io, Ref(x.$n))))
        else
            push!(out.args, :(write(io, x.$n)))
        end
    end
    return out
end


#=
@generated function write_header(io, x::T) where {T}
    generate_write_header(T)
end

function generate_write_header(::Type{T}) where {T}
    if T <: NIfTI1Header
        out = Expr(:block, :(write(io, Int32(348))))
    else
        out = Expr(:block, :(write(io, Int32(540))))
    end
    for i in 2:nfields(T)
        T_i = fieldtype(T, i)
        n = fieldname(T, i)
        if T_i <: Tuple
            push!(out.args, write(io, Ref(x.$n)))
        else
            push!(out.args, write(io, x.$n))
        end
    end
    return out
end
=#

_byteswap(x::UInt8) = x
_byteswap(x::Int8) = x
_byteswap(x) = bswap(x)
_byteswap(x::Tuple{Vararg{UInt8}}) = x
_byteswap(x::Tuple{Vararg{Int8}}) = x
_byteswap(x::Tuple{Vararg{Any}}) = map(bswap, x)

define_packed(NIfTI1Header)
define_packed(NIfTI2Header)


niopen(file::AbstractString, mode::AbstractString) = niopen(open(file, mode))
@inline function niopen(io)
    b1 = read(io, UInt8)
    b2 = read(io, UInt8)
    if b1 === 0x1F && b2 === 0x8B
        seek(io, 0)
        return niopen(GzipDecompressorStream(io))
    else
        b3 = read(io, UInt8)
        b4 = read(io, UInt8)
        if (b1,b2,b3,b4) === (0x5c,0x01,0x00,0x00)
            return io, read_header(io, NIfTI1Header), false
        elseif (b1,b2,b3,b4) === (0x1c,0x02,0x00,0x00)
            return io, read_header(io, NIfTI2Header), false
        elseif (b1,b2,b3,b4) === (0x00,0x00,0x01,0x5c)
            return io, read_header_swap(io, NIfTI1Header), true
        elseif (b1,b2,b3,b4) === (0x00,0x00,0x02,0x1c)
            return io, read_header_swap(io, NIfTI2Header), true
        else
            error("This is not a NIfTI-1 file")
        end
    end
end

