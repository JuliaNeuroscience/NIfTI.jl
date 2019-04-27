using MappedArrays
import Base: read, read!, write

export SwapStream, stream

const LittleEndian = 0x01020304
const BigEndian = 0x04030201

"""

```jldoctest
julia> s = SwapStream(IOBuffer(), true);

julia> write(s, [1:10...])

julia> seek(s, 0);

julia> read!(stream(s), Vector{Int}(undef, 10))
10-element Array{Int64,1}:
  72057594037927936
 144115188075855872
 216172782113783808
 288230376151711744
 360287970189639680
 432345564227567616
 504403158265495552
 576460752303423488
 648518346341351424
 720575940379279360

julia> seek(s, 0);

julia> read!(s, Vector{Int}(undef, 10))
10-element Array{Int64,1}:
  1
  2
  3
  4
  5
  6
  7
  8
  9
 10


```
"""
struct SwapStream{S,IOType<:IO} <: IO
    io::IOType

    SwapStream{S}(io::IOType) where {S,IOType<:IO} = new{S,IOType}(io)
end

SwapStream(io::IO, endianness::UInt32) =
    SwapStream{ENDIAN_BOM != endianness}(io)
SwapStream(io::IO, needswap::Bool=false) = SwapStream{needswap}(io)


FileIO.stream(s::SwapStream{S,IOType}) where {S,IOType} = s.io::IOType
Base.seek(s::SwapStream, n::Integer) = seek(stream(s), n)
Base.position(s::SwapStream)  = position(stream(s))
Base.skip(s::SwapStream, n::Integer) = skip(stream(s), n)
Base.eof(s::SwapStream) = eof(stream(s))
Base.isreadonly(s::SwapStream) = isreadonly(stream(s))
Base.isreadable(s::SwapStream) = isreadable(stream(s))
Base.iswritable(s::SwapStream) = iswritable(stream(s))
Base.stat(s::SwapStream) = stat(stream(s))
Base.close(s::SwapStream) = close(stream(s))
Base.isopen(s::SwapStream) = isopen(stream(s))
Base.ismarked(s::SwapStream) = ismarked(stream(s))
Base.mark(s::SwapStream) = mark(stream(s))
Base.unmark(s::SwapStream) = unmark(stream(s))
Base.reset(s::SwapStream) = reset(stream(s))
Base.seekend(s::SwapStream) = seekend(stream(s))

read(s::SwapStream, n::Int) = read(stream(s), n)
read(s::SwapStream{true}, x::Type{<:AbstractArray}) = bswap.(read(stream(s), x))
read(s::SwapStream{false}, x::Type{<:AbstractArray}) = read(stream(s), x)
read(s::SwapStream, x::Type{Int8}) = read(stream(s), Int8)

read!(s::SwapStream{false}, a::Array) = read!(stream(s), a)
read!(s::SwapStream{true}, a::Array) = bswap.(read!(stream(s), a))

function read(s::SwapStream{false},
              T::Union{Type{Int16},Type{UInt16},Type{Int32},Type{UInt32},Type{Int64},Type{UInt64},Type{Int128},Type{UInt128},Type{Float16},Type{Float32},Type{Float64}}
             )
    return read!(stream(s), Ref{T}(0))[]::T
end

function read(s::SwapStream{true},
              T::Union{Type{Int16},Type{UInt16},Type{Int32},Type{UInt32},Type{Int64},Type{UInt64},Type{Int128},Type{UInt128},Type{Float16},Type{Float32},Type{Float64}}
             )
    return bswap(read!(stream(s), Ref{T}(0))[]::T)
end

read!(s::SwapStream{true}, ref::Base.RefValue{<:NTuple{N,T}}) where {N,T} = bswap.(read!(stream(s), ref)[])
read!(s::SwapStream{false}, ref::Base.RefValue{<:NTuple{N,T}}) where {N,T} = read!(stream(s), ref)[]

write(s::SwapStream{true,<:IO}, x::Array) = write(stream(s), mappedarray(ntoh, hton, x))
write(s::SwapStream{true,<:IO}, x::T) where T = write(stream(s), bswap.(x))
write(s::SwapStream{false,<:IO}, x::T) where T = write(stream(s), x)
