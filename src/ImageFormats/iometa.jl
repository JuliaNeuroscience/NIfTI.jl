struct IOMeta{F,IOType<:IO,S}
    io::IOType
    properties::ImageProperties{F}
end

function myendian()
    if ENDIAN_BOM == 0x04030201
        return "little"
    elseif ENDIAN_BOM == 0x01020304
        return "big"
    end
end


IOMeta(io::IOType, d::AbstractDict{String,Any}; needswap::Bool=false) where {IOType<:IO} =
    IOMeta(io, ImageProperties(d); needswap=needswap)
IOMeta(io::IOType, d::ImageProperties{F}; needswap::Bool=false) where {IOType<:IO, F} =
    IOMeta{F,IOType,needswap}(io, d)
IOMeta{F}(io::IOType; needswap::Bool=false) where {IOType<:IO, F} =
    IOMeta{F,IOType,needswap}(io, ImageProperties{F}())
IOMeta(io::IOType; needswap::Bool=false) where {IOType<:IO, F} =
    IOMeta{:NOTHING,IOType,needswap}(io, ImageProperties{:NOTHING}())


Base.setindex!(d::IOMeta, X, propname::AbstractString) = setindex!(d.properties, X, propname)
Base.getindex(d::IOMeta, propname::AbstractString) = d.properties[propname]

#Base.copy(d::IOMeta) = ImageFormat{F}(deepcopy(d.io))
Base.delete!(d::IOMeta, key) = (delete!(d.properties, key); d)
Base.empty!(d::IOMeta) = (empty!(d.properties); d)
#Base.empty(d::IOMeta{F}) where F = IOMeta{F}()
Base.isempty(d::IOMeta) = isempty(d.properties)

Base.in(item, d::IOMeta) = in(item, d.properties)

Base.pop!(d::IOMeta, key, default) = pop!(d.properties, key, default)
Base.pop!(d::IOMeta, key) = pop!(d.properties, key, Base.secret_table_token)
Base.push!(d::IOMeta, kv::Pair) = insert!(d, kv[1], kv[2])
Base.push!(d::IOMeta, kv) = insert!(d.properties, kv[1], kv[2])

Base.haskey(d::IOMeta, k::String) = haskey(d.properties, k)
Base.keys(d::IOMeta) = keys(d.io)
Base.getkey(d::IOMeta, key, default) = getkey(d.properties, key, default)
Base.get!(f::Function, d::IOMeta, key::String) = get!(f, d.properties, key)

#Base.iterate(d::IOMeta) = iterate(d.properties)
#Base.iterate(d::IOMeta, state) = iterate(d.properties, state)
#Base.length(d::IOMeta) = length(d.properties)

Base.filter!(f, d::IOMeta) = filter!(f, d.properties)

# I/O Interface #
FileIO.stream(s::IOMeta{F,IOType}) where {F,IOType} = s.io::IOType

Base.seek(s::IOMeta, n::Integer) = seek(stream(s), n)
Base.position(s::IOMeta)  = position(stream(s))
Base.skip(s::IOMeta, n::Integer) = skip(stream(s), n)
Base.eof(s::IOMeta) = eof(stream(s))
Base.isreadonly(s::IOMeta) = isreadonly(stream(s))
Base.isreadable(s::IOMeta) = isreadable(stream(s))
Base.iswritable(s::IOMeta) = iswritable(stream(s))
Base.stat(s::IOMeta) = stat(stream(s))
Base.close(s::IOMeta) = close(stream(s))
Base.isopen(s::IOMeta) = isopen(stream(s))
Base.ismarked(s::IOMeta) = ismarked(stream(s))
Base.mark(s::IOMeta) = mark(stream(s))
Base.unmark(s::IOMeta) = unmark(stream(s))
Base.reset(s::IOMeta) = reset(stream(s))
Base.seekend(s::IOMeta) = seekend(stream(s))
Base.peek(s::IOMeta) = peek(stream(s))

read(s::IOMeta, n::Int) = read(stream(s), n)
read(s::IOMeta{F,IOType,true}, x::Type{<:AbstractArray}) where {F,IOType} = bswap.(read(stream(s), x))
read(s::IOMeta{F,IOType,false}, x::Type{<:AbstractArray}) where {F,IOType} = read(stream(s), x)
read(s::IOMeta, x::Type{Int8}) = read(stream(s), Int8)

read!(s::IOMeta{<:Any,<:IO,false}, a::Array) = read!(stream(s), a)
read!(s::IOMeta{<:Any,<:IO,true}, a::Array) = bswap.(read!(stream(s), a))

function read(s::IOMeta{F,IOType,false},
              T::Union{Type{Int16},Type{UInt16},Type{Int32},Type{UInt32},Type{Int64},Type{UInt64},Type{Int128},Type{UInt128},Type{Float16},Type{Float32},Type{Float64}}
             ) where {F,IOType}
    return read!(stream(s), Ref{T}(0))[]::T
end

function read(s::IOMeta{F,IOType,true},
              T::Union{Type{Int16},Type{UInt16},Type{Int32},Type{UInt32},Type{Int64},Type{UInt64},Type{Int128},Type{UInt128},Type{Float16},Type{Float32},Type{Float64}}
             ) where {F,IOType}
    return bswap(read!(stream(s), Ref{T}(0))[]::T)
end


write(s::IOMeta{F,IOType,true}, x::AbstractArray) where {F,IOType} = write(stream(s), mappedarray((ntoh, hton), x))
write(s::IOMeta{F,IOType,true}, x::T) where {F,IOType,T} = write(stream(s), bswap.(x))
write(s::IOMeta{F,IOType,false}, x::T) where {F,IOType,T} = write(stream(s), x)
