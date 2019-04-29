"""
    IOMeta

```jldoctest
julia> io = IOMeta(IOBuffer())
```
"""
struct IOMeta{IOType<:IO,P<:AbstractDict{String,Any}} <:IO
    io::IOType
    properties::P

    IOMeta(io::IO, d::AbstractDict{String,Any}) where {IOType<:IO} = new{typeof(io),typeof(d)}(io, d)
end

IOMeta(io::IOType) where {IOType<:IO, F} = IOMeta(io, Dict{String,Any}())


ImageMetadata.properties(io::IOMeta) = io.properties

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

#Base.iterate(d::IOMeta) = iterate(d.properties)
#Base.iterate(d::IOMeta, state) = iterate(d.properties, state)
#Base.length(d::IOMeta) = length(d.properties)

Base.filter!(f, d::IOMeta) = filter!(f, d.properties)

# I/O Interface #
FileIO.stream(s::IOMeta{IOType}) where IOType = s.io::IOType

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
read(s::IOMeta, x::Type{<:AbstractArray}) = read(stream(s), x)
read(s::IOMeta, x::Type{Int8}) = read(stream(s), Int8)


for T in (Bool,UInt128,UInt16,UInt32,UInt64,UInt8,BigInt,Int128,Int16,Int32,Int64,Int8,
          BigFloat,Float16,Float32,Float64)
    @eval begin
        read!(s::IOMeta, a::Array{$T}) = read!(stream(s), a)
    end
end


#=
function read!(s::IOMeta{<:IO,<:AbstractDict},
    a::Array{<:Union{Type{Int16},Type{UInt16},Type{Int32},Type{UInt32},Type{Int64},Type{UInt64},
                     Type{Int128},Type{UInt128},Type{Float16},Type{Float32},Type{Float64}}})
    read!(stream(s), a)
end
=#

function read(s::IOMeta,
    T::Union{Type{Int16},Type{UInt16},Type{Int32},Type{UInt32},Type{Int64},Type{UInt64},
             Type{Int128},Type{UInt128},Type{Float16},Type{Float32},Type{Float64}})
    read!(stream(s), Ref{T}(0))[]::T
end

read!(s::IOMeta, ref::Base.RefValue{<:NTuple{N,T}}) where {N,T} = read!(stream(s), ref)[]

write(s::IOMeta{<:IO,<:AbstractDict}, x::Array) = write(stream(s), x)
write(s::IOMeta{<:IO,<:AbstractDict}, x::T) where T = write(stream(s), x)
