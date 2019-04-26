mutable struct ImageStream{T,I<:Tuple,IOType}
    io::IOType
    indices::I
end

ImageStream{T}(io::IOType, indices::I) where {T,IOType,I} = ImageStream{T,I,IOType}(io, indices)

# partial dict like interface
ImageStreamMeta{T,I,F,IOType,S} = ImageStream{T,I,IOMeta{F,IOType,S}}

Base.setindex!(d::ImageStreamMeta, val, key::String) = setindex!(d.io, val, key)
Base.getindex(d::ImageStreamMeta, key::String) = d.io[key]

#Base.copy(d::ImageStreamMeta) = ImageFormat{F}(deepcopy(d.io))
##Base.empty(d::ImageStreamMeta{F}) where F = ImageStreamMeta{F}()

Base.delete!(d::ImageStreamMeta, k::String) = (delete!(d.io, k); d)
Base.empty!(d::ImageStreamMeta) = (empty!(d.io); d)
Base.isempty(d::ImageStreamMeta) = isempty(d.io)

Base.in(item, d::ImageStreamMeta) = in(item, d.io)

Base.pop!(d::ImageStreamMeta, key, default) = pop!(d.io, key, default)
Base.pop!(d::ImageStreamMeta, key) = pop!(d.io, key)
Base.push!(d::ImageStreamMeta, kv::Pair) = insert!(d.io, kv)
Base.push!(d::ImageStreamMeta, kv) = push!(d.io, kv)

Base.haskey(d::ImageStreamMeta, k::String) = haskey(d.io, k)
Base.keys(d::ImageStreamMeta) = keys(d.io)
Base.getkey(d::ImageStreamMeta, key, default) = getkey(d.io, key, default)
Base.get!(f::Function, d::ImageStreamMeta, key) = get!(f, d.io, key)


# array like interface
Base.ndims(s::ImageStream{T,I}) where {T,I} = length(I.parameters)
Base.eltype(s::ImageStream{T,I}) where {T,I} = length(I.parameters)

Base.axes(s::ImageStream) = s.indices
Base.axes(s::ImageStream, i::Int) = s.indices[i]

Base.size(s::ImageStream) = length.(Base.axes(s))
Base.size(s::ImageStream, i::Int) = length(s.indices[i])

Base.length(s::ImageStream) = prod(size(s))

# I/O Interface #
FileIO.stream(s::ImageStream{T,I,IOType}) where {T,I,IOType} = s.io::IOType

Base.seek(s::ImageStream, n::Integer) = seek(stream(s), n)
Base.position(s::ImageStream)  = position(stream(s))
Base.skip(s::ImageStream, n::Integer) = skip(stream(s), n)
Base.eof(s::ImageStream) = eof(stream(s))
Base.isreadonly(s::ImageStream) = isreadonly(stream(s))
Base.isreadable(s::ImageStream) = isreadable(stream(s))
Base.iswritable(s::ImageStream) = iswritable(stream(s))
Base.stat(s::ImageStream) = stat(stream(s))
Base.close(s::ImageStream) = close(stream(s))
Base.isopen(s::ImageStream) = isopen(stream(s))
Base.ismarked(s::ImageStream) = ismarked(stream(s))
Base.mark(s::ImageStream) = mark(stream(s))
Base.unmark(s::ImageStream) = unmark(stream(s))
Base.reset(s::ImageStream) = reset(stream(s))
Base.seekend(s::ImageStream) = seekend(stream(s))
Base.peek(s::ImageStream) = peek(stream(s))

function read(s::ImageStream{T}, sink::Type{A}; mmap::Bool=false) where {T,A<:Array}
    if mmap
        Mmap.mmap(stream(s), Array{T}, size(s))
    else
        read!(stream(s), Array{T}(undef, size(s)))
    end
end

read!(s::ImageStream{T}, sink::Array{T}) where T = read!(stream(s), sink)

@inline function read(s::ImageStream{T}, sink::Type{A}; mmap::Bool=false) where {T,A<:StaticArray}
    SA = similar_type(A, T, Size(size(s)))
    if mmap
        SA(Mmap.mmap(s, Array{T}, size(s)))
    else
        read(stream(s), SA)
    end
end

read(s::ImageStream, sink::Type{A}; kwargs...) where A<:ImageMeta =
    ImageMeta(read(s, fieldtype(A, :data); kwargs...), s.io.properties)

read(s::ImageStream, sink::Type{A}; kwargs...) where A<:AxisArray =
    AxisArray(read(s, fieldtype(A, :data); kwargs...), Base.axes(s))

const AxisSArray{T,N,Ax} = AxisArray{T,N,<:StaticArray,Ax}
const AxisDArray{T,N,Ax} = AxisArray{T,N,<:Array,Ax}

const ImageDAxes{T,N,Ax} = ImageMeta{T,N,<:AxisDArray{T,N,Ax}}
const ImageSAxes{T,N,Ax} = ImageMeta{T,N,<:AxisSArray{T,N,Ax}}

read(s::ImageStream) = read(s, ImageDAxes)

function readcat(f::Vector{<:AbstractString}; dims::Int)
end


