mutable struct ImageStream{T,I<:Tuple,IOType}
    io::IOType
    indices::I
    properties::ImageProperties
end

ImageStream{T}(io::IOMeta{IOType}, indices::I) where {T,I,IOType} =
    ImageStream{T,I,IOType}(stream(io), indices, properties(io))

function ImageStream(f, A::AbstractArray{T,N}; mode="w+", copyprops::Bool=false) where {T,N}
    ImageStream{T}(open(f, mode), AxisArrays.axes(A), ImageProperties(A; copyprops=copyprops))
end

ImageMetadata.properties(s::ImageStream) = s.properties
ImageMetadata.copyproperties(s::ImageStream, io::IOType) where IOType =
     ImageStream(io, deepcopy(properties(s)))
ImageMetadata.shareproperties(s::ImageStream, io::IOType) where IOType =
    ImageStream(io, properties(s))


# array like interface
Base.ndims(s::ImageStream{T,I}) where {T,I} = length(I.parameters)
Base.eltype(s::ImageStream{T,I}) where {T,I} = length(I.parameters)

Base.axes(s::ImageStream) = s.indices
Base.axes(s::ImageStream, i::Int) = s.indices[i]

Base.size(s::ImageStream) = length.(Base.axes(s))
Base.size(s::ImageStream, i::Int) = length(s.indices[i])

Base.length(s::ImageStream) = prod(size(s))

# AbstractDict Interface #
Base.setindex!(d::ImageStream, val, key::String) = setindex!(properties(d), val, key)
Base.getindex(d::ImageStream, key::String) = properties(d)[key]

#Base.copy(d::ImageStream) = ImageFormat{F}(deepcopy(properties(d)))
##Base.empty(d::ImageStream{F}) where F = ImageStream{F}()

Base.delete!(d::ImageStream, k::String) = (delete!(properties(d), k); d)
Base.empty!(d::ImageStream) = (empty!(properties(d)); d)
Base.isempty(d::ImageStream) = isempty(properties(d))

Base.in(item, d::ImageStream) = in(item, properties(d))

Base.pop!(d::ImageStream, key, default) = pop!(properties(d), key, default)
Base.pop!(d::ImageStream, key) = pop!(properties(d), key)
Base.push!(d::ImageStream, kv::Pair) = insert!(properties(d), kv)
Base.push!(d::ImageStream, kv) = push!(properties(d), kv)

Base.haskey(d::ImageStream, k::String) = haskey(properties(d), k)
Base.keys(d::ImageStream) = keys(properties(d))
Base.getkey(d::ImageStream, key, default) = getkey(properties(d), key, default)
Base.get!(f::Function, d::ImageStream, key) = get!(f, properties(d), key)

# I/O Interface #
FileIO.stream(s::ImageStream) = s.io

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
        seek(s, data_offset(s))
        Mmap.mmap(stream(s), Array{T}, size(s))
    else
        read!(stream(s), Array{T}(undef, size(s)))
    end
end

# use seek to ensure that there is not over/under shooting
function read!(s::ImageStream{T}, sink::Array{T}) where T
    seek(s, data_offset(s))
    read!(stream(s), sink)
end

@inline function read(s::ImageStream{T}, sink::Type{A}; mmap::Bool=false) where {T,A<:StaticArray}
    SA = similar_type(A, T, Size(size(s)))
    if mmap
        seek(s, data_offset(s))
        SA(Mmap.mmap(s, Array{T}, size(s)))
    else
        seek(s, data_offset(s))
        read(stream(s), SA)
    end
end

read(s::ImageStream, sink::Type{A}; kwargs...) where A<:ImageMeta =
    ImageMeta(read(s, fieldtype(A, :data); kwargs...), properties(s))

read(s::ImageStream, sink::Type{A}; kwargs...) where A<:AxisArray =
    AxisArray(read(s, fieldtype(A, :data); kwargs...), Base.axes(s))

const AxisSArray{T,N,Ax} = AxisArray{T,N,<:StaticArray,Ax}
const AxisDArray{T,N,Ax} = AxisArray{T,N,<:Array,Ax}

const ImageDAxes{T,N,Ax} = ImageMeta{T,N,<:AxisDArray{T,N,Ax}}
const ImageSAxes{T,N,Ax} = ImageMeta{T,N,<:AxisSArray{T,N,Ax}}

read(s::ImageStream) = read(s, ImageDAxes)

# this will help with batch reading or possible implementations of CLI such as fslmerge
#function readcat(f::Vector{<:AbstractString}; dims::Int) end

# this will allow chunkwise reading of images
# function readchunk() end

#= TODO: These are what I'd like supported here
ImageAxes.timeaxis(img::ImageMetaAxis) = timeaxis(data(img))
ImageAxes.timedim(img::ImageMetaAxis) = timedim(data(img))
ImageAxes.colordim(img::ImageMetaAxis) = colordim(data(img))
ImageCore.pixelspacing(img::ImageMeta) = pixelspacing(data(img))

istimeaxis,
TimeAxis,
HasTimeAxis,
sdims
coords_spatial
nimages
size_spatial
indices_spatial
assert_timedim_last
spatialorder
=#
