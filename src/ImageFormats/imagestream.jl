mutable struct ImageStream{T,N,I<:Tuple,IOType,F}
    io::IOType
    indices::I
    properties::ImageProperties{F}
end

ImageStream{T}(io::IOMeta{IOType}, indices::Ax) where {T,Ax,IOType} =
    ImageStream{T,length(indices),Ax,IOType}(stream(io), indices, properties(io))

ImageStream{T}(io::IOType, indices::Ax, properties::AbstractDict{String,Any}; copyprops::Bool=false) where {T,Ax,IOType<:IO} =
    ImageStream{T,length(indices),Ax,IOType}(io, indices, ImageProperties(properties; copyprops=copyprops))

ImageStream(io::IOType, A::AbstractArray{T,N}; copyprops::Bool=false) where {T,N,IOType} =
    ImageStream{T,N,typeof(AxisArrays.axes(A)),IOType}(io, AxisArrays.axes(A), ImageProperties(A; copyprops=copyprops))

function ImageStream(f::AbstractString, A::AbstractArray{T,N}; mode="w+", copyprops::Bool=false) where {T,N}
    open(f, mode) do io
        ImageStream{T,N,typeof(AxisArrays.axes(A)),typeof(io)}(io, AxisArrays.axes(A), ImageProperties(A; copyprops=copyprops))
    end
end
ImageStream{T,N,Ax,IOType}(io::IOType, indices::Ax, props::ImageProperties{F}) where {T,N,Ax,IOType,F} =
    ImageStream{T,N,Ax,IOType,F}(io::IOType, indices::Ax, props::ImageProperties{F})

ImageMetadata.properties(s::ImageStream) = s.properties
ImageMetadata.copyproperties(s::ImageStream, io::IOType) where IOType =
     ImageStream(io, deepcopy(properties(s)))
ImageMetadata.shareproperties(s::ImageStream, io::IOType) where IOType =
    ImageStream(io, properties(s))
ImageCore.spacedirections(s::ImageStream) = @get s "spacedirections" ImageCore._spacedirections(s)


# array like interface
Base.ndims(s::ImageStream{T,N}) where {T,N} = N
Base.eltype(s::ImageStream{T,N}) where {T,N} = T

Base.axes(s::ImageStream) = s.indices
Base.axes(s::ImageStream, i::Int) = s.indices[i]

Base.size(s::ImageStream) = length.(axes(s))
Base.size(s::ImageStream, i::Int) = length(s.indices[i])

Base.length(s::ImageStream) = prod(size(s))

AxisArrays.axisnames(s::ImageStream{T,N,I}) where {T,N,I} = axisnames(I)
AxisArrays.axisnames(s::ImageStream{T,N,I}, i::Int) where {T,N,I} = axisnames(I)[i]
AxisArrays.axisvalues(s::ImageStream) = axisvalues(axes(s)...)

AxisArrays.axistype(s::ImageStream, i::Int) = eltype(axes(s, i))

ImageCore.spatialorder(s::ImageStream) = ImageAxes.filter_space_axes(axes(s), axisnames(s))
ImageCore.size_spatial(s::ImageStream)    = ImageAxes.filter_space_axes(axes(s), size(s))
ImageCore.indices_spatial(s::ImageStream) = ImageAxes.filter_space_axes(axes(s), axes(s))

ImageCore.coords_spatial(s::ImageStream{T,N}) where {T,N} =
    ImageAxes.filter_space_axes(axes(s), ntuple(identity, Val(N)))
ImageCore.sdims(s::ImageStream) = length(coords_spatial(s))
ImageCore.pixelspacing(s::ImageStream) = map(step, ImageAxes.filter_space_axes(axes(s), axisvalues(s)))


getheader(s::ImageStream, k::String, default) = getheader(properties(s), k, default)

# TODO: should probably not export underscores in future
ImageAxes.timeaxis(s::ImageStream) = ImageAxes._timeaxis(axes(s)...)
ImageAxes.nimages(s::ImageStream) = ImageAxes._nimages(timeaxis(s))
function ImageAxes.assert_timedim_last(s::ImageStream)
    istimeaxis(axes(s)[end]) || error("time dimension is not last")
end

ImageAxes.timedim(s::ImageStream{T,N}) where {T,N} =
    ImageAxes._timedim(ImageAxes.filter_time_axis(axes(s), ntuple(identity, Val(N))))

function ImageAxes.colordim(s::ImageStream)
    d = ImageAxes._colordim(1, axes(s))
    d > ndims(s) ? 0 : d
end

axesoffsets(img::Union{AbstractArray,ImageStream}) = map(i -> _firstindex(i) - 1, axes(img))
axesoffsets(img::Union{AbstractArray,ImageStream}, i::Int) = _firstindex(axes(img, i)) - 1

_firstindex(x::Axis) = firstindex(x.val)
_firstindex(x::AbstractArray) = firstindex(x)


#=
TimeAxis,
HasTimeAxis,

=#

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

function read(s::ImageStream{T,N}, sink::Type{A}; mmap::Bool=false) where {T,N,A<:Array}
    if mmap
        seek(s, data_offset(s))
        Mmap.mmap(stream(s), Array{T,N}, size(s))
    else
        read!(stream(s), Array{T,N}(undef, size(s)))
    end
end

# use seek to ensure that there is not over/under shooting
function read!(s::ImageStream{T}, sink::Array{T}) where T
    seek(s, data_offset(s))
    read!(stream(s), sink)
end

@inline function read(s::ImageStream{T,N}, sink::Type{A}; mmap::Bool=false) where {T,N,A<:StaticArray}
    SA = similar_type(A, T, Size(size(s)))
    if mmap
        seek(s, data_offset(s))
        SA(Mmap.mmap(s, Array{T,N}, size(s)))
    else
        seek(s, data_offset(s))
        read(stream(s), SA)
    end
end

read(s::ImageStream, sink::Type{A}; kwargs...) where A<:ImageMeta =
    ImageMeta(read(s, fieldtype(A, :data); kwargs...), properties(s))

read(s::ImageStream, sink::Type{A}; kwargs...) where A<:AxisArray =
    AxisArray(read(s, fieldtype(A, :data); kwargs...), axes(s))

const AxisSArray{T,N,Ax} = AxisArray{T,N,<:StaticArray,Ax}
const AxisDArray{T,N,Ax} = AxisArray{T,N,<:Array,Ax}

const ImageDAxes{T,N,Ax} = ImageMeta{T,N,<:AxisDArray{T,N,Ax}}
const ImageSAxes{T,N,Ax} = ImageMeta{T,N,<:AxisSArray{T,N,Ax}}

read(s::ImageStream) = read(s, ImageDAxes)

write(s::ImageStream, img::ImageMeta) = write(s, img.data)
write(s::ImageStream, a::AxisArray) = write(s, a.data)
write(s::ImageStream, x::AbstractArray) = write(stream(s), x)

write(s::ImageStream, x::Real) where T = write(stream(s), x)
write(s::ImageStream, x::String) = write(stream(s), x)

# this is a separate function to handle writing the actual image data versus metadata
#function imagewrite() end


# this will help with batch reading or possible implementations of CLI such as fslmerge
#function readcat(f::Vector{<:AbstractString}; dims::Int) end

# this will allow chunkwise reading of images
# function readchunk() end
