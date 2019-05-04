mutable struct ImageStream{T,I<:Tuple,IOType}
    io::IOType
    indices::I
    properties::ImageProperties
end

ImageStream{T}(io::IOMeta{IOType}, indices::I) where {T,I,IOType} =
    ImageStream{T,I,IOType}(stream(io), indices, properties(io))

ImageStream{T}(io::IOType, indices::I, properties::AbstractDict{String,Any}; copyprops::Bool=false) where {T,I,IOType<:IO} =
    ImageStream{T,I,IOType}(io, indices, ImageProperties(properties; copyprops=copyprops))

ImageStream(io::IOType, A::AbstractArray{T}; copyprops::Bool=false) where {T,I,IOType} =
    ImageStream{T,typeof(AxisArrays.axes(A)),IOType}(io, AxisArrays.axes(A), ImageProperties(A; copyprops=copyprops))

function ImageStream(f::AbstractString, A::AbstractArray{T,N}; mode="w+", copyprops::Bool=false) where {T,N}
    open(f, mode) do io
        ImageStream{T,typeof(AxisArrays.axes(A)),typeof(io)}(io, AxisArrays.axes(A), ImageProperties(A; copyprops=copyprops))
    end
end

ImageMetadata.properties(s::ImageStream) = s.properties
ImageMetadata.copyproperties(s::ImageStream, io::IOType) where IOType =
     ImageStream(io, deepcopy(properties(s)))
ImageMetadata.shareproperties(s::ImageStream, io::IOType) where IOType =
    ImageStream(io, properties(s))
ImageCore.spacedirections(s::ImageStream) = @get s "spacedirections" ImageCore._spacedirections(s)


# array like interface
Base.ndims(s::ImageStream{T,I}) where {T,I} = length(I.parameters)
Base.eltype(s::ImageStream{T,I}) where {T,I} = length(I.parameters)

Base.axes(s::ImageStream) = s.indices
Base.axes(s::ImageStream, i::Int) = s.indices[i]

Base.size(s::ImageStream) = length.(axes(s))
Base.size(s::ImageStream, i::Int) = length(s.indices[i])

Base.length(s::ImageStream) = prod(size(s))

AxisArrays.axisnames(s::ImageStream{T,I}) where {T,I} = axisnames(I)
AxisArrays.axisvalues(s::ImageStream) = axisvalues(axes(s)...)

ImageCore.spatialorder(s::ImageStream) = ImageAxes.filter_space_axes(axes(s), axisnames(s))
ImageCore.size_spatial(s::ImageStream)    = ImageAxes.filter_space_axes(axes(s), size(s))
ImageCore.indices_spatial(s::ImageStream) = ImageAxes.filter_space_axes(axes(s), axes(s))

ImageCore.coords_spatial(s::ImageStream{T,I}) where {T,I} =
    ImageAxes.filter_space_axes(axes(s), ntuple(identity, Val(length(I.parameters))))
ImageCore.sdims(s::ImageStream) = length(coords_spatial(s))
ImageCore.pixelspacing(s::ImageStream) = map(step, ImageAxes.filter_space_axes(axes(s), axisvalues(s)))

# TODO: should probably not export underscores in future
ImageAxes.timeaxis(s::ImageStream) = ImageAxes._timeaxis(axes(s)...)
ImageAxes.nimages(s::ImageStream) = ImageAxes._nimages(timeaxis(s))
function ImageAxes.assert_timedim_last(s::ImageStream)
    istimeaxis(axes(s)[end]) || error("time dimension is not last")
end

ImageAxes.timedim(s::ImageStream{T,I}) where {T,I} =
    ImageAxes._timedim(ImageAxes.filter_time_axis(axes(s), ntuple(identity, Val(length(I.parameters)))))

function ImageAxes.colordim(s::ImageStream)
    d = ImageAxes._colordim(1, axes(s))
    d > ndims(s) ? 0 : d
end

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


