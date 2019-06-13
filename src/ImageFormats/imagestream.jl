mutable struct ImageStream{T,N,Ax,P,IOType}
    io::IOType
    info::ImageInfo{T,N,Ax,P}

    function ImageStream(io::IOType, info::ImageInfo{T,N,Ax,P}) where {T,N,Ax,P,IOType}
        new{T,N,Ax,P,IOType}(io, info)
    end
end

function ImageStream{T}(io::IOType, indices::Ax, properties::AbstractDict{String,Any}; copyprops::Bool=false) where {T,Ax,IOType<:IO}
    ImageStream(io, ImageInfo{T}(indices, ImageProperties(properties; copyprops=copyprops)))
end

function ImageStream(io::IOType, A::AbstractArray{T,N}; copyprops::Bool=false) where {T,N,IOType}
    ImageStream(io, ImageInfo{T}(AxisArrays.axes(A), ImageProperties(A; copyprops=copyprops)))
end

function ImageStream(f::AbstractString, A::AbstractArray{T,N}; mode="w+", copyprops::Bool=false) where {T,N}
    open(f, mode) do io
        ImageStream(io, ImageInfo{T}(AxisArrays.axes(A), ImageProperties(A; copyprops=copyprops)))
    end
end
getinfo(img::ImageStream) = img.info

# array like interface
Base.ndims(s::ImageStream{T,N}) where {T,N} = N
Base.eltype(s::ImageStream{T,N}) where {T,N} = T

#AxisArrays.axisnames(s::ImageStream{T,N,I}) where {T,N,I} = axisnames(I)
#AxisArrays.axisnames(s::ImageStream{T,N,I}, i::Int) where {T,N,I} = axisnames(I)[i]
#AxisArrays.axisvalues(s::ImageStream) = axisvalues(axes(s)...)

getheader(s::ImageStream, k::String, default) = getheader(properties(s), k, default)

inherit_imageinfo(ImageStream)

#=
TimeAxis,
HasTimeAxis,

=#

#Base.copy(d::ImageStream) = ImageFormat{F}(deepcopy(properties(d)))
##Base.empty(d::ImageStream{F}) where F = ImageStream{F}()

Base.empty!(d::ImageStream) = (empty!(properties(d)); d)
Base.isempty(d::ImageStream) = isempty(properties(d))

Base.in(item, d::ImageStream) = in(item, properties(d))

Base.pop!(d::ImageStream, key, default) = pop!(properties(d), key, default)
Base.pop!(d::ImageStream, key) = pop!(properties(d), key)
Base.push!(d::ImageStream, kv::Pair) = insert!(properties(d), kv)
Base.push!(d::ImageStream, kv) = push!(properties(d), kv)

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




function read(s::ImageStream{T,N}, sink::Type{A};
              mmap::Bool=false, grow::Bool=true, shared::Bool=true) where {T,N,A<:Array}
    if mmap
        seek(s, data_offset(s))
        Mmap.mmap(stream(s), Array{T,N}, size(s), grow=grow, shared=shared)
    else
        read!(stream(s), Array{T,N}(undef, size(s)))
    end
end

# use seek to ensure that there is not over/under shooting
function read!(s::ImageStream{T}, sink::Array{T}) where T
    seek(s, data_offset(s))
    read!(stream(s), sink)
end

function read(s::ImageStream{T,N}, sink::Type{A}; mmap::Bool=false) where {T,N,A<:StaticArray}
    SA = similar_type(A, T, Size(size(s)))
    if mmap
        seek(s, data_offset(s))
        SA(Mmap.mmap(s, Array{T,N}, size(s)))
    else
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

read(s::ImageStream; kwarg...) = read(s, ImageDAxes; kwargs...)

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
