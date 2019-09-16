struct ArrayStream{T,N,Names,Ax,P,IOType}
    io::IOType
    info::ArrayInfo{T,N,Names,Ax,P}

    function ArrayStream(io::IOType, info::ArrayInfo{T,N,Names,Ax,P}) where {T,N,Names,Ax,P,IOType<:IO}
        new{T,N,Names,Ax,P,IOType}(io, info)
    end
end

ImageCore.HasDimNames(::Type{T}) where T<:ArrayStream = HasDimNames{true}()
Base.names(a::ArrayStream{T,N,Names}) where {T,N,Names} = Names
Base.axes(a::ArrayStream) = axes(a.info)

ImageCore.HasProperties(::Type{T}) where T<:ArrayStream = HasProperties{true}()
ImageMetadata.properties(a::ArrayStream) = properties(a.info)

ImageCore.pixelspacing(a::ArrayStream) = pixelspacing(a.info)
ImageCore.spacedirections(a::ArrayStream) = spacedirections(a.info)

function ArrayStream(io::IOType, a::A, copyprops::Bool=false) where {A,IOType}
    _arraystream(HasProperties(A), io, a, copyprops)
end

function _arraystream(::HasProperties{true}, io::IOType, a::AbstractArray{T}, copyprops::Bool) where {IOType,T}
    if copyprops
        ArrayStream(io, ArrayInfo{T}(namedaxes(a), deepcopy(properties(a))))
    else
        ArrayStream(io, ArrayInfo{T}(namedaxes(a), properties(a)))
    end
end

function _arraystream(::HasProperties{false}, io::IOType, a::AbstractArray{T}, copyprops::Bool) where {IOType,T}
    ArrayStream(io, ArrayInfo{T}(axes(a), Dict{String,Any}))
end

function ArrayStream{T}(io::IOType,
                        indices::Ax,
                        properties::AbstractDict{String,Any};
                        copyprops::Bool=false,
                       ) where {T,Ax,IOType<:IO}

    return ArrayStream(io, ArrayInfo{T}(indices, ImageProperties{typeof(format)}(properties; copyprops=copyprops)))
end

getinfo(a::ArrayStream) = getfield(a, :info)

# array like interface
Base.ndims(s::ArrayStream{T,N}) where {T,N} = N
Base.eltype(s::ArrayStream{T,N}) where {T,N} = T
Base.size(s::ArrayStream) = size(s.info)
Base.size(s::ArrayStream, i) = size(s.info, i)


getheader(s::ArrayStream, k::String, default) = getheader(properties(s), k, default)

#inherit_imageinfo(ArrayStream)

#Base.copy(d::ArrayStream) = ImageFormat{F}(deepcopy(properties(d)))
##Base.empty(d::ArrayStream{F}) where F = ArrayStream{F}()

Base.empty!(d::ArrayStream) = (empty!(properties(d)); d)
Base.isempty(d::ArrayStream) = isempty(properties(d))

Base.in(item, d::ArrayStream) = in(item, properties(d))

Base.pop!(d::ArrayStream, key, default) = pop!(properties(d), key, default)

Base.pop!(d::ArrayStream, key) = pop!(properties(d), key)

Base.push!(d::ArrayStream, kv::Pair) = insert!(properties(d), kv)

Base.push!(d::ArrayStream, kv) = push!(properties(d), kv)

# I/O Interface #
FileIO.stream(s::ArrayStream) = s.io

Base.seek(s::ArrayStream, n::Integer) = seek(stream(s), n)

Base.position(s::ArrayStream)  = position(stream(s))

Base.skip(s::ArrayStream, n::Integer) = skip(stream(s), n)

Base.eof(s::ArrayStream) = eof(stream(s))

Base.isreadonly(s::ArrayStream) = isreadonly(stream(s))

Base.isreadable(s::ArrayStream) = isreadable(stream(s))

Base.iswritable(s::ArrayStream) = iswritable(stream(s))

Base.stat(s::ArrayStream) = stat(stream(s))

Base.close(s::ArrayStream) = close(stream(s))

Base.isopen(s::ArrayStream) = isopen(stream(s))

Base.ismarked(s::ArrayStream) = ismarked(stream(s))

Base.mark(s::ArrayStream) = mark(stream(s))

Base.unmark(s::ArrayStream) = unmark(stream(s))

Base.reset(s::ArrayStream) = reset(stream(s))

Base.seekend(s::ArrayStream) = seekend(stream(s))

Base.peek(s::ArrayStream) = peek(stream(s))

function read(s::ArrayStream{T,N}, sink::Type{A};
              mmap::Bool=false, grow::Bool=true, shared::Bool=true) where {T,N,A<:Array}
    if mmap
        seek(s, dataoffset(s))
        Mmap.mmap(stream(s), Array{T,N}, size(s), grow=grow, shared=shared)
    else
        read!(stream(s), Array{T,N}(undef, size(s)))
    end
end

# use seek to ensure that there is not over/under shooting
function read!(s::ArrayStream{T}, sink::Array{T}) where T
    seek(s, dataoffset(s))
    read!(stream(s), sink)
end

function read(s::ArrayStream{T,N}, sink::Type{A}; mmap::Bool=false) where {T,N,A<:StaticArray}
    SA = similar_type(A, T, Size(size(s)))
    if mmap
        seek(s, dataoffset(s))
        SA(Mmap.mmap(s, Array{T,N}, size(s)))
    else
        read(stream(s), SA)
    end
end

function read(s::ArrayStream, sink::Type{A}; kwargs...) where A<:ImageMeta
    ImageMeta(read(s, fieldtype(A, :data); kwargs...), properties(s))
end

function read(s::ArrayStream, sink::Type{A}; kwargs...) where A<:AxisArray
    AxisArray(read(s, fieldtype(A, :data); kwargs...), namedaxes2axisarray(namedaxes(s)))
end

const AxisSArray{T,N,Ax} = AxisArray{T,N,<:StaticArray,Ax}
const AxisDArray{T,N,Ax} = AxisArray{T,N,<:Array,Ax}

const ImageDAxes{T,N,Ax} = ImageMeta{T,N,<:AxisDArray{T,N,Ax}}
const ImageSAxes{T,N,Ax} = ImageMeta{T,N,<:AxisSArray{T,N,Ax}}

read(s::ArrayStream; kwargs...) = read(s, ImageDAxes; kwargs...)

write(s::ArrayStream, a::ImageMeta) = write(s, a.data)
write(s::ArrayStream, a::AxisArray) = write(s, a.data)
write(s::ArrayStream, x::AbstractArray) = write(stream(s), x)

write(s::ArrayStream, x::Real) where T = write(stream(s), x)
write(s::ArrayStream, x::String) = write(stream(s), x)

# this is a separate function to handle writing the actual image data versus metadata
#function imagewrite() end


# this will help with batch reading or possible implementations of CLI such as fslmerge
#function readcat(f::Vector{<:AbstractString}; dims::Int) end

# this will allow chunkwise reading of images
# function readchunk() end
