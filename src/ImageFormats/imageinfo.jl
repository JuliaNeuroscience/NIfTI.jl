

"""
    ImageInfo

# Examples
```jldoctest
```
"""
struct ImageInfo{T,N,Ax,P}
    axes::Ax
    properties::P
end

function ImageInfo{T}(indices::Ax;
                      format::DataFormat=DataFormat{:Nothing}()
                     ) where {T,Ax}

    return ImageInfo{T,length(indices),Ax}(indices, ImageProperties{typeof(format)}())
end

function ImageInfo{T}(indices::Ax,
                      properties::AbstractDict{String,Any};
                      copyprops::Bool=false,
                      format::DataFormat=DataFormat{:Nothing}()
                     ) where {T,Ax}

    return ImageInfo{T,length(indices),Ax}(indices, ImageProperties{typeof(format)}(properties; copyprops=copyprops))
end

function ImageInfo(A::AbstractArray{T,N};
                   copyprops::Bool=false,
                   format::DataFormat=DataFormat{:Nothing}()
                  ) where {T,N}

    return ImageInfo{T,N,typeof(AxisArrays.axes(A))}(AxisArrays.axes(A), ImageProperties{typeof(format)}(A; copyprops=copyprops))
end

function ImageInfo(img::ImageInfo{T,N,Ax,ImageProperties{F1}};
                   copyprops::Bool=false,
                   format::DataFormat=F1()
                  ) where {T,N,Ax,F1}

    if format == F1
        return img
    else
        return ImageInfo(axes(img), ImageProperties{typeof(format)}(properties(img)))
    end
end

function ImageInfo{T,N,Ax}(indices::Ax, props::P) where {T,N,Ax,P}
    ImageInfo{T,N,Ax,P}(indices, props)
end

ImageMetadata.properties(img::ImageInfo{T,N,Ax,P}) where {T,N,Ax,P} = img.properties::P

# array like interface
Base.ndims(img::ImageInfo{T,N}) where {T,N} = N

Base.eltype(img::ImageInfo{T}) where T = T

Base.axes(img::ImageInfo{T,N,Ax}) where {T,N,Ax} = img.axes::Ax

Base.axes(img::ImageInfo, i::Integer) = axes(img)[i]

Base.size(img::ImageInfo) = length.(axes(img))

Base.size(img::ImageInfo, i::Integer) = length(axes(img, i))

Base.length(img::ImageInfo) = prod(size(img))

# ImageMetadata interface
ImageMetadata.copyproperties(img::ImageInfo) = copyproperties(properties(s))

function ImageAxes.colordim(x::ImageInfo)
    d = ImageAxes._colordim(1, axes(x))
    d > ndims(s) ? 0 : d
end

function ImageAxes.assert_timedim_last(s::ImageInfo)
    istimeaxis(axes(s)[end]) || error("time dimension is not last")
end

ImageAxes.timedim(s::ImageInfo{T,N}) where {T,N} =
    ImageAxes._timedim(ImageAxes.filter_time_axis(axes(s), ntuple(identity, Val(N))))

ImageAxes.timeaxis(img::ImageInfo) = ImageAxes._timeaxis(axes(img)...)
ImageAxes.nimages(img::ImageInfo) = ImageAxes._nimages(timeaxis(img))



# FIXME
ImageCore.spacedirections(s::ImageInfo) = @get s "spacedirections" ImageCore._spacedirections(s)

# AxisArrays interface

#Base.names(img::ImageInfo{T,N,L}) where {T,N,L} = L
#Base.names(img::ImageInfo{T,N,L}, i::Int) where {T,N,L} = L[i]
#AxisArrays.axisvalues(s::ImageStream) = axisvalues(axes(s)...)

# ImageCore interface
ImageCore.spatialorder(img::ImageInfo) = ImageAxes.filter_space_axes(axes(img), axisnames(img))

ImageCore.size_spatial(img::ImageInfo) = ImageAxes.filter_space_axes(axes(img), size(img))

ImageCore.indices_spatial(img::ImageInfo) = ImageAxes.filter_space_axes(axes(img), axes(img))

ImageCore.coords_spatial(img::ImageInfo{T,N}) where {T,N} =
    ImageAxes.filter_space_axes(axes(img), ntuple(identity, Val(N)))

ImageCore.sdims(img::ImageInfo) = length(coords_spatial(img))

ImageCore.pixelspacing(img::ImageInfo) =
    map(step, ImageAxes.filter_space_axes(axes(img), axisvalues(img)))

getimageinfo(x::AbstractArray) = nothing

AxisArrays.axisnames(s::ImageInfo{T,N,Ax}) where {T,N,Ax} = axisnames(Ax)
AxisArrays.axisnames(s::ImageInfo{T,N,Ax}, i::Int) where {T,N,Ax} = axisnames(Ax)[i]
AxisArrays.axistype(img::ImageInfo, i::Int) = eltype(axes(img, i))
AxisArrays.axisvalues(img::ImageInfo) = axisvalues(axes(img)...)

axesoffsets(img::Union{AbstractArray,ImageInfo}) = map(i -> _firstindex(i) - 1, axes(img))
axesoffsets(img::Union{AbstractArray,ImageInfo}, i::Int) = _firstindex(axes(img, i)) - 1

_firstindex(x::Axis) = firstindex(x.val)
_firstindex(x::AbstractArray) = firstindex(x)

function inherit_imageinfo(::Type{T}) where T
    @eval begin
        axesoffsets(x::$T) = axesoffsets(getinfo(x))
        axesoffsets(x::$T, i::Int) = axesoffsets(getinfo(x), i)

        ImageCore.spacedirections(x::$T) = spacedirections(getinfo(x))
        ImageCore.pixelspacing(x::$T) = pixelspacing(getinfo(x))
        ImageCore.sdims(x::$T) = sdims(getinfo(x))
        ImageCore.coords_spatial(x::$T) = coords_spatial(getinfo(x))
        ImageCore.spatialorder(x::$T) = spatialorder(getinfo(x))
        ImageCore.indices_spatial(x::$T) = indices_spatial(getinfo(x))

        ImageMetadata.properties(x::$T) = properties(getinfo(x))
        ImageMetadata.copyproperties(x::$T) = copyproperties(getinfo(x))


        AxisArrays.axisnames(x::$T) = axisnames(getinfo(x))
        AxisArrays.axisnames(x::$T, i::Int) = axisnames(getinfo(x), i)
        AxisArrays.axisvalues(x::$T) = axisvalues(getinfo(x))
        AxisArrays.axistype(x::$T) = axistype(getinfo(x))
        AxisArrays.axistype(x::$T, i::Int) = AxisArrays.axistype(getinfo(x), i)

        Base.names(x::$T) = names(getinfo(x))
#       Base.ndims(x::$T) = ndims(getinfo(x))
#       Base.eltype(x::$T) = eltype(getinfo(x))

        Base.axes(x::$T) = axes(getinfo(x))
        Base.axes(x::$T, i::Integer) = axes(getinfo(x), i)

        Base.size(x::$T, i::Integer) = size(getinfo(x), i)
        Base.size(x::$T) = size(getinfo(x))

        Base.length(x::$T) = length(getinfo(x))

        ImageAxes.timeaxis(x::$T) = timeaxis(getinfo(x))
        ImageAxes.nimages(x::$T) = timeaxis(getinfo(x))
        ImageAxes.colordim(x::$T) = colordim(getinfo(x))
        ImageAxes.assert_timedim_last(x::$T) = assert_timedim_last(getinfo(x))
        ImageAxes.timedim(x::$T) = timedim(getinfo(x))
    end
    return nothing
end

#=

axes_first(x::Tuple{Vararg{<:Axis,N}}) where N = map(i->first(i.val), x)

function axes_first(x)
    HasAxes(x) || error("no axes")
    axes_first(getaxes(x))
end

struct AxesStruct{Ax}
    axes::Ax
end

HasAxes(x::AxesStruct) = true
getaxes(x::AxesStruct{Ax}) where Ax = x.axes::Ax

axs = (Axis{:x}(1:10), Axis{:y}(1:3:10))

axstruct = AxesStruct(axs)

Base.empty!(d::ImageInfo) = (empty!(properties(d)); d)
Base.isempty(d::ImageInfo) = isempty(properties(d))

Base.in(item, d::ImageInfo) = in(item, properties(d))

Base.pop!(d::ImageInfo, key, default) = pop!(properties(d), key, default)
Base.pop!(d::ImageInfo, key) = pop!(properties(d), key)
Base.push!(d::ImageInfo, kv::Pair) = insert!(properties(d), kv)
Base.push!(d::ImageInfo, kv) = push!(properties(d), kv)

Base.get!(f::Function, d::ImageInfo, key) = get!(f, properties(d), key)


# TODO: should probably not export underscores in future
=#
#=
function d
    ImageMetadata.copyproperties(x::$T) = copyproperties(getinfo(x))
    ImageMetadata.properties(x::$T) = properties(getinfo(x))
end
=#

