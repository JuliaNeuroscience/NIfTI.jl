"""
    ArrayInfo

# Examples
```jldoctest
```
"""
struct ArrayInfo{T,N,Sym,Ax,P}
    namedaxes::NamedTuple{Sym,Ax}
    properties::P
end

function ArrayInfo{T}(indices::NamedTuple{Sym,Ax}, properties::AbstractDict{String}, copyprops::Bool=false) where {T,Sym,Ax}
    if copyprops
        return ArrayInfo{T,length(indices),Sym,Ax,typeof(properties)}(indices, deepcopy(properties))
    else
        return ArrayInfo{T,length(indices),Sym,Ax,typeof(properties)}(indices, properties)
    end
end

ArrayInfo(a::AbstractArray, copyprops::Bool=false) = ArrayInfo(HasProperties(a), a, copyprops)

ArrayInfo(::HasProperties{true}, a::AbstractArray{T}, copyprops) where {T} = ArrayInfo{T}(namedaxes(a), properties(a), copyprops)

ArrayInfo(::HasProperties{false}, a::AbstractArray{T}, copyprops) where {T} = ArrayInfo{T}(namedaxes(a), Dict{String,Any}(), false)

ImageCore.namedaxes(a::ArrayInfo) = getproperty(a, :namedaxes)

Base.axes(a::ArrayInfo) = Tuple(namedaxes(a))

axestype(::ArrayInfo{T,N,Sym,Ax,P}) where {T,N,Sym,Ax,P} = Ax

ImageCore.HasProperties(::Type{T}) where T<:ArrayInfo = HasProperties{true}()

ImageMetadata.properties(a::ArrayInfo) = getproperty(a, :properties)

ImageCore.HasDimNames(::Type{T}) where T<:ArrayInfo = HasDimNames{true}()

Base.names(::ArrayInfo{T,N,Sym,Ax,P}) where {T,N,Sym,Ax,P} = Sym
Base.names(::ArrayInfo{T,N,Sym,Ax,P}, i::Int) where {T,N,Sym,Ax,P} = Sym[i]

# array like interface
Base.ndims(img::ArrayInfo{T,N}) where {T,N} = N

Base.eltype(img::ArrayInfo{T}) where T = T

Base.axes(img::ArrayInfo, i::Integer) = axes(img)[i]

Base.size(img::ArrayInfo) = length.(axes(img))

Base.size(img::ArrayInfo, i::Integer) = length(axes(img, i))

Base.length(img::ArrayInfo) = prod(size(img))

# ImageMetadata interface
ImageMetadata.copyproperties(img::ArrayInfo) = copyproperties(properties(s))

function ImageAxes.colordim(x::ArrayInfo)
    d = ImageAxes._colordim(1, axes(x))
    d > ndims(s) ? 0 : d
end

function ImageAxes.assert_timedim_last(s::ArrayInfo)
    istimeaxis(axes(s)[end]) || error("time dimension is not last")
end
