# TODO: these are things that should probably be in ImageCore but I need now
# TODO: ensure index properties on these (unique indexes etc.)


const NamedAxes{Names,N} = NamedTuple{Names,Tuple{Vararg{<:AbstractVector,N}}}

const NamedAxis{Name,Ax<:AbstractVector} = NamedTuple{Name,Tuple{Ax}}

NamedAxis{Name}(x::AbstractVector) where {Name} = NamedTuple{(Name,)}((x,))

ImageCore.namedaxes(x::NamedTuple{(),Tuple{}}, i::Int) = throw(BoundsError(x, i))
function ImageCore.namedaxes(x::NamedAxes{Sym,N}, i::Int) where {Sym,N}
    @boundscheck if 1 > 1 > N
        error(BoundsError(x, i))
    end
    @inbounds NamedTuple{(sym[i],)}((x[i],))
end

"""
    firstaxis(x) -> NamedAxis
"""
firstaxis(x::NamedTuple) = NamedTuple{(first(keys(x)),)}((first(x),))

"""
    filteraxes(f, x) -> Tuple{Vararg{AbstractVector}}

The function `f` iterates over each named axis of `x` returning the
corresponding axis if `f(namedaxis)` returns `true`.
"""
filteraxes(f::F, x) where {F} = mapfilteraxes(f, first, namedaxes(x))


"""
    mapfilteraxes(f, m, x)

`f` is a function that determines whether to use a given axis from `x` and `mf`
is a function that acts on select axes to return some new value.
"""
mapfilteraxes(f::F, m::M, x::Any) where {F,M} = _mapfilteraxes(f, m, namedaxes(x))
_mapfilteraxes(f::F, m::M, x::NamedTuple) where {F,M} = _mapfilteraxes(f, m, firstaxis(x), tail(x))
@inline function _mapfilteraxes(f::F, m::M, fa::NamedAxis, x::NamedTuple{Sym}) where {F,M,Sym}
    if f(fa)
        (m(fa), _mapfilteraxes(f, m, firstaxis(x), tail(x))...)
    else
        _mapfilteraxes(f, m, x)
    end
end
@inline function _mapfilteraxes(f::F, m::M, fa::NamedAxis, x::NamedTuple{(),Tuple{}}) where {F,M}
    if f(fa)
        return (m(fa),)
    else
        return ()
    end
end


"""
    filterdims(f, x) -> NTuple{N,Int}

The function `f` iterates over each named axis of `x` returning the
corresponding dimension position.
"""
filterdims(f::F, x) where {F} = _filterdims(f, namedaxes(x), 1)
_filterdims(f::F, x::NamedAxes, i::Int) where {F} = _filterdims(f, firstaxis(x), tail(x), i)
@inline function _filterdims(f::F, fa::NamedAxis, x::NamedAxes, i::Int) where {F}
    if f(fa)
        (i, _filteraxes(f, x, i+1)...)
    else
        _filteraxes(f, x, i+1)
    end
end
@inline function _filterdims(f::F, x::NamedAxis, i::Int) where {F}
    if f(x)
        (i,)
    else
        ()
    end
end

"""
    istimeaxis(ax) -> Bool

Test whether the axis `ax` corresponds to time.
"""
istimeaxis(::Type{<:NamedTuple{(:time,)}}) = true
istimeaxis(::Type{<:NamedTuple{Name}}) where {Name} = false


"""
    isspatialaxis(::NamedAxis) -> Bool

Determines whether a given axis refers to a spatial dimension. Default is true.
"""
isspatialaxis(::T) where {T} = isspatialaxis(T)
isspatialaxis(::Type{<:NamedAxis{Name}}) where {Name} = true
isspatialaxis(::Type{<:NamedAxis{(:time,)}}) = false
isspatialaxis(::Type{<:NamedAxis{(:color,)}}) = false

ImageCore.coords_spatial(a::ArrayInfo) = filterdims(isspatialaxis, a)

ImageCore.size_spatial(a::ArrayInfo) = mapfilteraxes(isspatialaxis, i -> length(first(i)), a)

"""
    spatialoffset(img)

Provides the offset of each dimension (i.e., where each spatial axis starts).
"""
spatialoffset(x) = first.(axes(x))

"""
    spatialunits(img)

Returns the units (i.e. Unitful.unit) that each spatial axis is measured in. If not
available `nothing` is returned for each spatial axis.
"""
spatialunits(x::Any) = mapfilteraxes(isspatialaxis, i->unit(eltype(first(i))), x)

spatialunits(x::ImageMeta) = map(i->unit(eltype(i)), pixelspacing(x))


"""
    timefirst(x)

Returns the first time point of `x`.
"""
timefirst(x::Any) = first(timeaxis(x))

"""
    timeunits(img)

Returns the units (i.e. Unitful.unit) the time axis is measured in. If not available
`nothing` is returned.
"""
timeunits(x::Any) = unit(eltype(timeaxis(x)))




