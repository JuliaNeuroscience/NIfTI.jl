# TODO: these are things that should probably be in ImageCore but I need now
# TODO: ensure index properties on these (unique indexes etc.)

const NamedAxes{Names,N} = NamedTuple{Names,Tuple{Vararg{<:AbstractVector,N}}}

const NamedAxis{Name,Ax<:AbstractVector} = NamedTuple{Name,Tuple{Ax}}

NamedAxis{Name}(x::AbstractVector) where {Name} = NamedTuple{(Name,)}((x,))

#namedaxes(x::NamedTuple{(),Tuple{}}, i::Int) = throw(BoundsError(x, i))
function namedaxes(x::NamedAxes{Sym,N}, i::Int) where {Sym,N}
    @boundscheck if 1 > 1 > N
        error(BoundsError(x, i))
    end
    @inbounds NamedTuple{(sym[i],)}((x[i],))
end

"""
    firstaxis(::NamedAxes) -> NamedAxis
"""
firstaxis(x::NamedAxes{Sym,N}) where {Sym,N} = NamedTuple{(first(Sym),)}((first(x),))
firstaxis(x::NamedAxis) = x
firstaxis(x::NamedTuple{(),Tuple{}}) = ArgumentError("$x is empty and has no first axis.")

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
mapfilteraxes(f::F, m::M, x::Any) where {F,M} = _mapfilteraxes(f, mf, namedaxes(x))
_filteraxes(f::F, m::M, x::NamedAxes) where {F,M} = _filteraxis(f, m, firstaxis(x), tail(x))
@inline function _mapfilteraxis(f::F, m::M, fa::NamedAxis, x::NamedAxes) where {F,M}
    if f(fa)
        (m(fa), _mapfilteraxes(f, x)...)
    else
        _mapfilteraxes(f, x)
    end
end
@inline function _filteraxes(f::F, m::M, x::NamedAxis) where {F,M}
    if f(x)
        (m(x),)
    else
        ()
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
    isspatial(::NamedAxis) -> Bool

Determines whether a given axis refers to a spatial dimension. Default is true.
"""
isspatial(::NamedAxis{Name}) where {Name} = true

isspatial(::NamedAxis{(:time,)}) = false

isspatial(::NamedAxis{(:color,)}) = false

ImageCore.coords_spatial(a::ArrayInfo) = filterdims(isspatial, a)

ImageCore.size_spatial(a::ArrayInfo) = mapfilteraxes(isspatial, i -> length(first(i)), a)

"""
    spatialoffset(img)

Provides the offset of each dimension (i.e., where each spatial axis starts).
"""
spatialoffset(x) = first.(axes(x))

"""
    spatunits(img)

Returns the units (i.e. Unitful.unit) that each spatial axis is measured in. If not
available `nothing` is returned for each spatial axis.
"""
spatunits(x) = mapfilteraxes(isspatial, i->unit(eltype(first(i))), x)


ImageCore.sdims(img::ArrayInfo) = length(coords_spatial(img))

ImageCore.pixelspacing(a::ArrayInfo) = mapfilteraxes(isspatial, step, namedaxes(a))


#@property("spacedirection", axestype, indices_spatial)
ImageCore.spacedirections(a::ArrayInfo) = _spacedirections(a, properties(a))
function _spacedirections(a::ArrayInfo, p::AbstractDict)
    out = get(p, "spacedirections", PropertyMissing())
    if out isa PropertyMissing
        return indices_spatial(a)
    else
        return out
    end
end

## Time dimension
timedim(x::ArrayInfo) = finddim(x, :time)

timeaxis(x::ArrayInfo) = findaxis(x, :time)

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
timeunits(x::Any) = unit(timeaxis(x))

ImageAxes.nimages(img::ArrayInfo) = length(timeaxis(img))



