"""
    hasproperties(img) -> Bool

Returns true if `img` has a `properties` method.
"""
hasproperties(x::T) where T = hasproperties(T)
hasproperties(::Type{T}) where T = false
hasproperties(::Type{<:ImageMeta}) = true
hasproperties(::Type{<:ImageStream}) = true
hasproperties(::Type{<:ImageInfo}) = true

function noprops_error(::T) where T
    error("Type $T does not appear to have a `properties` method defined.
           Ensure that `hasproperties` and `properties` are defined for this Type.")
end

"""
    auxfiles(img) -> Vector{String}

Retrieves string for auxiliary file associated with the image.
"""

function auxfiles(p::ImageProperties)::Vector{String}
    out = get(p, "auxfiles", nothing)
    if isnothing(out)
        [""]
    else
        out
    end
end
function auxfiles(x::T) where T
    if hasproperties(x)
        auxfiles(properties(x))
    else
        [""]
    end
end

"""
    data_offset(img)

Returns offset to data given an image derived from a certain filetype. This is intended for
use by those developing and manipulating image IO.
"""
data_offset(p::ImageProperties) = get(p, "data_offset", 0)::Int
function data_offset(x::T) where T
    if hasproperties(x)
        data_offset(properties(x))
    else
        0
    end
end

"data_offset!(img, n) - Change data_offset defined in an `ImageProperties` type."
function data_offset!(p::ImageProperties, n::Int)
    p["data_offset"] = n
    return data_offset(p)
end
data_offset!(img::ImageMeta{T,N,A,<:ImageProperties}, n::Int) where {T,N,A} =
    data_offset!(properties(img), n)
function data_offset!(x::T, n::Int) where T
    if hasproperties(x)
        data_offset!(properties(x), n)
    else
        noprops_error(x)
    end
end

"""
    filename(img) -> String

File names where the image was derived from. Returns `""` by default.
"""

function FileIO.filename(p::ImageProperties)::String
    out = get(p, "filename", nothing)
    if isnothing(out)
        ""
    else
        out
    end
end

FileIO.filename(img::ImageMeta{T,N,A,<:ImageProperties}) where {T,N,A} = FileIO.filename(properties(img))
FileIO.filename(A::AbstractArray) = ""

"filename!(img, file) - Change filename defined in an `ImageProperties` type."
function filename!(p::ImageProperties, file::String)
    p["filename"] = file
    return file
end
function filename!(x::T, file::String) where T
    if hasproperties(x)
        filename!(properties(x), file)
    else
        noprops_error(x)
    end
end

"""
    modality(p)

Returns image modality that corresponds to a given `ImageProperties` instance.
"""
modality(p::ImageProperties) = get(p, "modality", "")::String
function modality(x::T) where T
    if hasproperties(x)
        modality(properties(x))
    else
        ""
    end
end

"""
    header(img)

Returns header information stored in the image properties.
"""
header(p::ImageProperties) = get(p, "header", nothing)
function header(x::T) where T
    if hasproperties(x)
        header(properties(x))
    else
        nothing
    end
end

"popheader!(p, key, default)"
popheader!(p::ImageProperties, key::String, default) = pop!(header(d), key, default)
popheader!(img::ImageMeta{T,N,A,<:ImageProperties}, key::String, default) where {T,N,A} =
    popheader!(properties(img), key, default)


"pushheader!(d, kv) - `push!` specifically for the header of ImageProperties"
pushheader!(d::ImageProperties, kv::Pair{String}) = push!(header())
pushheader!(d::ImageProperties, kv) = insert!(header(d), kv[1], kv[2])
pushheader!(img::ImageMeta{T,N,A,<:ImageProperties}, kv) where {T,N,A} =
    pushheader!(properties(img), kv)
pushheader!() where {T,N,A} =
    pushheader!(properties(img), kv)

"description(p) - Retrieves description field that may say whatever you like."
function description(p::Union{ImageProperties,ImageMeta})::String
    out = get(p, "description", nothing)
    if isnothing(out)
        ""
    else
        out
    end
end
function description(x::T) where T
    if hasproperties(x)
        description(properties(x))
    else
        ""
    end
end


"description!(p, descrip) - Change description defined in an `ImageProperties` type."
function description!(p::ImageProperties, descrip::String)
    p["descrip"] = descrip
end
description!(img::ImageMeta{T,N,A,<:ImageProperties}, descrip::String) where {T,N,A} =
    description!(properties(img), descrip)

function description!(x::T, descrip::String) where T
    if hasproperties(x)
        description!(properties(x), val)
    else
        noprops_error(x)
    end
end

"calmax(p) - Specifies maximum element for display puproses"
calmax(p::ImageProperties) = get(p, "calmax", 1)
function calmax(x::AbstractArray{T}) where T
    if hasproperties(x)
        calmax(properties(x))
    else
        maximum(x)
    end
end
function calmax(x::T) where T
    if hasproperties(x)
        calmax(properties(x))
    else
        noprops_error(x)
    end
end

"calmax!(p, val) - Change calmax defined in an `ImageProperties` type."
function calmax!(p::ImageProperties, val)
    p["calmax"] = val
end
function calmax!(x::T, val) where T
    if hasproperties(x)
        calmax!(properties(x), val)
    else
        noprops_error(x)
    end
end


"calmin - Specifies minimum element for display puproses."
calmin(d::ImageProperties) = get(d, "calmin", -1)
function calmin(x::AbstractArray{T}) where T
    if hasproperties(x)
        calmin(properties(x))
    else
        minimum(a)
    end
end
function calmin(x::T) where T
    if hasproperties(x)
        calmin(properties(x))
    else
        noprops_error(x)
    end
end

"calmin!(p, val) - Change calmin defined in an `ImageProperties` type."
function calmin!(p::ImageProperties, val)
    p["calmin"] = val
    return val
end

function calmin!(x::T, val) where T
    if hasproperties(x)
        calmin!(properties(x), val)
    else
        noprops_error(x)
    end
end





function forward_properties(::Type{A}) where A
    # TODO
    hasmethod(properties, Tuple{A}) || error("$A must have the method `properties` defined
                                             to forward properties methods.")
    @eval begin
        Base.get(d::$A, k::String, v) = get(properties(d), k, v)
        Base.get!(f::Function, d::$A, key::String) = get!(f, properties(d), key)

        Base.setindex!(d::$A, val, key::String) = setindex!(properties(d), val, key)
        Base.getindex(d::$A, key::String) = getindex(properties(d), key)

        Base.haskey(d::$A, k::String) = haskey(properties(d), k)
        Base.keys(d::$A) = keys(properties(d))
        Base.getkey(d::$A, key, default) = getkey(properties(d), key, default)

        Base.delete!(d::$A, k::String) = (delete!(properties(d), k); d)


        FileIO.filename(d::$A) = filename(properties(d))
    end
end

#forward_properties(IOMeta)
forward_properties(ImageStream)
forward_properties(ImageInfo)


"""
    spataxes(img)

Returns the axis associated with each spatial dimension. This differs from `indices_spatial`
in that it returns the AxisArrays related axes instead of the axes of the parent array.
"""
spataxes(img::AbstractArray) = ImageAxes.filter_space_axes(AxisArrays.axes(img), AxisArrays.axes(img))
spataxes(s::Union{ImageStream,ImageInfo}) = map(i->axes(s, i), coords_spatial(s))

"""
    spatunits(img)

Returns the units (i.e. Unitful.unit) that each spatial axis is measured in. If not
available `nothing` is returned for each spatial axis.
"""
spatunits(a::Union{AbstractArray,ImageStream,ImageInfo}) =
    map(i->unit(i[1]), spataxes(a))  # TODO: handle non unitful

"""
    timeunits(img)

Returns the units (i.e. Unitful.unit) the time axis is measured in. If not available
`nothing` is returned.
"""
function timeunits(a::Union{AbstractArray,ImageStream,ImageInfo})
    ta = timeaxis(a)
    if ta == nothing
        return nothing
    else
        return unit(ta[1])
    end
end
