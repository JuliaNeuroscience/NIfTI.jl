"""
    auxfiles(p)

Retrieves string for auxiliary file associated with the image.
"""
auxfiles(p::ImageProperties) = @get p "auxfiles" [""]
auxfiles(A::AbstractArray) = [""]

"""
    data_offset(p)
"""
data_offset(p::ImageProperties) = get(p, "data_offset", 0)::Int
data_offset(A::AbstractArray) = 0

"""
    data_offset!(p, n)

Change data_offset defined in an `ImageProperties` type.
"""
function data_offset!(p::ImageProperties, n::Int)
    p["data_offset"] = n
end

"""
    filename(p)

File names used to read in image.
"""
filename(p::ImageProperties) = @get p "filename" ""
filename(A::AbstractArray) = ""

"""
    filename!(p, file)

Change filename defined in an `ImageProperties` type.
"""
function filename!(p::ImageProperties, file::String)
    p["filename"] = file
    return file
end

"""
    modality(p)

Returns image modality that corresponds to a given `ImageProperties` instance.
"""
modality(p::ImageProperties) = get(p, "modality", "")::String
modality(a::AbstractArray) = ""

"""
    header(p)
"""
header(p::ImageProperties) = get(p, "header", nothing)
header(a::AbstractArray) = nothing

"""
    popheader!(p, key, default)
"""
popheader!(p::ImageProperties, key::String, default) = pop!(header(d), key, default)
popheader!(p::ImageProperties, key::String) = pop!(header(p), key)

"""
    pushheader(d, kv)
"""
pushheader!(d::ImageProperties, kv::Pair{String}) = push!(header())
pushheader!(d::ImageProperties, kv) = insert!(header(d), kv[1], kv[2])

"""
    description(p)

Retrieves description field that may say whatever you like.
"""
description(p::Union{ImageProperties,ImageMeta}) = get(p, "description", "")
description(a::AbstractArray) = ""

"""
    description!(p, descrip)

Change description defined in an `ImageProperties` type.
"""
function description!(p::ImageProperties, descrip::String)
    p["descrip"] = descrip
end

"""
    calmax(p)

Specifies maximum element for display puproses
"""
calmax(p::ImageProperties) = get(p, "calmax", 1)
calmax(a::AbstractArray{T}) where T = maximum(a)

"""
    calmax!(p, val)

Change calmax defined in an `ImageProperties` type.
"""
function calmax!(p::ImageProperties, val)
    p["calmax"] = val
end

"""
    calmin

Specifies minimum element for display puproses
"""
calmin(d::ImageProperties) = get(d, "calmin", -1)
calmin(a::AbstractArray{T}) where T = minimum(a)

"""
    calmin!(p, val)

Change calmin defined in an `ImageProperties` type.
"""
function calmin!(p::ImageProperties, val)
    p["calmin"] = val
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


        auxfiles(d::$A) = auxfiles(properties(d))
        header(d::$A) = header(properties(d))

        description(d::$A) = description(properties(d))
        description!(d::$A, descrip::String) = description!(properties(d), descript)


        filename(d::$A) = filename(properties(d))
        filename!(d::$A, file::String) = filename!(properties(d), file)

        data_offset(d::$A) = data_offset(properties(d))
        data_offset!(d::$A, n::Int) = data_offset!(properties(d), n)

        modality(d::$A) = modality(properties(d))

        calmax(d::$A) = calmax(properties(d))
        calmax!(d::$A, val) = calmax!(properties(d), val)

        calmin(d::$A) = calmin(properties(d))
        calmin!(d::$A, val) = calmin!(properties(d), val)

    end
end

forward_properties(IOMeta)
forward_properties(ImageStream)

"""
    spataxes(img)

Returns the axis associated with each spatial dimension.
"""
spataxes(a::AbstractArray) = map(i->AxisArrays.axes(a, i), coords_spatial(a))
spataxes(s::ImageStream) = map(i->axes(s, i), coords_spatial(s))


"""
    spatunits(img)

Returns the units (i.e. Unitful.unit) that each spatial axis is measured in. If
not available `nothing` is returned for each spatial axis.
"""
spatunits(a::Union{AbstractArray,ImageStream}) =
    map(i->unit(i[1]), spataxes(a))  # TODO: handle non unitful

"""
    timeunits(img)

Returns the units (i.e. Unitful.unit) the time axis is measured in. If not
available `nothing` is returned.
"""
function timeunits(a::Union{AbstractArray,ImageStream})
    ta = timeaxis(a)
    if ta == nothing
        return nothing
    else
        return unit(ta[1])
    end
end


