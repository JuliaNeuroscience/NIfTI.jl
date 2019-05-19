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
end

"""
    modality(d)

Returns image modality that corresponds to a given `ImageProperties` instance.
"""
modality(d::ImageProperties) = get(d, "modality", "")::String
modality(a::AbstractArray) = ""

"""
    header(d)
"""
header(d::ImageProperties) = get(d, "header", nothing)
header(a::AbstractArray) = nothing

"""
    description(p)

Retrieves description field that may say whatever you like.
"""
description(p::ImageProperties) = get(p, "description", "")
description(a::AbstractArray) = ""

"""
    description!(p, descrip)

Change description defined in an `ImageProperties` type.
"""
function description!(p::ImageProperties, descrip::String)
    p["descrip"] = descrip
end

"""
    calmax

Specifies maximum element for display puproses
"""
calmax(d::ImageProperties) = get(d, "calmax", 1)
calmax(a::AbstractArray{T}) where T = maximum(a)

"""
    calmin

Specifies minimum element for display puproses
"""
calmin(d::ImageProperties) = get(d, "calmin", -1)
calmin(a::AbstractArray{T}) where T = minimum(a)

# this allows 
function metafy(::Type{T}) where T
    @eval begin
        Base.get(d::$T, k::String, v) = get(properties(d), k, v)
        Base.get!(f::Function, d::$T, key::String) = get!(f, properties(d), key)
        Base.setindex!(d::$T, val, key::String) = setindex!(properties(d), val, key)
        Base.getindex(d::$T, key::String) = getindex(properties(d), key)
        Base.haskey(d::$T, k::String) = haskey(properties(d), k)
        Base.keys(d::$T) = keys(properties(d))


        auxfiles(d::$T) = auxfiles(properties(d))
        header(d::$T) = header(properties(d))
        description(d::$T) = description(properties(d))
        data_offset(d::$T) = data_offset(properties(d))
        modality(d::$T) = modality(properties(d))
        calmax(d::$T) = calmax(properties(d))
        calmin(d::$T) = calmin(properties(d))
    end
end

metafy(IOMeta)
metafy(ImageStream)

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
spatunits(a::Union{AbstractArray,ImageStream}) = map(i->unit(i.val[1]), spataxes(a)) # TODO: handle non unitful

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
