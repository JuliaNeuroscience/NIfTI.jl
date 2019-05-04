"""
auxfile(d)

Retrieves string for auxiliary file associated with the image.
"""
auxfile(d::ImageProperties) = @get d "auxfile" ""
auxfile(A::AbstractArray) = ""

"""
    data_offset(d)
"""
data_offset(d::ImageProperties) = get(d, "data_offset", 0)::Int
data_offset(A::AbstractArray) = 0

"""
    filenames(d)

File names used to read in image.
"""
filenames(d::ImageProperties) = get(d, "filenames", String[])
filenames(A::AbstractArray) = String[]

"""
    modality(d)
"""
modality(d::ImageProperties) = get(d, "modality", "")::String
modality(a::AbstractArray) = ""

"""
    header(d)
"""
header(d::ImageProperties) = get(d, "header", nothing)
header(a::AbstractArray) = nothing

"""
description(d)

Retrieves description field that may say whatever you like.
"""
description(d::ImageProperties) = get(d, "description", "")
description(a::AbstractArray) = ""

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


        auxfile(d::$T) = auxfile(properties(d))
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
