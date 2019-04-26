"""
auxfile(d)

Retrieves string for auxiliary file associated with the image.
"""
auxfile(d::ImageProperties) = get(d, "auxfile", "")
auxfile(A::AbstractArray) = ""

"""
    data_offset(d)
"""
data_offset(d::ImageProperties) = get(d, "data_offset", 0)
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
modality(d::ImageProperties) = get(d, "modality", "")
modality(A::AbstractArray) = ""

"""
    header(d)
"""
header(d::ImageProperties) = get(d, "header", :NOTHING)
header(A::AbstractArray) = :NOTHING

"""
description(d)

Retrieves description field that may say whatever you like.
"""
description(d::ImageProperties) = get(d, "description", "")
description(A::AbstractArray) = ""

"""
    calmax

Specifies maximum element for display puproses

"""
calmax(d::ImageProperties) = get(d, "calmax", 1)
calmax(A::AbstractArray{T}) where T = typemax(t)

"""
    calmin

Specifies minimum element for display puproses
"""
calmin(img::ImageProperties) = get(d, "calmax", -1)
calmin(A::AbstractArray{T}) where T = typemin(t)

function metafy(::Type{T}) where T
    @eval begin
        if any(fieldnames(Dict) .== :properties)
            auxfile(d::$T) = auxfile(getfield(d, :properties))
            header(d::$T) = header(getfield(d, :properties))
            description(d::$T) = description(getfield(d, :properties))
            data_offset(d::$T) = data_offset(getfield(d, :properties))
            modality(d::$T) = modality(getfield(d, :properties))
            calmax(d::$T) = calmax(getfield(d, :properties))
            calmin(d::$T) = calmin(getfield(d, :properties))
        else
            auxfile(d::$T) = d["auxfile"]
            header(d::$T) = d["header"]
            description(d::$T) = d["description"]
            data_offset(d::$T) = d["data_offset"]
            modality(d::$T) = d["modality"]
            calmax(d::$T) = d["calmax"]
            calmin(d::$T) = d["calmin"]
        end
    end
end

metafy(IOMeta)
metafy(ImageStream)

"""
    spataxes(img)

Returns the axis associated with each spatial dimension.
"""
spataxes(A::AbstractArray) = map(i->AxisArrays.axes(a, i), coords_spatial(A))

"""
    spatunits(img)

Returns the units (i.e. Unitful.unit) that each spatial axis is measured in. If
not available `nothing` is returned for each spatial axis.
"""
spatunits(A::AbstractArray) = map(i->unit(i.val[1]), spataxes(img))

"""
    timeunits(img)

Returns the units (i.e. Unitful.unit) the time axis is measured in. If not
available `nothing` is returned.
"""
function timeunits(A::AbstractArray)
    ta = timeaxis(A)
    if ta == nothing
        return nothing
    else
        return unit(ta[1])
    end
end
