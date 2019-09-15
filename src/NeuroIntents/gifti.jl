# DISCLAIMER!!! Everything in this file should be considered pre-alpha and subject to change
"""
    GiftiGeneric

A generic GIfTI format
"""
struct GiftiGeneric end
intent2ext(::GiftiGeneric) = ".gii"

"""
    GiftiCoordinates
"""
struct GiftiCoordinates end
intent2ext(::GiftiCoordinates)= ".coord.gii"

"""
    GiftiFunctional
"""
struct GiftiFunctional end
intent2ext(::GiftiFunctional) = ".func.gii"


"""
    Labels
"""
struct GiftiLabels end
intent2ext(::GiftiLabels) = ".label.gii"


struct GiftiRGB end
struct GiftiRGBA end
"""
    GiftiColor

Extension for GiftiRGB or GiftiRGBA
"""
const GiftiColor = Union{GiftiRGB,GiftiRGBA}
intent2ext(::GiftiRGB) = ".rgba.gii"
intentaxis(::Type{<:GiftiColor}) = :color

"""
    GiftiShape
"""
struct GiftiShape end
intent2ext(::GiftiShape) = ".shape.gii"
intentaxis(::Type{<:GiftiShape}) = :shape

"""
    GiftiSurface
"""
struct GiftiSurface end
intent2ext(::GiftiSurface) = ".surf.gii"
intentaxis(::Type{<:GiftiSurface}) = :surface

"""
Tensors
"""
struct GiftiTensors end
intent2ext(::GiftiTensors) = ".tensor.gii"

"""
    GiftiTimeSeries
"""
struct GiftiTimeSeries end
intent2ext(::GiftiTimeSeries) = ".time.gii"
intentaxis(::Type{<:GiftiTimeSeries}) = :time

"""
    GiftiTopology
"""
struct GiftiTopology end
intent2ext(::GiftiTopology) = ".topo.gii"
intentaxis(::Type{<:GiftiTopology}) = :topology

"""
    GiftiVector
"""
struct GiftiVector end
intent2ext(::GiftiVector) = ".vector.gii"
intentaxis(::Type{<:GiftiVector}) = :vector

function giftiintent(i::Integer)
    if i == 2001
        return TimeSeries
    elseif i == 2002
        return Vector{Point}
    elseif i == 2003
        return GiftiRGB
    elseif i == 2004
        return GiftiRGBA
    elseif i == 2005
        return Polygon
    elseif i == 2006
       return FSLDisplacementVector # TODO
    elseif i == 2007
        return FSLCubicSplineCoefficient  # TODO
    elseif i == 2008
        return FSLDCTCoefficients  # TODO
    elseif i == 2009
        return FSLQuadraticSplineCoefficients  # TODO
    elseif i == 2016
        return FSLTopupCubicSplineCoefficients  # TODO
    elseif i == 2017
        return TopupQuadraticSplineCoefficients # TODO
    elseif i == 2018
        return TopupField  # TODO
    else
        return UnkownIntent
    end
end

