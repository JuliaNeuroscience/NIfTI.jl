module ImageFormats

using FileIO, StaticArrays, AxisArrays, ImageMetadata, ImageAxes, MappedArrays
using Unitful: unit
using GeometryTypes: Point
using Rotations: Quat
using Distributions: TDist,Chi,Chisq,Poisson,FDist,Beta,Binomial,Gamma,Normal,
                     NoncentralT,NoncentralChisq,Logistic,Uniform,NoncentralF,
                     GeneralizedExtremeValue,Distribution

import Base: read, read!, write

export ImageFormat,
       ImageProperties,
       IOMeta,
       ImageStream,
       # properties
       timeunits,
       spatunits,
       data_offset,
       spataxes,
       description,
       auxfile,
       calmin,
       calmax,
       header,
       modality,
       data_offset


include("imageproperties.jl")
include("iometa.jl")
include("imagestream.jl")
include("traits.jl")
include("orientation.jl")
# TODO: I'm shelfing this until we get NamedDims worked out
#include("transforms.jl")


const ImageFormat{T,N,A,Ax,F} = ImageMeta{T,N,AxisArray{T,N,A,Ax},ImageProperties{F}}

end
