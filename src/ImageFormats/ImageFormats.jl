module ImageFormats

using FileIO, StaticArrays, ImageMetadata, ImageAxes, MappedArrays, Mmap, GZip
using Unitful: unit, ustrip
using GeometryTypes: Point
using LinearAlgebra
using Rotations: Quat
using CoordinateTransformations: AffineMap
using Distributions: TDist,Chi,Chisq,Poisson,FDist,Beta,Binomial,Gamma,Normal,
                     NoncentralT,NoncentralChisq,Logistic,Uniform,NoncentralF,
                     GeneralizedExtremeValue,Distribution
import AxisArrays
import AxisArrays: AxisArray, Axis, axisnames, axisvalues, axisnames
import Base: read, read!, write

export ImageProperties,
       ImageInfo,
       ImageStream,
       SwapStream,
       timeunits,
       spatunits,
       spataxes,
       description,
       description!,
       auxfiles,
       auxfiles!,
       calmin,
       calmin!,
       calmax,
       calmax!,
       header,
       modality,
       data_offset,
       data_offset!,
       getheader,
       axesoffsets

include("swapstreams.jl")
include("imageproperties.jl")
include("imageinfo.jl")
include("imagestream.jl")
include("traits.jl")
include("orientation.jl")
# TODO: I'm shelfing this until we get NamedDims worked out
#include("transforms.jl")


const ImageFormat{T,N,A,Ax,F} = ImageMeta{T,N,AxisArray{T,N,A,Ax},ImageProperties{F}}

end
