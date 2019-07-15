module NIfTI

using GZip, Mmap, ImageMetadata, ImageAxes, ImageCore, ColorTypes,
      Unitful, FileIO, Distributions, LinearAlgebra, StaticArrays,
      Distributions, MappedArrays

import Base: read, write

# using GeometryTypes
using TranscodingStreams, CodecZlib, Mmap, ImageMetadata, ImageAxes, ImageCore, ColorTypes,
      Unitful, FileIO, Distributions, LinearAlgebra, StaticArrays, Distributions, MappedArrays

import GeometryTypes: Triangle, Point
import Rotations: Quat
import AxisArrays: axisnames, permutation, AxisArray, axistype

import Base64

include("ImageFormats/ImageFormats.jl")
using .ImageFormats
using .ImageFormats: @get

const NiftiStream{T,N,Ax,IOType} = ImageStream{T,N,Ax,ImageProperties{DataFormat{:NII}},IOType}
const NiftiInfo{T,N,Ax} = ImageInfo{T,N,Ax,ImageProperties{DataFormat{:NII}}}
const NiftiImage{T,N,Ax,D} = ImageMeta{T,N,AxisArray{T,N,D,Ax},ImageProperties{format"NII"}}

NiftiFormat{T,N,Ax} = Union{NiftiStream{T,N,Ax},NiftiInfo{T,N,Ax},NiftiImage{T,N,Ax}}

include("dictionaries.jl")
include("traits.jl")
include("orientation.jl")
include("extension.jl")
include("intent.jl")
include("read.jl")
include("write.jl")
#include("fileio.jl")

export ImageFormats,
       ImageProperties,
       ImageInfo,
       ImageStream,
       # properties
       timeunits,
       spatunits,
       data_offset,
       spataxes,
       description,
       auxfiles,
       calmin,
       calmax,
       header,
       modality


end

# TODO
# - should probably implement IndirectArrays for label images
# - SymmetricMatrix
# - check shape for Matrix intents
# - Distribution documentation
# - write intent methods
# - scl_slope, scl_inter
