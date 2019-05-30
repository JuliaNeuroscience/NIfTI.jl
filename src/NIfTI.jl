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

NiftiFormat = Union{ImageStream,ImageMeta{T,N,A,ImageProperties{format"NII"}}} where {T,N,A}

include("dictionaries.jl")
include("traits.jl")
include("orientation.jl")
include("extension.jl")
include("intent.jl")
include("read.jl")
include("write.jl")
include("fileio.jl")


export niread,
       niwrite,
       phasedim,
       frequencydim,
       slicedim,
       slicecode,
       sliceduration,
       sliceend,
       qform,
       sform,
       # ImageFormats
       ImageFormat,
       ImageProperties,
       IOMeta,
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
