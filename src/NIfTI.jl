module NIfTI

using GZip, Mmap, ImageMetadata, ImageAxes, ImageCore, ColorTypes,
      Unitful, FileIO, Distributions, LinearAlgebra, StaticArrays,
      Distributions, MappedArrays

import Base: read, write

using GeometryTypes: Triangle, Point
using Rotations: Quat
using AxisArrays: axisnames, permutation, AxisArray



# using GeometryTypes
using TranscodingStreams, CodecZlib, Mmap, ImageMetadata, ImageAxes, ImageCore, ColorTypes,
      Unitful, FileIO, Distributions, LinearAlgebra, StaticArrays, Distributions, MappedArrays

import GeometryTypes: Triangle, Point
import Rotations: Quat
import AxisArrays: axisnames, permutation, AxisArray

import Base64


include("ImageFormats/ImageFormats.jl")
using .ImageFormats

include("dictionaries.jl")
include("traits.jl")
include("orientation.jl")
include("extension.jl")
include("intent.jl")
include("traits.jl")
include("read.jl")
include("write.jl")
include("fileio.jl")

export niread,
       niwrite,
       phasedim,
       frequencydim,
       slicedim,
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
       auxfile,
       calmin,
       calmax,
       header,
       modality,
       data_offset


end

# TODO
# - should probably implement IndirectArrays for label images
# - SymmetricMatrix
# - check shape for Matrix intents
# - Distribution documentation
# - write intent methods
# - scl_slope, scl_inter
