module NIfTI

# using GeometryTypes
using TranscodingStreams, CodecZlib, Mmap, ImageMetadata, ImageAxes, ImageCore, ColorTypes,
      Unitful, FileIO, Distributions, LinearAlgebra, StaticArrays, Distributions, MappedArrays

using ImageMetadata: @get
using GeometryTypes: Triangle, Point
using Rotations: Quat
using AxisArrays: axisnames, permutation, AxisArray

import Base64

include("dictionaries.jl")
include("header.jl")
include("orientation.jl")
include("extension.jl")
include("traits.jl")
include("transform.jl")
include("read.jl")
include("write.jl")
#include("fileio.jl")

export NiftiSchema, sliceinfo, slicedim, spatunits, timeunits,
       frequencydim, phasedim, description, auxfile, niread

end

# TODO
# - should probably implement IndirectArrays for label images
# - SymmetricMatrix
# - check shape for Matrix intents
# - Distribution documentation
# - write intent methods
# - scl_slope, scl_inter
