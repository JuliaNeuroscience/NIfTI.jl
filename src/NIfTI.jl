module NIfTI

using GZip, Mmap, ImageMetadata, ImageAxes, ImageCore, ColorTypes,
      Unitful, FileIO, Distributions, LinearAlgebra, StaticArrays,
      Distributions, MappedArrays

import Base: read, write

# using GeometryTypes
using TranscodingStreams, CodecZlib, Mmap, ImageMetadata, ImageAxes, ImageCore, ColorTypes,
      Unitful, FileIO, Distributions, LinearAlgebra, StaticArrays, Distributions, MappedArrays

import GeometryBasics: Triangle, Point, Polygon, Mesh
import Rotations: Quat
import AxisArrays: axisnames, permutation, AxisArray, axistype

import Base64

primitive type Float128 <: AbstractFloat 128 end
const ComplexF128 = Complex{Float128}

BitTypes = Union{Integer,AbstractFloat,ComplexF128,ComplexF32}

include("./ArrayFormats/ArrayFormats.jl")
using .ArrayFormats

include("./ImageFormats/ImageFormats.jl")
#using .ImageFormats

include("./NeuroIntents/NeuroIntents.jl")

#const NiftiStream{T,N,Ax,IOType} = ImageStream{T,N,Ax,ImageProperties{DataFormat{:NII}},IOType}
#const NiftiInfo{T,N,Ax} = ImageInfo{T,N,Ax,ImageProperties{DataFormat{:NII}}}
#const NiftiImage{T,N,Ax,D} = ImageMeta{T,N,AxisArray{T,N,D,Ax},ImageProperties{format"NII"}}

#NiftiFormat{T,N,Ax} = Union{NiftiStream{T,N,Ax},NiftiInfo{T,N,Ax},NiftiImage{T,N,Ax}}

include("dictionaries.jl")
include("traits.jl")
include("extension.jl")
include("load.jl")
include("readhdr1.jl")
include("readhdr2.jl")
include("write.jl")
#include("fileio.jl")

end

# TODO
# - should probably implement IndirectArrays for label images
# - check shape for Matrix intents
# - Distribution documentation
# - write intent methods
# - scl_slope, scl_inter
