module ImageFormats

using FileIO, StaticArrays, ImageMetadata, ImageAxes, MappedArrays, Mmap, GZip
using ImageCore
using Unitful: unit, ustrip
using GeometryBasics: Point, Mesh
using LinearAlgebra
using Rotations: Quat
using CoordinateTransformations: AffineMap

import AxisArrays
import AxisArrays: AxisArray, Axis
import ImageCore: HasDimNames, HasProperties, namedaxes
import Base: read, read!, write

export ArrayInfo,
       ImageStream,
       SwapStream,
       CoordinateSpace,
       AbstractLabel,
       NodeIndex,
       # constants
       ScannerSpace,
       AnatomicalSpace,
       TailarachSpace,
       UnkownSpace,
       MNI152Space,
       orientationspace,
       # methods
       getter,
       setter!,
       timeunits,
       spatunits,
       spataxes,
       description,
       description!,
       auxfiles,
       auxfiles!,
       srcfile,
       srcfile!,
       calmin,
       calmin!,
       calmax,
       calmax!,
       freqdim,
       freqdim!,
       slicedim,
       slicedim!,
       phasedim,
       phasedim!,
       slicestart,
       slicestart!,
       sliceend,
       sliceend!,
       sliceduration,
       sliceduration!,
       magicbytes,
       magicbytes!,
       modality,
       dataoffset,
       dataoffset!,
       axesoffsets,
       sformcode,
       sformcode!,
       qformcode,
       qformcode!,
       sform,
       sform!,
       qform,
       qform!,
       orientation



include("semanticpositions.jl")
using .SemanticPositions

include("coordinatespace.jl")
include("utils.jl")
include("property.jl")
include("docproperties.jl")
include("swapstreams.jl")
include("imageinfo.jl")
include("imagestream.jl")
include("axes.jl")
include("orientation.jl")
include("nodeindex.jl")
include("labels.jl")
# TODO: I'm shelfing this until we get NamedDims worked out
#include("transforms.jl")


#const ImageFormat{T,N,A,Ax,F} = ImageMeta{T,N,AxisArray{T,N,A,Ax},ImageProperties{F}}

end
