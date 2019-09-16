module ArrayFormats

using StaticArrays, MappedArrays, Mmap, FileIO, Unitful
using ImageCore, ImageMetadata, ImageAxes

import ImageAxes: timedim, timeaxis, istimeaxis
import AxisArrays
import AxisArrays: AxisArray, Axis, axisvalues, axisnames
import Base: read, read!, write, tail

export ArrayInfo,
       ArrayStream,
       SwapStream,
       # property traits
       dataoffset,
       dataoffset!,
       description,
       description!,
       srcfile,
       srcfile!,
       auxfiles,
       auxfiles!,
       # axis traits
       firstaxis,
       filteraxes,
       mapfilteraxes,
       filterdims,
       isspatialaxis,
       spatialoffset,
       spatialunits,
       timedim,
       timeunits,
       getter,
       setter!



include("property.jl")
include("traits.jl")
include("arrayinfo.jl")
include("compat.jl")
include("arraystream.jl")
include("swapstreams.jl")
include("imageaxes.jl")

end
