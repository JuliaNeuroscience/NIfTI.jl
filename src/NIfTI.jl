module NIfTI

# using GeometryTypes
using GZip, Mmap, ImageMetadata, ImageAxes, Unitful, FileIO, Distributions, LinearAlgebra
import Base.getindex, Base.size, Base.ndims, Base.length, Base.write, Base64

include("dictionaries.jl")
include("orientation.jl")
include("header.jl")
include("extension.jl")
include("imagevector.jl")
include("imagestats.jl")
include("imagelabel.jl")
include("gifti.jl")
include("imageintent.jl")
include("image.jl")
include("io.jl")

export NiftiHeader,
       sdims,
       pixelspacing,
       nimages,
       timedim
       getaffine

end
