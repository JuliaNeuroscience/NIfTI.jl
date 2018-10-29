module NIfTI

# using GeometryTypes
using GZip, Mmap, ImageMetadata, Unitful, FileIO, Distributions
import Base.getindex, Base.size, Base.ndims, Base.length, Base.write, Base64

include("dictionaries.jl")
include("orientation.jl")
include("header.jl")
include("extension.jl")
include("imagevector.jl")
include("imagestats.jl")
include("imagelabel.jl")
include("gifti.jl")
include("image.jl")
include("io.jl")

#export load, niread, niwrite, voxel_size, time_step, vox, getaffine, setaffine

end
