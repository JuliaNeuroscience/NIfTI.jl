@property("dataoffset", x -> Int, x -> 1)

@property("auxfiles", x -> Vector{String}, x -> [""])

@property("imagefile", x -> String, x -> "")

@property("description", x -> String, x -> "")

@property("freqdim", x -> Int, x -> 1)

@property("slicedim", x -> Int, ndims)

@property("phasedim", x -> Int, x -> 1)

# TODO slice funciton names
# These should probably be:
# - slicefirstindex
# - slicelastindex
# b/c we aren't technically returning the slice but the index of the slice
@property("slicestart", x -> Int, x -> 1)

@property("sliceend", x -> Int, x -> 1)

@property("sliceduration", x -> Float64, x -> one(Float64))

_caltype(x::AbstractArray{T}) where {T} = T
_caltype(x::Any) = Float64

_calmax(x::AbstractArray) = maximum(x)
_calmax(x::Any) = zero(Float64)
@property("calmax", _caltype, _calmax)

_calmin(x::AbstractArray) = minimum(x)
_calmin(x::Any) = zero(Float64)
@property("calmin", _caltype, _calmin)



const SFormCodes = Union{CoordinateSpace{:unknown},CoordinateSpace{:anatomical},CoordinateSpace{:tailarach},CoordinateSpace{:MNI152Space}}
const QFormCodes = Union{CoordinateSpace{:unknown},CoordinateSpace{:scanner}}

@property("sformcode", x -> CoordinateSpace , x -> UnkownSpace)

@property("qformcode", x -> CoordinateSpace, x -> UnkownSpace)


