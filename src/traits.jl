
const NiftiSliceCodes = Dict{Int,String}(
    0 => "Unkown",
    1 => "Sequential+Increasing",
    2 => "Sequential+Decreasing",
    3 => "Alternating+Increasing",
    4 => "Alternating+Decreasing",
    5 => "Alternating+Increasing#2",
    6 => "Alternating+Decreasing#2")

const NiftiSliceCodesReverse = Dict{String,Int16}()
for (k, v) in NiftiSliceCodes
    NiftiSliceCodesReverse[v] = k
end

"""
    scaleslope(img) -> Float64

The values stored in each voxel can be scaled, allowing storage of voxels as
smaller datatypes (`scaleslope * stored_value + scaleintercept -> actual_value`).
These values are ignored for RGB(A) data types.
"""
scaleslope(img::NiftiFormat) = scaleslope(properties(img))
scaleslope(p::ImageProperties) = getheader(p, "scaleslope", zero(Float64))::Float64
scaleslope(A::AbstractArray) = 0.0

"""
    scaleintercept -> Float64

The values stored in each voxel can be scaled, allowing storage of voxels as
smaller datatypes (`scaleslope * stored_value + scaleintercept -> actual_value`).
These values are ignored for RGB(A) data types.
"""
scaleintercept(img::NiftiFormat) = scaleintercept(properties(img))
scaleintercept(p::ImageProperties) = getheader(p, "scaleintercept", zero(Float64))::Float64
scaleintercept(A::AbstractArray) = 0.0

# dimension info for nifti header
diminfo(img::NiftiFormat) = diminfo(properties(img))
diminfo(img::ImageInfo) = diminfo(properties(img))
diminfo(A::AbstractArray) = Int8(0)
diminfo(p::ImageProperties) = getheader(p, "diminfo", zero(Int8))

"""
    frequencydim(x) -> Int8

Which spatial dimension (1, 2, or 3) corresponds to phase acquisition. If not
applicable to scan type defaults to `0`.
"""
frequencydim(s) = diminfo(s)::Int8 & Int8(3) + 1

"""
    phasedim(x) -> Int8

Which spatial dimension (1, 2, or 3) corresponds to phase acquisition. If not
applicable to scan type defaults to `0`.
"""
phasedim(s) = (diminfo(s)::Int8 >> 2) & Int8(3) + 1

"""
    slicedim(x) -> Int8

Which dimension slices where acquired at throughout MRI acquisition. Defaults
to size of last dimension.
"""
slicedim(s) = (diminfo(s)::Int8 >> 4) + 1

"""
    slicecode(x) -> String

Indicates the timing and pattern of slice acquisition. The following codes are
defined:

* "Unkown",
* "Sequential+Increasing"
* "Sequential+Decreasing"
* "Alternating+Increasing"
* "Alternating+Decreasing"
* "Alternating+Increasing#2"
* "Alternating+Decreasing#2"

"""
slicecode(img::NiftiFormat) = slicecode(properties(img))
slicecode(p::ImageProperties) = getheader(p, "slicecode", "Unkown")
slicecode(A::AbstractArray) = "Unkown"

"""
    sliceduration(x) -> Float64

The amount of time necessary to acquire each slice throughout the MRI
acquisition. Defaults to `0.0`.
"""
sliceduration(img::NiftiFormat) = sliceduration(properties(img))
sliceduration(p::ImageProperties) = getheader(p, "sliceduration", 0.0)
sliceduration(A::AbstractArray) = 0.0

"""
    slicestart(x) -> Int

Which slice corresponds to the first slice acquired during MRI acquisition
(i.e. not padded slices). Defaults to `1`.
"""
slicestart(img::NiftiFormat) = slicestart(properties(img))
slicestart(p::ImageProperties) = getheader(p, "slicestart", 1)
slicestart(A::AbstractArray) = 1

"""
    sliceend(x)

Which slice corresponds to the last slice acquired during MRI acquisition
(i.e. not padded slices). Defaults to `size(x, slicedim(x))`.
"""
sliceend(img::NiftiFormat) = _sliceend(properties(img), size(img, slicedim(img)))
sliceend(img::ImageInfo) = _sliceend(properties(img), size(img, slicedim(img)))
_sliceend(p::ImageProperties, slice_dim_size::Int) = getheader(p, "sliceend", slice_dim_size)
sliceend(p::ImageProperties) = getheader(p, "sliceend", 1)
sliceend(A::AbstractArray) = size(A, ndims(A))



# these are for internal use in writing nifti headers

function voxoffset(s::ImageStream, v::Val{1})
    if isempty(extension(s))
        SIZEOF_HDR1 + 4
    else
        mapreduce(esize, +, extension(s)) + SIZEOF_HDR1 + 4
    end
end

function voxoffset(s::ImageStream, v::Val{2})
    if isempty(extension(s))
        SIZEOF_HDR2+4
    else
        mapreduce(esize, +, extension(s))+SIZEOF_HDR2+4
    end
end

# Gets dim to be used in header
nidim(x::ImageStream{T,N}) where {T,N} = [N, size(x)..., fill(1, 8-(N+1))...]

pixdim(s::ImageStream{T,N}, qfac::T2) where {T,N,T2} =
    [qfac, map(i->T2(ustrip(step(i.val))), axes(s))..., fill(T2(1), 8-(N+1))...]

function xyztunits(s::ImageStream)
    (get(NiftiUnitsReverse, spatunits(s)[1], Int16(0)) & 0x07) |
    (get(NiftiUnitsReverse,    timeunits(s), Int16(0)) & 0x38)
end

toffset(s::ImageStream{T,1}) where T = 0
toffset(s::ImageStream{T,2}) where T = 0
toffset(s::ImageStream{T,3}) where T = 0
toffset(s::ImageStream{T,N}) where {T,N} =
    ustrip(first(ImageFormats.axisvalues(s)[4]))
