# TODO: change slice values to reflect 1 based indexing
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
    scaleslope --> Float64

The values stored in each voxel can be scaled, allowing storage of voxels as
smaller datatypes (`scaleslope * stored_value + scaleintercept -> actual_value`).
These values are ignored for RGB(A) data types.
"""
scaleslope(img::ImageMeta{T,N,A,ImageProperties{:NII}}) where {T,N,A} = scaleslope(properties(img))
scaleslope(s::ImageStream) = scaleslope(properties(s))
scaleslope(p::ImageProperties) = getheader(p, "scaleslope", zero(Float64))::Float64
scaleslope(A::AbstractArray) = 0.0

"""
    scaleintercept -> Float64

The values stored in each voxel can be scaled, allowing storage of voxels as
smaller datatypes (`scaleslope * stored_value + scaleintercept -> actual_value`).
These values are ignored for RGB(A) data types.
"""
scaleintercept(img::ImageMeta{T,N,A,ImageProperties{:NII}}) where {T,N,A} = scaleintercept(properties(img))
scaleintercept(s::ImageStream) = scaleintercept(properties(s))
scaleintercept(p::ImageProperties) = getheader(p, "scaleintercept", zero(Float64))::Float64
scaleintercept(A::AbstractArray) = 0.0

# dimension info for nifti header
diminfo(img::ImageMeta{T,N,A,ImageProperties{:NII}}) where {T,N,A} =
    _diminfo(properties(img), size(img, N))
diminfo(s::ImageStream) = _diminfo(properties(s), size(s, ndims(s)))::Int8
_diminfo(p::ImageProperties, last_size::Int) =
    getheader(p, "diminfo", 0x00 | 0x00 | (Int8(last_size-1) << 4))
diminfo(A::AbstractArray) = 0x00 | 0x00 | (Int8(size(A, ndims(A))-1) << 4)
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
slicecode(img::ImageMeta{T,N,A,ImageProperties{:NII}}) where {T,N,A} = slicecode(properties(img))
slicecode(s::ImageStream) = slicecode(properties(s))
slicecode(p::ImageProperties) = getheader(p, "slicecode", "Unkown")
slicecode(A::AbstractArray) = "Unkown"

"""
    sliceduration(x) -> Float64

The amount of time necessary to acquire each slice throughout the MRI
acquisition. Defaults to `0.0`.
"""
sliceduration(img::ImageMeta{T,N,A,ImageProperties{:NII}}) where {T,N,A} = sliceduration(properties(img))
sliceduration(s::ImageStream) = sliceduration(properties(s))
sliceduration(p::ImageProperties) = getheader(p, "sliceduration", 0.0)
sliceduration(A::AbstractArray) = 0.0

"""
    slicestart(x) -> Int

Which slice corresponds to the first slice acquired during MRI acquisition
(i.e. not padded slices). Defaults to `1`.
"""
slicestart(img::ImageMeta{T,N,A,ImageProperties{:NII}}) where {T,N,A} = slicestart(properties(img))
slicestart(s::ImageStream) = slicestart(properties(s))
slicestart(p::ImageProperties) = getheader(p, "slicestart", 1)
slicestart(A::AbstractArray) = 1

"""
    sliceend(x)

Which slice corresponds to the last slice acquired during MRI acquisition
(i.e. not padded slices). Defaults to `size(x, slicedim(x))`.
"""
sliceend(img::ImageMeta{T,N,A,ImageProperties{:NII}}) where {T,N,A} = _sliceend(properties(img), size(img, slice_dim(img)))
sliceend(s::ImageStream) = _sliceend(properties(s), size(s, slicedim(s)))
_sliceend(p::ImageProperties, slice_dim_size::Int) = getheader(p, "sliceend", slice_dim_size)
sliceend(p::ImageProperties) = getheader(p, "sliceend", 1)
sliceend(A::AbstractArray) = size(A, ndims(A))

"""
    qform(img)
"""
qform(img::ImageMeta{T,N,A,ImageProperties{:NII}}) where {T,N,A} = qform(properties(img))
qform(s::ImageStream) = qform(properties(s))
qform(p::ImageProperties) = getheader(p, "qform", qform())
qform(A::AbstractArray) = qform()
function qform()
    SMatrix{4,4,Float64,16}(1.0, 0.0, 0.0, 0.0,
                            0.0, 1.0, 0.0, 0.0,
                            0.0, 0.0, 1.0, 0.0,
                            0.0, 0.0, 0.0, 1.0)
end

"""
    sform(A)
"""
# may just drop sform as property and always grab from spacedirections in future
sform(img::Union{AbstractArray,ImageStream}) = _sform(spacedirections(img))

function _sform(sd::Tuple{Ax1,Ax2,Ax3}) where {Ax1,Ax2,Ax3}
    SMatrix{4,4,eltype(Ax1),16}(sd[1]...,      0.0,
                               sd[2]...,      0.0,
                               sd[3]...,      0.0,
                               0.0, 0.0, 0.0, 1.0)
end

function _sform(sd::Tuple{Ax1,Ax2}) where {Ax1,Ax2}
    SMatrix{4,4,eltype(Ax1),16}(sd[1]..., 0.0, 0.0,
                                sd[2]..., 0.0, 0.0,
                                0.0, 0.0, 0.0, 0.0,
                                0.0, 0.0, 0.0, 1.0)
end

"""
qformcode
"""
qformcode(img::ImageMeta{T,N,A,ImageProperties{:NII}}) where {T,N,A} = qformcode(properties(img))
qformcode(s::ImageStream) = qformcode(properties(s))
qformcode(p::ImageProperties) = getheader(p, "qformcode", :Unkown)
qformcode(A::AbstractArray) = :Unkown

"""
    sformcode
"""
sformcode(img::ImageMeta{T,N,A,ImageProperties{:NII}}) where {T,N,A} = sformcode(properties(img))
sformcode(s::ImageStream) = sformcode(properties(s))
sformcode(p::ImageProperties) = getheader(p, "sformcode", :Unkown)
sformcode(A::AbstractArray) = :Unkown

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
function nidim(x::ImageStream)
    dim = ones(Int16, 8)
    dim[1] = ndims(x)
    dim[2:dim[1]+1] = [size(x)...]
    return dim
end

function pixdim(s::ImageStream)
    pd = fill(1.0, 8)
    # TODO: check on first pixdim measures
    pd[1] = -1
    pd[2:ndims(s)+1] = [map(i->ustrip(step(i.val)), axes(s))...]
    return pd
end

function xyztunits(s::ImageStream)
    (get(NiftiUnitsReverse, spatunits(s)[1], Int16(0)) & 0x07) |
    (get(NiftiUnitsReverse,    timeunits(s), Int16(0)) & 0x38)
end

function toffset(s::ImageStream{T,I}) where {T,I}
    if length(I.parameters) < 4
        return 0
    else
        ustrip(timeaxis(s)[1])
    end
end
