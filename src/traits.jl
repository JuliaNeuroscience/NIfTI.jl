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
scaleslope(img::ImageMeta{T,N,A,ImageProperties{format"NII"}}) where {T,N,A} = scaleslope(properties(img))
scaleslope(s::ImageStream) = scaleslope(properties(s))
scaleslope(p::ImageProperties) = getheader(p, "scaleslope", zero(Float64))::Float64
scaleslope(A::AbstractArray) = 0.0

"""
    scaleintercept -> Float64

The values stored in each voxel can be scaled, allowing storage of voxels as
smaller datatypes (`scaleslope * stored_value + scaleintercept -> actual_value`).
These values are ignored for RGB(A) data types.
"""
scaleintercept(img::ImageMeta{T,N,A,ImageProperties{format"NII"}}) where {T,N,A} = scaleintercept(properties(img))
scaleintercept(s::ImageStream) = scaleintercept(properties(s))
scaleintercept(p::ImageProperties) = getheader(p, "scaleintercept", zero(Float64))::Float64
scaleintercept(A::AbstractArray) = 0.0

# dimension info for nifti header
diminfo(img::ImageMeta{T,N,A,ImageProperties{format"NII"}}) where {T,N,A} = diminfo(properties(img))
diminfo(s::ImageStream) = Int8(0)
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
slicecode(img::ImageMeta{T,N,A,ImageProperties{format"NII"}}) where {T,N,A} = slicecode(properties(img))
slicecode(s::ImageStream) = slicecode(properties(s))
slicecode(p::ImageProperties) = getheader(p, "slicecode", "Unkown")
slicecode(A::AbstractArray) = "Unkown"

"""
    sliceduration(x) -> Float64

The amount of time necessary to acquire each slice throughout the MRI
acquisition. Defaults to `0.0`.
"""
sliceduration(img::ImageMeta{T,N,A,ImageProperties{format"NII"}}) where {T,N,A} = sliceduration(properties(img))
sliceduration(s::ImageStream) = sliceduration(properties(s))
sliceduration(p::ImageProperties) = getheader(p, "sliceduration", 0.0)
sliceduration(A::AbstractArray) = 0.0

"""
    slicestart(x) -> Int

Which slice corresponds to the first slice acquired during MRI acquisition
(i.e. not padded slices). Defaults to `1`.
"""
slicestart(img::ImageMeta{T,N,A,ImageProperties{format"NII"}}) where {T,N,A} = slicestart(properties(img))
slicestart(s::ImageStream) = slicestart(properties(s))
slicestart(p::ImageProperties) = getheader(p, "slicestart", 1)
slicestart(A::AbstractArray) = 1

"""
    sliceend(x)

Which slice corresponds to the last slice acquired during MRI acquisition
(i.e. not padded slices). Defaults to `size(x, slicedim(x))`.
"""
sliceend(img::ImageMeta{T,N,A,ImageProperties{format"NII"}}) where {T,N,A} =
    _sliceend(properties(img), size(img, slicedim(img)))
sliceend(s::ImageStream) = _sliceend(properties(s), size(s, slicedim(s)))
_sliceend(p::ImageProperties, slice_dim_size::Int) = getheader(p, "sliceend", slice_dim_size)
sliceend(p::ImageProperties) = getheader(p, "sliceend", 1)
sliceend(A::AbstractArray) = size(A, ndims(A))

"""
    qform(img)
"""
#=
qform(img::ImageMeta{T,N,A,ImageProperties{format"NII"}}) where {T,N,A} = qform(properties(img))
qform(s::ImageStream) = qform(properties(s))
qform(p::ImageProperties) = getheader(p, "qform", qform())
=#
function qform(s::Union{ImageStream,AbstractArray})
    qb, qc, qd, qfac = mat2quat(sform(s))
    return qform(qb, qc, qd, qoffsetx(s), qoffsety(s), qoffsetz(s),
                 ustrip.(pixelspacing(s))..., qfac)
end

# These are stored in the `properties["header"]["qoffset*"]` fields, so they can be used
# if desired but are not integrated into spacedirections because it's unlikely that we want
# to offset every single image axis by a couple of millimeters

qoffsetx(s::Union{ImageStream,AbstractArray}) =
    getheader(s, "qoffsetx", ustrip(axesoffsets(s, 1)))
qoffsety(s::Union{ImageStream,AbstractArray}) =
    getheader(s, "qoffsety", ustrip(axesoffsets(s, 2)))
qoffsetz(s::Union{ImageStream,AbstractArray}) =
    getheader(s, "qoffsetz", ustrip(axesoffsets(s, 3)))

"""
    sform(A)

The 4th column of the matrix is the offset of the affine matrix.
This is primarily included for the purpose of compatibility with DICOM formats, where the
"Image Position" stores the coordinates of the center of the first voxel
(see the [DICOM standard](http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.7.6.2.html#sect_C.7.6.2.1.1) for more details;
Note, these values should be in interpreted as 'mm').
"""
# may just drop sform as property and always grab from spacedirections in future
sform(img::Union{AbstractArray,ImageStream}) = _sform(img, spacedirections(img))

function _sform(img, sd::Tuple{Ax1,Ax2,Ax3}) where {Ax1,Ax2,Ax3}
    SMatrix{4,4,eltype(Ax1),16}([[sd[1]..., qoffsetx(img)]'
                                 [sd[2]..., qoffsety(img)]'
                                 [sd[3]..., qoffsetz(img)]'
                                 [0.0, 0.0, 0.0, 1.0]'])
end

function _sform(img, sd::Tuple{Ax1,Ax2}) where {Ax1,Ax2}
    SMatrix{4,4,eltype(Ax1),16}([[sd[1]..., qoffsetx(img)]'
                                 [sd[2]..., qoffsety(img)]'
                                 [0.0, 0.0, 0.0, 0.0]'
                                 [0.0, 0.0, 0.0, 1.0]'])
end

"""
    qformcode(x)


Code describing the orientation of the image in the scanner.
May be any of the following:

* Unkown
* Scanner_anat
"""
qformcode(img::ImageMeta{T,N,A,ImageProperties{format"NII"}}) where {T,N,A} = qformcode(properties(img))
qformcode(s::ImageStream) = qformcode(properties(s))
qformcode(p::ImageProperties) = getheader(p, "qformcode", :Unkown)
qformcode(A::AbstractArray) = :Unkown

"""
    sformcode(x)

Code describing the orientation of the image.
May be any of the following:

* Unkown
* Aligned_anat
* Talairach
* MNI152
"""
sformcode(img::ImageMeta{T,N,A,ImageProperties{format"NII"}}) where {T,N,A} = sformcode(properties(img))
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
toffset(s::ImageStream{T,N}) where {T,N} = ustrip(firstindex(axes(s, )))
