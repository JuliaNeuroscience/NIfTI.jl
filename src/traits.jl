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
    scaleslope
"""
scaleslope(img::ImageMeta{T,N,A,ImageProperties{:NII}}) where {T,N,A} = scaleslope(properties(img))
scaleslope(s::ImageStream) = scaleslope(properties(s))
scaleslope(p::ImageProperties) = getheader(p, "scaleslope", zero(Float64))::Float64
scaleslope(A::AbstractArray) = 0.0

"""
    scaleintercept
"""
scaleintercept(img::ImageMeta{T,N,A,ImageProperties{:NII}}) where {T,N,A} = scaleintercept(properties(img))
scaleintercept(s::ImageStream) = scaleintercept(properties(s))
scaleintercept(p::ImageProperties) = getheader(p, "scaleintercept", zero(Float64))::Float64
scaleintercept(A::AbstractArray) = 0.0

# dimension info for nifti header
diminfo(img::ImageMeta{T,N,A,ImageProperties{:NII}}) where {T,N,A} = diminfo(properties(img))
diminfo(A::AbstractArray) = Int8(0)
diminfo(s::ImageStream) = diminfo(properties(s))::Int8
diminfo(p::ImageProperties) = getheader(p, "diminfo", zero(Int8))

"""
    frequencydim

```jldoctest
julia> using NIfTI, ImageMetadata

julia> p = ImageProperties{:NII}();

julia> p["header"] = Dict{String,Any}("diminfo" => Int8(57))

julia> img = ImageMeta(rand(4,4), p)

julia> frequencydim(img)
1

julia> phasedim(img)
2

julia> slicedim(img)
3
```
"""
frequencydim(s) = diminfo(s)::Int8 & Int8(3)

"""
    phasedim
"""
phasedim(s) = (diminfo(s)::Int8 >> 2) & Int8(3)

"""
    slicedim
"""
slicedim(s) = diminfo(s)::Int8 >> 4

"""
    slicecode
"""
slicecode(img::ImageMeta{T,N,A,ImageProperties{:NII}}) where {T,N,A} = slicecode(properties(img))
slicecode(s::ImageStream) = slicecode(properties(s))
slicecode(p::ImageProperties) = getheader(p, "slicecode", "Unkown")
slicecode(A::AbstractArray) = "Unkown"

"""
    sliceduration
"""
sliceduration(img::ImageMeta{T,N,A,ImageProperties{:NII}}) where {T,N,A} = sliceduration(properties(img))
sliceduration(s::ImageStream) = sliceduration(properties(s))
sliceduration(p::ImageProperties) = getheader(p, "sliceduration", 0.0)
sliceduration(A::AbstractArray) = 0.0

"""
    slicestart
"""
slicestart(img::ImageMeta{T,N,A,ImageProperties{:NII}}) where {T,N,A} = slicestart(properties(img))
slicestart(s::ImageStream) = slicestart(properties(s))
slicestart(p::ImageProperties) = getheader(p, "slicestart", 0)
slicestart(A::AbstractArray) = 0

"""
    sliceend
"""
sliceend(img::ImageMeta{T,N,A,ImageProperties{:NII}}) where {T,N,A} = sliceend(properties(img))
sliceend(s::ImageStream) = sliceend(properties(s))
sliceend(p::ImageProperties) = getheader(p, "sliceend", 0)
sliceend(A::AbstractArray) = 0

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
