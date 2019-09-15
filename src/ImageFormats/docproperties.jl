"""
    dataoffset(x) -> Int

Returns the IO stream offset to data given an type instance. Defaults to 0.
"""
dataoffset(x::Any) = getter(x, "dataoffset", Int, 1)

"""
    dataoffset!(x, val)

Set the the `data_offset` property.
"""
dataoffset!(x::Any, i::Integer) = setter!(x, "dataoffset", i, Int)

"""

    auxfiles(x) -> Vector{String}

Retrieves string for auxiliary file associated with the image.
"""
auxfiles(x::Any) = getter(x, "auxfiles", Vector{String}, [""])

"""
    auxfiles!(x, val)

Sets the `auxfiles` property. `val` should be a `String` or `Vector{String}`.
"""
auxfiles!(x::Any, i::Vector{String}) = setter!(x, "auxfiles", i, Vector{String})

"""
    srcfile(x) -> String

Retrieves the file name that the image comes from.
"""
srcfile(x::Any) = getter(x, "srcfile", String, "")

"""
    srcfile!(x, f::String)

Change `srcfile` property.
"""
srcfile!(x::Any, val::AbstractString) = setter!(x, "srcfile", val, String)

"""
    description(x) -> String

Retrieves description field that may say whatever you like.
"""
description(x::Any) = getter(x, "description", String, "")

"""
    description!(x, descrip::String)

Change description defined in an properties type.
"""
description!(x::Any, val::AbstractString) = setter!(x, "description", val, String)


"""
    freqdim(x) -> Int

Which spatial dimension (1, 2, or 3) corresponds to phase acquisition. If not
applicable to scan type defaults to `0`.
"""
freqdim(x::Any) = getter(x, "freqdim", Int, 0)

"""
    freqdim!(x, val)

Set the frequency dimension of `x` to `val`.
"""
freqdim!(x::Any, i::Integer) = setter!(x, "freqdim", i, Int)

"""
    slicedim(x) -> Int

Which dimension slices where acquired at throughout MRI acquisition.
"""
slicedim(x::Any) = getter(x, "slicedim", Int, 0)

"""
    slicedim!(x, val)

Set the slice dimension of `x` to `val`.
"""
slicedim!(x::Any, i::Integer) = setter!(x, "slicedim", i, Int)

"""
    phasedim(x) -> Int

Which spatial dimension (1, 2, or 3) corresponds to phase acquisition.
"""
phasedim(x::Any) = getter(x, "phasedim", Int, 0)

"""
    phasedim!(x, val)

Set the slice dimension of `x` to `val`.
"""
phasedim!(x::Any, i::Integer) = setter!(x, "phaseedim", i, Int)

"""
    slicestart(x) -> Int

Which slice corresponds to the first slice acquired during MRI acquisition
(i.e. not padded slices). Defaults to `1`.
"""
slicestart(x::Any) = getter(x, "slicestart", Int, 0)

"""
    slicestart!(x, val)

Set the first slice acquired during the MRI acquisition.
"""
slicestart!(x::Any, i::Integer) = setter!(x, "slicestart", i, Int)

"""
    sliceend(x) -> Int

Which slice corresponds to the last slice acquired during MRI acquisition
(i.e. not padded slices).
"""
sliceend(x::Any) = getter(x, "sliceend", Int, 1)

"""
    sliceend!(x, val)

Set the last slice acquired during the MRI acquisition.
"""
sliceend!(x::Any, i::Integer) = setter!(x, "sliceend", i, Int)

"""
    sliceduration(x) -> Float64

The amount of time necessary to acquire each slice throughout the MRI
acquisition.
"""
sliceduration(x::Any) = getter(x, "sliceduration", Float64, one(Float64))

"""
    sliceduration!(x, val)

Set the sliceduration of `x` to `val`.
"""
sliceduration!(x::Any, i::Real) = setter!(x, "sliceduration", i, Float64)

_caltype(x::AbstractArray{T}) where {T} = T
_caltype(x::Any) = Float64

"""
    calmax(x)

Specifies maximum element for display puproses. Defaults to the maximum of `x`.
"""
calmax(x::Any) = getter(x, "calmax", i -> _caltype(x), i -> maximum(i))

"""
    calmax!(x, val)

Change calmax defined in an `ImageProperties` type.
"""
calmax!(x::Any, val::Any) = setter!(x, "calmax", val, i -> _caltype(i))

"""

    calmin(x)

Specifies minimum element for display puproses. Defaults to the minimum of `x`.
"""
calmin(x::Any) = getter(x, "calmin", i -> _caltype(i), i -> minimum(i))

"""
    calmin!(x, val)

Set the calmin property.
"""
calmin!(x::Any, val::Any) = setter!(x, "calmin", val, i -> _caltype(i))

"""
    magicbytes(x) -> Vector{UInt8}

Retrieves the magicbytes associated with the file that produced an image
instance. Defaults to `[0x00]` if not found.
"""
magicbytes(x::Any) = getter(x, "magicbytes", Vector{UInt8}, [0x00])

"""
    magicbytes(x, val)

Sets the `magicbytes` property.
"""
magicbytes!(x::Any, val::Vector{UInt8}) = setter!(x, "magicbytes", val, Vector{UInt8})

"""
    sformcode(x) -> CoordinateSpace

Code describing the orientation of the image.

Should only be one of following (although others are allowed):

* UnkownSpace
* AnatomicalSpace
* TalairachSpace
* MNI152Space
"""
sformcode(x::Any) = getter(x, "sformcode", CoordinateSpace, UnkownSpace)

"""
    sformcode!(x, val)

Set the `sform` coordinate space of `x` to `val`.
"""
sformcode!(x::Any, val::CoordinateSpace) = setter!(x, "sformcode", val, CoordinateSpace)

"""
    qformcode(x) -> CoordinateSpace

Code describing the orientation of the image in the scanner.

Should only be one of the following (although others are allowed):

* UnkownSpace
* ScannerSpace
"""
qformcode(x::Any) = getter(x, "sformcode", CoordinateSpace, UnkownSpace)

"""
    qformcode!(x, val)

Set the `qfrom` coordinate space of `x` to `val`.
"""
qformcode!(x::Any, val::CoordinateSpace) = setter!(x, "qformcode", val, CoordinateSpace)

"""
    qform(x) -> MMatrix{4,4,Float64,16}
"""
qform(x::Any) = getter(x, "qform", MMatrix{4,4,Float64,16}, i->default_affinemat(i))

qform!(x::Any, val::AbstractMatrix) = setter!(x, "qform", val, MMatrix{4,4,Float64,16})

"""
    sform(x) -> MMatrix{4,4,Float64,16}

The 4th column of the matrix is the offset of the affine matrix.
This is primarily included for the purpose of compatibility with DICOM formats, where the
"Image Position" stores the coordinates of the center of the first voxel
(see the [DICOM standard](http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.7.6.2.html#sect_C.7.6.2.1.1) for more details;
Note, these values should be in interpreted as 'mm').
"""
sform(x::Any) = getter(x, "sform", MMatrix{4,4,Float64,16}, i->default_affinemat(i))

sform!(x::Any, val::AbstractMatrix) = setter!(x, "sform", val, MMatrix{4,4,Float64,16})


function default_affinemat(x::Any)
    if sformcode(x) === UnkownSpace
        if qformcode(x) === UnkownSpace
            return _default_affinemat(pixelspacing(x))
        else
            return qform(x)
        end
    else
        return sform(x)
    end
    default_affinemat(x)
end

function _default_affinemat(x::Any) where {N}
    _default_affinemat(spacedirections(x), pixelspacing(x))
end

function _default_affinemat(sd::NTuple{2,NTuple{2,T}}, ps::NTuple{2,T}) where T<:AbstractFloat
    MMatrix{4,4,Float64,16}(sd[1][1], sd[2][1], 0, 0,
                            sd[1][2], sd[2][2], 0, 0,
                                   0,        0, 0, 0,
                               ps[1],    ps[2], 0, 0)
end

function _default_affinemat(sd::NTuple{3,NTuple{3,T}}, ps::NTuple{3,T}) where T<:AbstractFloat
    MMatrix{4,4,Float64,16}(sd[1][1], sd[2][1], sd[3][1], 0,
                            sd[1][2], sd[2][2], sd[3][2], 0,
                            sd[1][3], sd[2][3], sd[3][3], 0,
                               ps[1],    ps[2],    ps[3], 0)
end

"""
    affinematrix(x) -> MMatrix{4,4,Float64,16}

Returns affine matrix. For an instance that returns `spacedirections` this is
the corresponding tuple converted to an array.
"""
affinematrix(x::Any) = getter(x, "affinematrix", MMatrix{4,4,Float64,16}, i -> default_affinemat(i))


affinematrix!(x::Any, val::AbstractMatrix) = setter!(x, "affinematrix", val, MMatrix{4,4,Float64,16})
