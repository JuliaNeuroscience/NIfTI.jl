# dim_info() -> diminfo() in order to avoid confusion with dim_info field
# voxel_size() -> pixelspacing() to for consistency with Image.jl

function define_packed(ty::DataType)
    packed_offsets = cumsum([sizeof(x) for x in ty.types])
    sz = pop!(packed_offsets)
    pushfirst!(packed_offsets, 0)

    @eval begin
        function Base.read(io::IO, ::Type{$ty}, needswap)
            bytes = read!(io, Array{UInt8}(undef, $sz...))
            hdr = $(Expr(:new, ty, [:(unsafe_load(convert(Ptr{$(ty.types[i])}, pointer(bytes)+$(packed_offsets[i])))) for i = 1:length(packed_offsets)]...,))
            if needswap
                return byteswap(hdr)
            end
            hdr
        end
        function Base.write(io::IO, x::$ty)
            bytes = UInt8[]
            for name in fieldnames($ty)
                append!(bytes, reinterpret(UInt8, [getfield(x,name)]))
            end
            write(io, bytes)
            $sz
        end
    end
    nothing
end

abstract type NiftiHeader end

mutable struct Nifti1Header <: NiftiHeader
    sizeof_hdr::Int32
    data_type::NTuple{10,UInt8}
    db_name::NTuple{18,UInt8}
    extents::Int32
    session_error::Int16
    regular::Int8
    dim_info::Int8
    dim::NTuple{8,Int16}
    intent_p1::Float32
    intent_p2::Float32
    intent_p3::Float32
    intent_code::Int16
    datatype::Int16
    bitpix::Int16
    slice_start::Int16
    pixdim::NTuple{8,Float32}
    vox_offset::Float32
    scl_slope::Float32
    scl_inter::Float32
    slice_end::Int16
    slice_code::Int8
    xyzt_units::Int8
    cal_max::Float32
    cal_min::Float32
    slice_duration::Float32
    toffset::Float32
    glmax::Int32
    glmin::Int32
    descrip::NTuple{80,UInt8}
    aux_file::NTuple{24,UInt8}
    qform_code::Int16
    sform_code::Int16
    quatern_b::Float32
    quatern_c::Float32
    quatern_d::Float32
    qoffset_x::Float32
    qoffset_y::Float32
    qoffset_z::Float32
    srow_x::NTuple{4,Float32}
    srow_y::NTuple{4,Float32}
    srow_z::NTuple{4,Float32}
    intent_name::NTuple{16,UInt8}
    magic::NTuple{4,UInt8}
end

define_packed(Nifti1Header)

mutable struct Nifti2Header <: NiftiHeader
    # must be 540
    sizeof_hdr::Int32
    # This placement is intentionally different from Nifti1
    magic::NTuple{8,UInt8}
    datatype::Int16
    bitpix::Int16
    dim::NTuple{8,Int64}
    intent_p1::Float64
    intent_p2::Float64
    intent_p3::Float64
    pixdim::NTuple{8,Float64}
    vox_offset::Int64
    scl_slope::Float64
    scl_inter::Float64
    cal_max::Float64
    cal_min::Float64
    slice_duration::Float64
    toffset::Float64
    slice_start::Int64
    slice_end::Int64
    descrip::NTuple{80,UInt8}
    aux_file::NTuple{24,UInt8}
    qform_code::Int32
    sform_code::Int32
    quatern_b::Float64
    quatern_c::Float64
    quatern_d::Float64
    qoffset_x::Float64
    qoffset_y::Float64
    qoffset_z::Float64
    srow_x::NTuple{4,Float64}
    srow_y::NTuple{4,Float64}
    srow_z::NTuple{4,Float64}
    slice_code::Int32
    xyzt_units::UInt32
    intent_code::Int32
    intent_name::NTuple{16,UInt8}
    dim_info::Int8
    unused_str::NTuple{15,UInt8}
end

define_packed(Nifti2Header)

# byteswapping
function byteswap(hdr::NiftiHeader)
    for fn in fieldnames(typeof(hdr))
        val = getfield(hdr, fn)
        if isa(val, Number) && sizeof(val) > 1
            setfield!(hdr, fn, ntoh(val))
        elseif isa(val, NTuple) && sizeof(eltype(val)) > 1
            setfield!(hdr, fn, map(ntoh, val))
        end
    end
    hdr
end

# NiftiHeader API
# - sdim
# - pixelspacing
# - nimages
# - timedim

sdim(hdr::NiftiHeader) = min(hdr.dim[1], 3)

function pixelspacing(hdr::NiftiHeader, standardize::Bool=false)
    [hdr.pixdim[i] * NIFTI_UNITS[hdr.xyzt_units & 0x07] for i = 2:sdim(hdr)+1]
end

nimages(hdr::NiftiHeader) = hdr.dim[1] < 5 ? 0 : hdr.dim[5]

function timedim(hdr::NiftiHeader)
    uu = NIFTI_UNITS[hdr.xyzt_units & 0x38]
    range(0.7, step=hdr.pixdim[5], length=hdr.dim[5]) * uu
end

"""
    diminfo(hdr::NiftiHeader) -> NTuple{N,Int}

Returns or sets dim_info as a tuple whose values are:
* frequency dimensions
* phase dimensions
* slice dimensions
"""
function to_dim_info(dim_info::Tuple{Integer,Integer,Integer})
    if dim_info[1] > 3 || dim_info[1] < 0
        error("Invalid frequency dimension $(dim_info[1])")
    elseif dim_info[2] > 3 || dim_info[2] < 0
        error("Invalid phase dimension $(dim_info[2])")
    elseif dim_info[3] > 3 || dim_info[3] < 0
        error("Invalid slice dimension $(dim_info[3])")
    end
    Int8(dim_info[1] | (dim_info[2] << 2) | (dim_info[3] << 4))
end

diminfo(hdr::NiftiHeader) = (hdr.dim_info & int8(3),
                             (hdr.dim_info >> 2) & int8(3),
                             (hdr.dim_info >> 4) & int8(3))

diminfo(hdr::NiftiHeader, dim_info::Tuple{T, T, T}) where {T<:Integer} =
    hdr.dim_info = to_diminfo(dim_info)

"""
    getdatatype(hdr::NiftiHeader) -> ::Type

Interpret NiftiHeader data_type code into Julia Type
"""
function getdatatype(hdr::NiftiHeader)
    t = get(NIFTI_DT_BITSTYPES, hdr.data_type, nothing)
    if t == nothing
        error("Unsupported data type $(hdr.data_type)")
    end
    t
end

### tools for setting NiftiHeader ###
"""
    setdim!(img) -> hdr::NiftiHeader

Set the dimensions of a NiftiHeader based on some AbstractArray
"""
function setdim!(hdr::NiftiHeader, img::AbstractArray)
    hdr.dim = ones(Int16, 8)
    hdr.dim[1] = ndims(img)
    hdr.dim[2:dim[1]+1] = [size(img)...]
end

"""
    setdatatype!(hdr::NiftiHeader, img::AbstractArray)

Sets datatype code for NIfTI header based on eltype of img
"""
# TODO double check data_type vs datatype is set properly here and elsewhere
function setdatatype!(hdr::NiftiHeader, img::AbstractArray)
    t = eltype(img)
    code = get(NIFTI_DT_BITSTYPES_REVERSE, t, nothing)
    if t == NIFTI_DT_BITSTYPES_REVERSE[t]
        error("Unsupported data type $t")
    else
        hdr.data_type = code
    end
end

# TODO
function setpixdim!(hdr::NiftiHeader, img::ImageMeta)
    hdr.pixdim = ones(Int16, 8)
    hdr.pixdim[1:ndims(img)] = pixelspacing(img)
end

"""
    setintent!(hdr::NiftiHeader, img::ImageMeta)

"""
function setintent!(hdr::NiftiHeader, img::ImageMeta)
    hdr.intent_code = get(NIFTI_INTENT_REVERS,
                          img.properties["header"]["intent_name"], 0)
    hdr.intent_name = img.properties["header"]["intent_name"]
    hdr.intent_p1 = img.properties["header"]["intent_param"][1]
    hdr.intent_p2 = img.properties["header"]["intent_param"][2]
    hdr.intent_p3 = img.properties["header"]["intent_param"][3]
end

# TODO test
# only handles cases where all spatial dimensions are same units right now
function setunits!(hdr::NiftiHeader, img::ImageMeta)
    su = unit(pixelspacing(img)[1])
    tu = timaxis(img) == nothing ? 0 : unit(timedim(img))
    hdr.xyzt_units = ((NIFTI_UNITS_REVERSE[su] & 0x07) &
                      (NIFTI_UNITS_REVERSE[tu] & 0x38))
end

# TODO I imagine these may be consolidated into a more useful parameter
# specific to Julia. Until it's figured out it just copies them from 
# property["header"]
# - reverse slice_code
function setslicing!(hdr::NiftiHeader, img::ImageMeta)
    hdr.slice_start = img.properties["header"]["slice_end"]
    hdr.slice_end = img.properties["header"]["slice_end"]
    hdr.slice_code = img.properties["header"]["slice_code"]
    hdr.slice_duration = img.properties["header"]["slice_duration"]
    hdr.toffset = img.properties["header"]["toffset"]
end

# TODO should get this from spacedirections() and image orientation
function setaffine!(h::NiftiHeader, img::ImageMeta)
    sform = spacedirections(img)
    #size(affine, 1) == size(affine, 2) == 4 ||
    #    error("affine matrix must be 4x4")
    #affine[4, 1] == affine[4, 2] == affine[4, 3] == 0 && affine[4, 4] == 1 ||
    #    error("last row of affine matrix must be [0 0 0 1]")
    h.qform_code = one(Int16)
    h.sform_code = one(Int16)
    h.quatern_b = zero(Float32)
    h.quatern_c = zero(Float32)
    h.quatern_d = zero(Float32)
    h.qoffset_x = zero(Float32)
    h.qoffset_y = zero(Float32)
    h.qoffset_z = zero(Float32)
    # not just Tuple{Vararg{Float32}} b/c of NIfTI-2
    h.srow_x = convert(typeof(h.srow_x), sform[1])
    h.srow_y = convert(typeof(h.srow_y), sform[2])
    h.srow_z = convert(typeof(h.srow_z), sform[3])
    h
end

# Gets the size of a type in bits
nibitpix(t::Type) = Int16(sizeof(t)*8)

"""
    getintent(hdr::NiftiHeader) -> Type <: AbstractIntent

Gets intent type from intent_code in NIfTI header
"""
getintent(hdr::NiftiHeader) = get(NIFTI_INTENT, hdr.intent_code, 0)
