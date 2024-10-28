
function define_packed(ty::DataType)
    packed_offsets = cumsum([sizeof(x) for x in ty.types])
    sz = pop!(packed_offsets)
    pushfirst!(packed_offsets, 0)

    @eval begin
        function Base.read(io::IO, ::Type{$ty})
            bytes = read!(io, Array{UInt8}(undef, $sz...))
            hdr = $(Expr(:new, ty, [:(unsafe_load(convert(Ptr{$(ty.types[i])}, pointer(bytes) + $(packed_offsets[i])))) for i = 1:length(packed_offsets)]...,))
            if hdr.sizeof_hdr == ntoh(Int32(348))
                return byteswap(hdr), true
            end
            hdr, false
        end
        function Base.write(io::IO, x::$ty)
            bytes = UInt8[]
            for name in fieldnames($ty)
                append!(bytes, reinterpret(UInt8, [getfield(x, name)]))
            end
            write(io, bytes)
            $sz
        end
    end
    nothing
end

abstract type NIfTIHeader end

mutable struct NIfTI1Header <: NIfTIHeader
    sizeof_hdr::Int32

    data_type::NTuple{10,UInt8}
    db_name::NTuple{18,UInt8}
    extents::Int32
    session_error::Int16
    regular::Int8

    dim_info::UInt8
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
    slice_code::UInt8
    xyzt_units::UInt8
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
define_packed(NIfTI1Header)

mutable struct NIfTI2Header <: NIfTIHeader
    sizeof_hdr::Int32
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
    xyzt_units::Int32
    intent_code::Int32
    intent_name::NTuple{16,UInt8}
    dim_info::UInt8
    unused_str::NTuple{15,UInt8}
end
define_packed(NIfTI2Header)


# byteswapping

function byteswap(hdr::NIfTIHeader)
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

"""
    NIfTI.slice_start(x)::Int

Which slice corresponds to the first slice acquired during MRI acquisition (i.e. not padded slices).
"""
slice_start(x::NIfTIHeader) = Int(getfield(x, :slice_start)) + 1
slice_start(x) = slice_start(header(x))

"""
    NIfTI.slice_end(x)::Int

Which slice corresponds to the last slice acquired during MRI acquisition (i.e. not padded slices).
"""
slice_end(x::NIfTIHeader) = Int(getfield(x, :slice_end)) + 1
slice_end(x) = slice_end(header(x))

"""
    NIfTI.slice_duration(x)

Time to acquire one slice.
"""
slice_duration(x::NIfTIHeader) = getfield(x, :slice_duration)
slice_duration(x) = slice_duration(header(x))

"""
    NIfTI.sdims(img)

Return the number of spatial dimensions in the image.
"""
sdims(x::NIfTIHeader) = min(Int(getfield(getfield(x, :dim), 1)), 3)
sdims(x) = sdims(header(x))

# Conversion factors to mm/ms
# http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/xyzt_units.html
# 1 => NIfTI_UNITS_METER
# 2 => NIfTI_UNITS_MM
# 3 => NIfTI_UNITS_MICRON
function to_spatial_multiplier(xyzt_units::UInt8)
    i = xyzt_units & 0x03
    if i === 0x01
        return 1000.0f0
    elseif i === 0x02
        return 1.0f0
    else
        return 0.001f0
    end
end

"""
    NIfTI.voxel_size(header::NIfTIHeader)

Get the voxel size **in mm** from a `NIfTIHeader`.
"""
function voxel_size(x::NIfTIHeader)
    sd = sdims(x)
    if sd === 0
        return ()
    else
        sconvert = to_spatial_multiplier(getfield(x, :xyzt_units))
        if sd === 3
            pd = getfield(x, :pixdim)
            return (
                getfield(pd, 2) * sconvert,
                getfield(pd, 3) * sconvert,
                getfield(pd, 4) * sconvert
            )
        elseif sd === 2
            pd = getfield(x, :pixdim)
            return (getfield(pd, 2) * sconvert, getfield(pd, 3) * sconvert,)
        else # sd === 1
            pd = getfield(x, :pixdim)
            return (getfield(pd, 2) * sconvert,)
        end
    end
end

function to_dim_info(dim_info::Tuple{Integer,Integer,Integer})
    if dim_info[1] > 3 || dim_info[1] < 0
        error("Invalid frequency dimension $(dim_info[1])")
    elseif dim_info[2] > 3 || dim_info[2] < 0
        error("Invalid phase dimension $(dim_info[2])")
    elseif dim_info[3] > 3 || dim_info[3] < 0
        error("Invalid slice dimension $(dim_info[3])")
    end

    return UInt8(dim_info[1] | (dim_info[2] << 2) | (dim_info[3] << 4))
end

# Returns or sets dim_info as a tuple whose values are the frequency, phase, and slice dimensions
function dim_info(hdr::NIfTIHeader)
    return (hdr.dim_info & 0x03, (hdr.dim_info >> 0x02) & 0x03, (hdr.dim_info >> 0x04) & 0x03)
end
function dim_info(header::NIfTIHeader, dim_info::Tuple{T,T,T}) where {T<:Integer}
    header.dim_info = to_dim_info(dim_info)
end

"""
    NIfTI.freqdim(img)::Int

Returns the frequency dimension associated with with `img`. `img` is an image or a collection of
image related metadata. 
"""
freqdim(x::NIfTIHeader) = Int(getfield(x, :dim_info) & 0x03)
freqdim(x) = freqdim(header(x))

"""
    NIfTI.phasedim(img)::Int

Returns the phase dimension associated with `img`. `img` is an image or a collection of
image related metadata.
"""
phasedim(x::NIfTIHeader) = Int((getfield(x, :dim_info) >> 0x02) & 0x03)
phasedim(x) = phasedim(header(x))

"""
    NIfTI.slicedim(img)::Int

Returns the slice dimension associated with `img`. `img` is an image or a collection of
image related metadata.
"""
slicedim(x::NIfTIHeader) = Int((getfield(x, :dim_info) >> 0x04) & 0x03)
slicedim(x) = slicedim(header(x))

