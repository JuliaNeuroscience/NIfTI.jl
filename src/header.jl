
"""
    NIfTIHeader

* `dim_info`: MRI slice ordering.
* `dim::NTuple{8}`: Data array dimensions.
* `datatype`: Defines data type!
* `bitpix`: Number bits/voxel.
* `slice_start`: First slice index.
* `pixdim[8]`: Grid spacings.
* `vox_offset`: Offset into .nii file
* `xyzt_units`: Units of pixdim[1..4]
* `slice_duration: Time for 1 slice.
* `toffset`: Time axis shift.
* `slice_end`: Last slice index.
* `slice_code`: Slice timing order.

* `descrip::NTuple{80,UInt8}`: any text you like
* `aux_file::NTuple{24,UInt8}`: auxiliary filename.

* `qform_code`: NIFTI_XFORM_* code.
* `sform_code`: NIFTI_XFORM_* code.
* `quatern_b`: Quaternion b param.
* `quatern_c`: Quaternion c param.
* `quatern_d`: Quaternion d param.
* `qoffset_x`: Quaternion x shift
* `qoffset_y`: Quaternion y shift.
* `qoffset_z`: Quaternion z shift.
* `srow_x`: 1st row affine transform
* `srow_y`: 2nd row affine transform
* `srow_z`: 3rd row affine transform

* `intent_name[16]`: name' or meaning of data.
* `intent_p1`: 1st intent parameter.
* `intent_p2`: 2nd intent parameter.
* `intent_p3`: 3rd intent parameter.
* `intent_code`: NIFTI_INTENT_* code.

# We don't use these
* `scl_slope`: Data scaling: slope.
* `scl_inter`: Data scaling: offset.
* `cal_max`: Max display intensity
* `cal_min`: Min display intensity

"""
abstract type NIfTIHeader end

struct NIfTI1Header <: NIfTIHeader
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
    slice_code::Int8
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

struct NIfTI2Header <: NIfTIHeader
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

sizeof_hdr(::NIfTI1Header) = Int32(348)

sizeof_hdr(::NIfTI2Header) = Int32(540)

function write_header(io, x::NIfTIHeader)
    write(io, sizeof_hdr(x))
    for n in fieldnames(NIfTI1Header)
        write_type(io, getfield(x, n))
    end
end

write_type(io, x) = write(io, x)
write_type(io, x::Tuple) = write(io, Ref(x))

@generated function _write_header(io, x::T) where {T}
    out = Expr(:block)
    for i in 1:nfields(T)
        T_i = fieldtype(T, i)
        n = fieldname(T, i)
        if T_i <: Tuple
            push!(out.args, :(write(io, Ref(x.$n))))
        else
            push!(out.args, :(write(io, x.$n)))
        end
    end
    return out
end

@generated function read_header(io, ::Type{T}) where {T}
    out = Expr(:new, T)
    for p in fieldtypes(T)
        if p <: Tuple
            t = Expr(:tuple)
            for p_i in p.parameters
                push!(t.args, :(read(io, $p_i)))
            end
            push!(out.args, t)
        else
            push!(out.args, :(read(io, $p)))
        end
    end
    return out
end

@generated function read_header_swap(io, ::Type{T}) where {T}
    out = Expr(:new, T)
    for p in fieldtypes(T)
        if p <: Tuple{Vararg{UInt}}
            t = Expr(:tuple)
            for p_i in p.parameters
                push!(t.args, :(read(io, $p_i)))
            end
            push!(out.args, t)
        elseif p <: Tuple
            t = Expr(:tuple)
            for p_i in p.parameters
                push!(t.args, :(bswap(read(io, $p_i))))
            end
            push!(out.args, t)
        else
            push!(out.args, :(bswap(read(io, $p))))
        end
    end
    return out
end

niopen(file::AbstractString, mode::AbstractString) = niopen(open(file, mode))
@inline function niopen(io)
    b1, b2 = read(io, 2)
    if (b1, b2) === (0x1F, 0x8B)
        seek(io, 0)
        return niopen(GzipDecompressorStream(io))
    else
        b3, b4 = read(io, 2)
        if (b1,b2,b3,b4) === (0x5c,0x01,0x00,0x00)
            return io, read_header(io, NIfTI1Header), false
        elseif (b1,b2,b3,b4) === (0x1c,0x02,0x00,0x00)
            return io, read_header(io, NIfTI2Header), false
        elseif (b1,b2,b3,b4) === (0x00,0x00,0x01,0x5c)
            return io, read_header_swap(io, NIfTI1Header), true
        elseif (b1,b2,b3,b4) === (0x00,0x00,0x02,0x1c)
            return io, read_header_swap(io, NIfTI2Header), true
        else
            error("This is not a NIfTI-1 file")
        end
    end
end

