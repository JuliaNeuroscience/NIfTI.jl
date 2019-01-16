#########
# Types #
#########
function string_tuple(x::String, n::Int)
    a = codeunits(x)
    padding = zeros(UInt8, n-length(a))
    (a..., padding...)
end
string_tuple(x::AbstractString) = string_tuple(bytestring(x))

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
# keep order of NIfTI-1/2 fields same
function Nifti1Header(dim_info::Int8, dim::NTuple{8,Int16}, intent_p1::Float32,
                      intent_p2::Float32, intent_p3::Float32,
                      intent_code::Int16, datatype::Int16, bitpix::Int16,
                      slice_start::Int16, pixdim::NTuple{8,Float32},
                      vox_offset::Float32, scl_slope::Float32,
                      scl_inter::Float32, slice_end::Int16, slice_code::Int8,
                      xyzt_units::Int8, cal_max::Float32, cal_min::Float32,
                      slice_duration::Float32, toffset::Float32,
                      descrip::NTuple{80,UInt8}, aux_file::NTuple{24,UInt8},
                      qform_code::Int16, sform_code::Int16, quatern_b::Float32,
                      quatern_c::Float32, quatern_d::Float32,
                      qoffset_x::Float32, qoffset_y::Float32,
                      qoffset_z::Float32, srow_x::NTuple{4,Float32},
                      srow_y::NTuple{4,Float32}, srow_z::NTuple{4,Float32},
                      intent_name::NTuple{16,UInt8})
    Nifti1Header(SIZEOF_HDR1, string_tuple("", 10), string_tuple("", 18),
                 Int32(0), Int16(0), Int8(0), dim_info, dim, intent_p1,
                 intent_p2, intent_p3, intent_code, datatype, bitpix,
                 slice_start, pixdim, vox_offset, scl_slope, scl_inter,
                 slice_end, slice_code, xyzt_units, cal_max, cal_min,
                 slice_duration, toffset, Int32(0), Int32(0), descrip,
                 aux_file, qform_code, sform_code, quatern_b, quatern_c,
                 quatern_d,qoffset_x, qoffset_y, qoffset_z, srow_x, srow_y,
                 srow_z, intent_name, NP1_MAGIC)
end

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

# keep order of NIfTI-1/2 fields same
function Nifti2Header(dim_info::Int8, dim::NTuple{8,Int16}, intent_p1::Float64,
                      intent_p2::Float64, intent_p3::Float64,
                      intent_code::Int32, datatype::Int16, bitpix::Int16,
                      slice_start::Int64, pixdim::NTuple{8,Float64},
                      vox_offset::Int64, scl_slope::Float64,
                      scl_inter::Float64, slice_end::Int64, slice_code::Int32,
                      xyzt_units::Int32, cal_max::Float64, cal_min::Float64,
                      slice_duration::Float64, toffset::Float64,
                      descrip::NTuple{80,UInt8}, aux_file::NTuple{24,UInt8},
                      qform_code::Int32, sform_code::Int32, quatern_b::Float64,
                      quatern_c::Float64, quatern_d::Float64,
                      qoffset_x::Float64, qoffset_y::Float64,
                      qoffset_z::Float64, srow_x::NTuple{4,Float64},
                      srow_y::NTuple{4,Float64}, srow_z::NTuple{4,Float64},
                      intent_name::NTuple{16,UInt8})
    Nifti2Header(SIZEOF_HDR2, NP2_MAGIC, datatype, bitpix, dim, intent_p1,
                 intent_p2, intent_p3, pixdim, vox_offset, scl_slope,
                 scl_inter, cal_max, cal_min, slice_duration, toffset,
                 slice_start, slice_end, descrip, aux_file, qform_code,
                 sform_code, quatern_b, quatern_c, quatern_d, qoffset_x,
                 qoffset_y, qoffset_z, srow_x, srow_y, srow_z, intent_code,
                 intent_name, dim_info, string_tuple("", NTuple{15,UInt8}))
end

# common entry point for NIfTI-1 and 2
function NiftiHeader(dim_info, dim, intent_p1, intent_p2, intent_p3,
                     intent_code, datatype, bitpix, slice_start, pixdim,
                     vox_offset, scl_slope, scl_inter, slice_end, slice_code,
                     xyzt_units, cal_max, cal_min, slice_duration, toffset,
                     descrip, aux_file, qform_code, sform_code, quatern_b,
                     quatern_c, quatern_d, qoffset_x, qoffset_y,
                     qoffset_z, srow_x, srow_y, srow_z, intent_name; v::Int=1)
    if v == 1
        Nifti1Header(convert(Int8, dim_info), convert(NTuple{8,Int16}, dim),
                     convert(Float32,intent_p1), convert(Float32,intent_p2),
                     convert(Float32,intent_p3), convert(Int16, intent_code),
                     convert(Int16, datatype), convert(Int16, bitpix),
                     convert(Int16, slice_start), convert(NTuple{8,Float32}, pixdim),
                     convert(Float32, vox_offset), convert(Float32, scl_slope),
                     convert(Float32, scl_inter), convert(Int16, slice_end),
                     convert(Int8, slice_code), convert(Int8, xyzt_units),
                     convert(Float32, cal_max), convert(Float32, cal_min),
                     convert(Float32, slice_duration), convert(Float32, toffset),
                     convert(NTuple{80,UInt8}, descrip), convert(NTuple{24,UInt8}, aux_file),
                     convert(Int16, qform_code), convert(Int16, sform_code),
                     convert(Float32, quatern_b), convert(Float32, quatern_c),
                     convert(Float32, quatern_d), convert(Float32, qoffset_x),
                     convert(Float32, qoffset_y), convert(Float32, qoffset_z),
                     convert(NTuple{4,Float32}, srow_x), convert(NTuple{4,Float32}, srow_y),
                     convert(NTuple{4,Float32}, srow_z), convert(NTuple{16,UInt8},intent_name))
    elseif v == 2
        Nifti2Header(convert(Int8, dim_info), convert(NTuple{8,Int16}, dim),
                     convert(Float64,intent_p1), convert(Float64,intent_p2),
                     convert(Float64,intent_p3), convert(Int32, intent_code),
                     convert(Int16, datatype), convert(Int16, bitpix),
                     convert(Int64, slice_start), convert(NTuple{8,Float64}, pixdim),
                     convert(Int64, vox_offset), convert(Float64, scl_slope),
                     convert(Float64, scl_inter), convert(Int64, slice_end),
                     convert(Int32, slice_code), convert(Int32, xyzt_units),
                     convert(Float64, cal_max), convert(Float64, cal_min),
                     convert(Float64, slice_duration), convert(Float64, toffset),
                     convert(NTuple{80,UInt8}, descrip), convert(NTuple{24,UInt8}, aux_file),
                     convert(Int32, qform_code), convert(Int32, sform_code),
                     convert(Float64, quatern_b), convert(Float64, quatern_c),
                     convert(Float64, quatern_d), convert(Float64, qoffset_x),
                     convert(Float64, qoffset_y), convert(Float64, qoffset_z),
                     convert(NTuple{4,Float64}, srow_x), convert(NTuple{4,Float64}, srow_y),
                     convert(NTuple{4,Float64}, srow_z), convert(NTuple{16,UInt8},intent_name))
    end
end

##############
# Read/Write #
##############
# TODO ensure iostream is at proper position
# function base.write(s::NiftiStream)
#    if position(s.io) == 0
#        write(s.io, s.hdr)
#    else
#        @error "It appears the provided IOStream already has written a NiftiHeader."
#    end
#end

function define_packed(ty::DataType)
    packed_offsets = cumsum([sizeof(x) for x in ty.types])
    sz = pop!(packed_offsets)
    pushfirst!(packed_offsets, 0)

    @eval begin
        function niread(io::IO, ::Type{$ty}, needswap)
            bytes = read!(io, Array{UInt8}(undef, $sz...))
            hdr = $(Expr(:new, ty, [:(unsafe_load(convert(Ptr{$(ty.types[i])}, pointer(bytes)+$(packed_offsets[i])))) for i = 1:length(packed_offsets)]...,))
            if needswap
                return byteswap(hdr)
            end
            hdr
        end
        function niwrite(io::IO, x::$ty)
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

define_packed(Nifti1Header)
define_packed(Nifti2Header)

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



