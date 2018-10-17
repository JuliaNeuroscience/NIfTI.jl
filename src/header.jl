function define_packed(ty::DataType)
    packed_offsets = cumsum([sizeof(x) for x in ty.types])
    sz = pop!(packed_offsets)
    pushfirst!(packed_offsets, 0)

    @eval begin
        function Base.read(io::IO, ::Type{$ty})
            bytes = read!(io, Array{UInt8}(undef, $sz...))
            hdr = $(Expr(:new, ty, [:(unsafe_load(convert(Ptr{$(ty.types[i])}, pointer(bytes)+$(packed_offsets[i])))) for i = 1:length(packed_offsets)]...,))
            if hdr.sizeof_hdr == ntoh(Int32(348))
                return byteswap(hdr), true
            end
            hdr, false
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

mutable struct NIfTI1Header
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
define_packed(NIfTI1Header)

# byteswapping

function byteswap(hdr::NIfTI1Header)
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

# Always in mm
voxel_size(header::NIfTI1Header) =
    [header.pixdim[i] * SPATIAL_UNIT_MULTIPLIERS[header.xyzt_units & Int8(3)] for i = 2:min(header.dim[1], 3)+1]

# Always in ms
time_step(header::NIfTI1Header) =
    header.pixdim[5] * TIME_UNIT_MULTIPLIERS[header.xyzt_units >> 3]


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

# Returns or sets dim_info as a tuple whose values are the frequency, phase, and slice dimensions
dim_info(header::NIfTI1Header) = (header.dim_info & int8(3), (header.dim_info >> 2) & int8(3),
    (header.dim_info >> 4) & int8(3))

dim_info(header::NIfTI1Header, dim_info::Tuple{T, T, T}) where {T<:Integer} =
    header.dim_info = to_dim_info(dim_info)

# Convert a NIfTI header to a 4x4 affine transformation matrix
function getaffine(h::NIfTI1Header)
    pixdim = convert(Tuple{Vararg{Float64}}, h.pixdim)
    # See documentation at
    # http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/qsform.html
    if h.sform_code > 0
        # Method 3
        return Float64[
            collect(h.srow_x)'
            collect(h.srow_y)'
            collect(h.srow_z)'
            0 0 0 1
        ]
    elseif h.qform_code > 0
        # Method 2
        # http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/quatern.html
        b::Float64 = h.quatern_b
        c::Float64 = h.quatern_c
        d::Float64 = h.quatern_d
        a = sqrt(1 - b*b - c*c - d*d)
        A = [
            a*a+b*b-c*c-d*d   2*b*c-2*a*d       2*b*d+2*a*c
            2*b*c+2*a*d       a*a+c*c-b*b-d*d   2*c*d-2*a*b
            2*b*d-2*a*c       2*c*d+2*a*b       a*a+d*d-c*c-b*b
        ]
        B = [
            pixdim[2]
            pixdim[3]
            pixdim[1]*pixdim[4]
        ]
        return vcat(hcat(A*B, [
            h.qoffset_x
            h.qoffset_y
            h.qoffset_z
        ]), [0 0 0 1])
    elseif h.qform_code == 0
        # Method 1
        return Float64[
             pixdim[2] 0         0          0 0
                       0 pixdim[3]          0 0
                       0         0  bitpix[4] 0
                       0         0          0 1
        ]
    end
end

# Set affine matrix of NIfTI header
function setaffine(h::NIfTI1Header, affine::Array{T,2}) where {T}
    size(affine, 1) == size(affine, 2) == 4 ||
        error("affine matrix must be 4x4")
    affine[4, 1] == affine[4, 2] == affine[4, 3] == 0 && affine[4, 4] == 1 ||
        error("last row of affine matrix must be [0 0 0 1]")
    h.qform_code = one(Int16)
    h.sform_code = one(Int16)
    h.pixdim = (zero(Float32), h.pixdim[2:end]...,)
    h.quatern_b = zero(Float32)
    h.quatern_c = zero(Float32)
    h.quatern_d = zero(Float32)
    h.qoffset_x = zero(Float32)
    h.qoffset_y = zero(Float32)
    h.qoffset_z = zero(Float32)
    h.srow_x = convert(Tuple{Vararg{Float32}}, (affine[1, :]...,))
    h.srow_y = convert(Tuple{Vararg{Float32}}, (affine[2, :]...,))
    h.srow_z = convert(Tuple{Vararg{Float32}}, (affine[3, :]...,))
    h
end
