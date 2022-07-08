

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

"""
    freqdim(x)::Int

Returns the frequency dimension.
"""
freqdim(x::NIfTI1Header) = Int(getfield(x, :dim_info) & 0x03)
freqdim(x) = freqdim(header(x))

"""
    phasedim(x)::Int

Returns the phase dimension.
"""
phasedim(x::NIfTI1Header) = Int((getfield(x, :dim_info) >> 2) & 0x03)
phasedim(x) = phasedim(header(x))

"""
    slicedim(x)::Int

Returns the slice dimension.
"""
slicedim(x::NIfTI1Header) = Int((getfield(x, :dim_info) >> 4) & 0x03)
slicedim(x) = slicedim(header(x))

"""
    slice_start(x)::Int

Which slice corresponds to the first slice acquired during MRI acquisition (i.e. not padded slices).
"""
slice_start(x::NIfTI1Header) = Int(getfield(x, :slice_start)) + 1
slice_start(x) = slice_start(header(x))

"""
    slice_end(x)::Int

Which slice corresponds to the last slice acquired during MRI acquisition (i.e. not padded slices).
"""
slice_end(x::NIfTI1Header) = Int(getfield(x, :slice_end)) + 1
slice_end(x) = slice_end(header(x))

"""
    slice_duration(x)

Time to acquire one slice
"""
slice_duration(x::NIfTI1Header) = getfield(x, :slice_duration)
slice_duration(x) = slice_duration(header(x))

