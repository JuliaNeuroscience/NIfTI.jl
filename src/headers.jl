
function NIfTI1Header(img::MetadataArray{T}, extensions::Vector{Extension}) where {T}
    hdr = get(metadata(x), :header, nothing)
    hdr isa NIfTI1Header ? hdr : NIfTI1Header(parent(img), extensions)
end

function NIfTI1Header(
    x;
    extensions::Vector{Extension}=extensions(x),
    intent_parameters=intent_parameters(x),
    intent_code=intent_code(x),
    intent_name=intent_name(x),
    slice_start=slice_start(x),
    slice_end=slice_start(x),
    slice_duration=slice_duration(x),
    slice_code=slice_code(x),
    orientation::Union{AbstractMatrix{Float32},Nothing}=nothing,
    dim_info=nothing,
    freqdim=nothing,
    slicedim=nothing,
    phasedim=nothing
)
    T = eltype(x)
    # Fields specified as UNUSED in NIfTI1 spec
    data_type = ""
    db_name = ""
    extents = Int32(0)
    session_error  = Int16(0)
    regular = Int8(0)
    glmax = Int32(0)
    glmin = Int16(0)
    # The frequency encoding, phase encoding and slice dimensions.
    dim_info = (0, 0, 0)


    # The size of each voxel and the time step. These are formulated in mm unless xyzt_units is
    # explicitly specified
    voxel_size = (1f0, 1f0, 1f0)
    time_step = 0f0
    xyzt_units = Int8(18)
    # Slope and intercept by which volume shoudl be scaled
    scl_slope = 1f0
    scl_inter = 0f0
    # These describe how data should be scaled when displayed on the screen. They are probably
    # rarely used
    cal_max = 0f0
    cal_min = 0f0
    # Indicates a non-zero start point for time axis
    toffset = 0f0
    # "Any text you like"
    descrip = ""
    # Name of auxiliary file
    aux_file = ""

    # Parameters for Method 2. See the NIfTI spec
    qfac=0f0
    quatern_b=0f0
    quatern_c=0f0
    quatern_d=0f0
    qoffset_x=0f0
    qoffset_y=0f0
    qoffset_z=0f0
    # Orientation matrix for Method 3

    method2 = qfac != 0 || quatern_b != 0 || quatern_c != 0 || quatern_d != 0 || qoffset_x != 0 ||
        qoffset_y != 0 || qoffset_z != 0
    method3 = orientation !== nothing

    if method2 && method3
        error("Orientation parameters for Method 2 and Method 3 are mutually exclusive")
    end

    if method3
        if size(orientation) != (3, 4)
            error("Orientation matrix must be of dimensions (3, 4)")
        end
    else
        orientation = zeros(Float32, 3, 4)
    end

    if slice_start == 0 && slice_end == 0 && dim_info[3] != 0
        slice_start = 0
        slice_end = size(img, dim_info[3]) - 1
    end
    if isempty(extensions)
        voxoffset = Int32(352) + SIZEOF_HDR1
    else
        voxoffset = SIZEOF_HDR1
        for e in extensions
            voxoffset += esize(e)
        end
    end

    if dim_info === nothing
        dim_info = to_dim_info(
            freqdim === nothing ? freqdim(x) : freqdim,
            phasedim === nothing ? phasedim(x) : phasedim,
            slicedim === nothing ? slicedim(x) : slicedim
        )
    end

    NIfTI1Header(
        SIZEOF_HDR1,
        string_tuple(data_type, 10),
        string_tuple(db_name, 18),
        extents,
        session_error,
        regular,
        dim_info,
        to_dim_i16(size(img)),
        intent_parameters,
        intent_code,
        eltype_to_int(T),
        nibitpix(T),
        slice_start,
        (qfac, voxel_size..., time_step, 0, 0, 0),
        voxoffset,
        scl_slope,
        scl_inter,
        slice_end,
        slice_code,
        xyzt_units,
        cal_max,
        cal_min,
        slice_duration,
        toffset,
        glmax,
        glmin,
        string_tuple(descrip, 80),
        string_tuple(aux_file, 24),
        (method2 || method3),
        method3,
        quatern_b,
        quatern_c,
        quatern_d,
        qoffset_x,
        qoffset_y,
        qoffset_z,
        (orientation[1, :]...,),
        (orientation[2, :]...,),
        (orientation[3, :]...,),
        string_tuple(intent_name, 16),
        NP1_MAGIC
    )
end


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

function to_dim_info(fdim, pdim, sdim)
    if fdim > 3 || fdim < 0
        error("Invalid frequency dimension $(fdim)")
    elseif pdim > 3 || pdim < 0
        error("Invalid phase dimension $(pdim)")
    elseif sdim > 3 || sdim < 0
        error("Invalid slice dimension $(sdim)")
    end
    return Int8(fdim | (pdim << 2) | (sdim << 4))
end

# Returns or sets dim_info as a tuple whose values are the frequency, phase, and slice dimensions
function dim_info(header::NIfTIHeader)
    return (
        header.dim_info & int8(3),
        (header.dim_info >> 2) & int8(3),
        (header.dim_info >> 4) & int8(3)
    )
end
function dim_info(header::NIfTIHeader, dim_info::Tuple{T, T, T}) where {T<:Integer}
    header.dim_info = to_dim_info(dim_info)
end


"""
    NIfTI.freqdim(x)::Int

Returns the frequency dimension.
"""
freqdim(x::NIfTIHeader) = Int(getfield(x, :dim_info) & 0x03)
freqdim(x) = freqdim(header(x))

"""
    NIfTI.phasedim(x) -> Int

Returns the phase dimension.
"""
phasedim(x::NIfTIHeader) = Int((getfield(x, :dim_info) >> 2) & 0x03)
phasedim(x) = phasedim(header(x))

"""
    NIfTI.slicedim(x) -> Int

Returns the slice dimension.
"""
slicedim(x::NIfTIHeader) = Int((getfield(x, :dim_info) >> 4) & 0x03)
slicedim(x) = Int(slicedim(header(x)))

"""
    NIfTI.slice_start(x) -> Int

Which slice corresponds to the first slice acquired during MRI acquisition (i.e. not padded slices).
Note that this corresponds to where the slice is located using 1-based indexing.
"""
slice_start(x::NIfTIHeader) = Int(getfield(x, :slice_start)) + 1
slice_start(x) = Int(metadata(x, :slice_start, 1))

"""
    NIfTI.slice_end(x) -> Int

Which slice corresponds to the last slice acquired during MRI acquisition (i.e. not padded slices).
Note that this corresponds to where the slice is located using 1-based indexing.
"""
slice_end(x::NIfTIHeader) = Int(getfield(x, :slice_end)) + 1
slice_end(x) = Int(metadata(x, :slice_end, 1))

"""
    NIfTI.slice_duration(x) -> Union{Float32, Float64}

Time to acquire one slice.
"""
slice_duration(x::NIfTIHeader) = getfield(x, :slice_duration)
function slice_duration(x)
    sd = metadata(x, :slice_duration, 0f0)
    sd isa Union{Float32, Float64} ? Float32(sd) : sd
end

# Describes data contained in the file; for valid values, see
# http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/group__NIfTI1__INTENT__CODES.html
# Information about which slices were acquired, in case the volume has been padded
"""
    NIfTI.slice_code(x) -> Union{Int8, Int32}

Code corresponding to slice timing order.
"""
slice_code(x::NIfTIHeader) = getfield(x, :slice_code)
function slice_code(x)
    out = metadata(x, :slice_code, Int8(0))
    out isa Union{Int8, Int32} ? out : Int32(out)
end

"""
    NIfTI.qform_code(data) -> Union{Int16, Int32}

Integer code corresponding to the quaternion based coordinate system used for `data`.
This is described in more detail [here](https://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/qsform.html/document_view).
"""
qform_code(x::NIfTIHeader) = getfield(x, :qform_code)
function qform_code(x)
    out = metadata(x, :qform_code, Int16(0))
    if out isa Union{Int16, Int32}
        return out
    elseif out isa Symbol
        return xform_to_code(out)
    else
        return Int16(out)
    end
end

"""
    NIfTI.sform_code(data) -> Union{Int16, Int32}

Integer code corresponding to the standard space based coordinate system used for `data`.
This is described in more detail [here](https://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/qsform.html/document_view).
"""
sform_code(x::NIfTIHeader) = getfield(x, :sform_code)
function sform_code(x)
    out = metadata(x, :sform_code, Int16(0))
    if out isa Union{Int16, Int32}
        return out
    elseif out isa Symbol
        return xform_to_code(out)
    else
        return Int16(out)
    end
end

@inline function xform_to_code(x::Symbol)
    if x === :ScannerSpace
        return Int16(1)
    elseif x === :AnatomicalSpace
        return Int16(2)
    elseif x === :TalairachSpace
        return Int16(3)
    elseif x === :MNI152Space
        return Int16(4)
    elseif x === :OtherTemplate
        return Int16(5)
    else
        return Int16(0)
    end
end

"""
    NIfTI.voxel_size(header::NIfTIHeader)

Get the voxel size **in mm** from a `NIfTIHeader`.
"""
function voxel_size(header::NIfTIHeader)
    [header.pixdim[i] * SPATIAL_UNIT_MULTIPLIERS[header.xyzt_units & Int8(3)] for i = 2:min(header.dim[1], 3)+1]
end

