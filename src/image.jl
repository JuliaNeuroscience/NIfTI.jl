mutable struct NIVolume{T<:Number,N,R} <: AbstractArray{T,N}
    header::NIfTI1Header
    extensions::Vector{NIfTI1Extension}
    raw::R
end

NIVolume(header::NIfTI1Header, extensions::Vector{NIfTI1Extension}, raw::R) where {R}=
    niupdate(new(header, extensions, raw))

NIVolume(header::NIfTI1Header, extensions::Vector{NIfTI1Extension}, raw::AbstractArray{T,N}) where {T<:Number,N} =
    NIVolume{typeof(one(T)*1f0+1f0),N,typeof(raw)}(header, extensions, raw)

NIVolume(header::NIfTI1Header, raw::AbstractArray{T,N}) where {T<:Number,N} =
    NIVolume{typeof(one(T)*1f0+1f0),N,typeof(raw)}(header, NIfTI1Extension[], raw)

function string_tuple(x::String, n::Int)
    a = codeunits(x)
    padding = zeros(UInt8, n-length(a))
    (a..., padding...)
end
string_tuple(x::AbstractString) = string_tuple(bytestring(x))

# Gets dim to be used in header
function nidim(x::AbstractArray)
    dim = ones(Int16, 8)
    dim[1] = ndims(x)
    dim[2:dim[1]+1] = [size(x)...]
    (dim...,)
end

# Gets datatype to be used in header
function nidatatype(t::Type)
    t = get(NIfTI_DT_BITSTYPES_REVERSE, t, nothing)
    if t == nothing
        error("Unsupported data type $T")
    end
    t
end

# Gets the size of a type in bits
nibitpix(t::Type) = Int16(sizeof(t)*8)

# Constructor
function NIVolume(
    # Optional MRI volume; if not given, an empty volume is used
    raw::AbstractArray{T}=Int16[],
    extensions::Union{Vector{NIfTI1Extension},Nothing}=nothing;

    # Fields specified as UNUSED in NIfTI1 spec
    data_type::AbstractString="", db_name::AbstractString="", extents::Integer=Int32(0),
    session_error::Integer=Int16(0), regular::Integer=Int8(0), glmax::Integer=Int32(0),
    glmin::Integer=Int16(0),

    # The frequency encoding, phase encoding and slice dimensions.
    dim_info::NTuple{3, Integer}=(0, 0, 0),
    # Describes data contained in the file; for valid values, see
    # http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/group__NIfTI1__INTENT__CODES.html
    intent_p1::Real=0f0, intent_p2::Real=0f0, intent_p3::Real=0f0,
    intent_code::Integer=Int16(0), intent_name::AbstractString="",
    # Information about which slices were acquired, in case the volume has been padded
    slice_start::Integer=Int16(0), slice_end::Integer=Int16(0), slice_code::Int8=Int8(0),
    # The size of each voxel and the time step. These are formulated in mm unless xyzt_units is
    # explicitly specified
    voxel_size::NTuple{3, Real}=(1f0, 1f0, 1f0), time_step::Real=0f0, xyzt_units::Int8=Int8(18),
    # Slope and intercept by which volume shoudl be scaled
    scl_slope::Real=1f0, scl_inter::Real=0f0,
    # These describe how data should be scaled when displayed on the screen. They are probably
    # rarely used
    cal_max::Real=0f0, cal_min::Real=0f0,
    # The amount of time taken to acquire a slice
    slice_duration::Real=0f0,
    # Indicates a non-zero start point for time axis
    toffset::Real=0f0,

    # "Any text you like"
    descrip::AbstractString="",
    # Name of auxiliary file
    aux_file::AbstractString="",

    # Parameters for Method 2. See the NIfTI spec
    qfac::Float32=0f0, quatern_b::Real=0f0, quatern_c::Real=0f0, quatern_d::Real=0f0,
    qoffset_x::Real=0f0, qoffset_y::Real=0f0, qoffset_z::Real=0f0,
    # Orientation matrix for Method 3
    orientation::Union{Matrix{Float32},Nothing}=nothing) where {T <: Number}
    local t
    if isempty(raw)
        raw = Int16[]
        t = Int16
    else
        t = T
    end

    if extensions == nothing
        extensions = NIfTI1Extension[]
    end

    method2 = qfac != 0 || quatern_b != 0 || quatern_c != 0 || quatern_d != 0 || qoffset_x != 0 ||
        qoffset_y != 0 || qoffset_z != 0
    method3 = orientation != nothing

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
        slice_end = size(raw, dim_info[3]) - 1
    end

    NIVolume(NIfTI1Header(SIZEOF_HDR, string_tuple(data_type, 10), string_tuple(db_name, 18), extents, session_error,
        regular, to_dim_info(dim_info), nidim(raw), intent_p1, intent_p2,
        intent_p3, intent_code, nidatatype(t), nibitpix(t),
        slice_start, (qfac, voxel_size..., time_step, 0, 0, 0), 352,
        scl_slope, scl_inter, slice_end, slice_code,
        xyzt_units, cal_max, cal_min, slice_duration,
        toffset, glmax, glmin, string_tuple(descrip, 80), string_tuple(aux_file, 24), (method2 || method3),
        method3, quatern_b, quatern_c, quatern_d,
        qoffset_x, qoffset_y, qoffset_z, (orientation[1, :]...,),
        (orientation[2, :]...,), (orientation[3, :]...,), string_tuple(intent_name, 16), NP1_MAGIC), extensions, raw)
end

# Validates the header of a volume and updates it to match the volume's contents
function niupdate(vol::NIVolume{T}) where {T}
    vol.header.sizeof_hdr = SIZEOF_HDR
    vol.header.dim = nidim(vol.raw)
    vol.header.datatype = nidatatype(T)
    vol.header.bitpix = nibitpix(T)
    vol.header.vox_offset = isempty(vol.extensions) ? Int32(352) :
        Int32(mapreduce(esize, +, vol.extensions) + SIZEOF_HDR)
    vol
end
