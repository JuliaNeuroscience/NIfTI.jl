
@inline function to_eltype(x::Int16)
    if x == 1
        return Bool
    elseif x == 2
        return UInt8
    elseif x == 4
        return Int16
    elseif x == 8
        return Int32
    elseif x == 16
        return Float32
    elseif x == 32
        return ComplexF32
    elseif x == 64
        return Float64
    elseif x == 256
        return Int8
    elseif x == 12
        return UInt16
    elseif x == 768
        return UInt32
    elseif x == 1024
        return Int64
    elseif x == 1280
        return UInt64
    elseif x == 1792
        return ComplexF64
    else
        error("Unsupported data type $T")
    end
end

function eltype_to_int(T::DataType)
    if T <: Bool
        return 1
    elseif T <: UInt8
        return 2
    elseif T <: Int16
        return 4
    elseif T <: Int32
        return 8
    elseif T <: Float32
        return 16
    elseif T <: ComplexF32
        return 32
    elseif T <: Float64
        return 64
    elseif T <: Int8
        return 256
    elseif T <: UInt16
        return 512
    elseif T <: UInt32
        return 768
    elseif T <: Int64
        return 1024
    elseif T <: UInt64
        return 1280
    elseif T <: ComplexF64
        return 1792
    end
end

function hdr_to_img(file::AbstractString)
    volume_name = replace(file, r"\.\w+(\.gz)?$" => "")*".img"
    if isfile(volume_name)
        return volume_name
    else
        volume_name *= ".gz"
        if isfile(volume_name)
            return volume_name
        else
            error("NIfTI file is dual file storage, but $volume_name does not exist")
        end
    end
end

function string_tuple(x::String, n::Int)
    a = codeunits(x)
    padding = zeros(UInt8, n-length(a))
    (a..., padding...)
end
string_tuple(x::AbstractString) = string_tuple(bytestring(x))

# Conversion factors to mm/ms
# http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/xyzt_units.html
const SPATIAL_UNIT_MULTIPLIERS = [
    1000,   # 1 => NIfTI_UNITS_METER
    1,      # 2 => NIfTI_UNITS_MM
    0.001   # 3 => NIfTI_UNITS_MICRON
]
const TIME_UNIT_MULTIPLIERS = [
    1000,   # NIfTI_UNITS_SEC
    1,      # NIfTI_UNITS_MSEC
    0.001,  # NIfTI_UNITS_USEC
    1,      # NIfTI_UNITS_HZ
    1,      # NIfTI_UNITS_PPM
    1       # NIfTI_UNITS_RADS
]

# Gets dim to be used in header
function to_dim_i16(x::NTuple{N}) where {N}
    return ntuple(Val(8)) do i
        if i === 1
            Int16(N)
        elseif i > (N + 1)
            zero(Int16)
        else
            @inbounds(getfield(x, i - 1))
        end
    end
end

#=
to_spatial_units(x::Int8) = _to_spatial_units(x & 0x07)
@inline function _to_spatial_units(x::UInt8)
    if x === UInt8(1)
        return m  # meter
    elseif x === UInt8(2)
        return mm  # millimeters
    elseif x === UInt8(3)
        return μm  # microns
    else
        return Unitful.FreeUnits{(),NoDims,nothing}()
    end
end

to_time_units(x::Int8) = _to_time_units(x & 0x38)
@inline function _to_time_units(x::Int8)
    if x === Int8(8)
        return s  # second
    elseif x === Int8(16)
        return ms  # milliseconds
    elseif x === Int8(24)
        return μs # microseconds
    elseif x === Int8(32)
        return Hz  # hertz
    #(UInt8(40), u"ppm"),  # parts per million
    elseif x === Int8(48)
        return rad/s  # radians per second
    else
        return Unitful.FreeUnits{(),NoDims,nothing}()
    end
end

units_to_int8(::typeof(Unitful.s)) = Int8(1)
units_to_int8(::typeof(Unitful.ms)) = Int8(16)
units_to_int8(::typeof(Unitful.μs)) = Int8(24)
units_to_int8(::typeof(Unitful.Hz)) = Int8(32)
units_to_int8(::typeof((Unitful.rad/Unitful.s))) = Int8(48)
units_to_int8(::Unitful.FreeUnits{(), NoDims, nothing}) = Int8(0)

function space_time_to_xyzt(ss, tt)
    return (units_to_int8(ss) & 0x07) | (units_to_int8(tt) & 0x38) 
end
=#

