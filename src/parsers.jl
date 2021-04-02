

const SIZEOF_HDR1 = Int32(348)
const SIZEOF_HDR2 = Int32(540)


const NP1_MAGIC = (0x6e,0x2b,0x31,0x00)
const NI1_MAGIC = (0x6e,0x69,0x31,0x00)
const NP2_MAGIC = (0x6e,0x2b,0x32,0x00,0x0d,0x0a,0x1a,0x0a)
const NI2_MAGIC = (0x6e,0x69,0x32,0x00,0x0d,0x0a,0x1a,0x0a)


function to_eltype_error(@nospecialize(x))
    throw("NIfTI parse erorr: $x does not have a corresponding eltype ")
end

eltype_to_int16(::Type{T}) where {T} = error("Unsupported data type $T")

for T in (Int16, Int32)
    symint = Symbol(lowercase(string(T)))
    @eval begin
        #= TODO handle intent codes
        @inline function to_intent(x::Int16)
            if i === $T(1001)
                return Estimate
            elseif i === $T(1002)
                return NeuroLabel
            elseif i === $T(1003)
                return NeuroNames
            elseif i === $T(1004)
                return SMatrix
            elseif i === $T(1005)
                return Symmetric
            elseif i === $T(1006)
                return DisplacementVector
            elseif i === $T(1007)
                return SVector
            elseif i === $T(1008)
                return Point
            elseif i === $T(1009)
                return Triangle
            elseif i === $T(1010)
                return Quat
            elseif i === $T(1011)
                return Dimensionless
            elseif i === $T(2001)
                return TimeSeries
            elseif i === $T(2002)
                return Vector{Point}
            elseif i === $T(2003)
                return GiftiRGB
            elseif i === $T(2004)
                return GiftiRGBA
            elseif i === $T(2005)
                return Polygon
            elseif i === $T(2006)
               return FSLDisplacementVector # TODO
            elseif i === $T(2007)
                return FSLCubicSplineCoefficient  # TODO
            elseif i === $T(2008)
                return FSLDCTCoefficients  # TODO
            elseif i === $T(2009)
                return FSLQuadraticSplineCoefficients  # TODO
            elseif i === $T(2016)
                return FSLTopupCubicSplineCoefficients  # TODO
            elseif i === $T(2017)
                return TopupQuadraticSplineCoefficients # TODO
            elseif i === $T(2018)
                return TopupField  # TODO
            elseif i == $T(3000)
                return ConnUnkown
            elseif i == $T(3001)
                return ConnDense
            elseif i == $T(3002)
                return ConnDenseSeries
            elseif i == $T(3003)
                return ConnParcels
            elseif i == $T(3004)
                return ConnParcelSeries
            elseif i == $T(3006)
                return ConnDenseScalar
            elseif i == $T(3007)
                return ConnDenseLabel
            elseif i == $T(3008)
                return ConnParcelScalr
            elseif i == $T(3009)
                return ConnParcelDense
            elseif i == $T(3010)
                return ConnDenseParcel
            elseif i == $T(3013)
                return ConnPPSc
            else
                return UnkownIntent
            end
        end
        =#

        @inline function $(Symbol(:xform_to_, symint))(x::Symbol)
            if x === :ScannerSpace
                return $T(1)
            elseif x === :AnatomicalSpace
                return $T(2)
            elseif x === :TalairachSpace
                return $T(3)
            elseif x === :MNI152Space
                return $T(4)
            elseif x === :OtherTemplate
                return $T(5)
            else
                return $T(0)
            end
        end

        @inline function to_xform(i::$T)
            if i === $T(1)
                return :ScannerSpace
            elseif i === $T(2)
                return :AnatomicalSpace
            elseif i === $T(3)
                return :TalairachSpace
            elseif i === $T(4)
                return :MNI152Space
            elseif i === $T(5)
                return :OtherTemplate
            else
                return :UnkownSpace
            end
        end

        @inline function to_eltype(x::$T)
            if x === $T(1)
                return Bool
            elseif x === $T(2)
                return UInt8
            elseif x === $T(4)
                return Int16
            elseif x === $T(8)
                return Int32
            elseif x === $T(16)
                return Float32
            elseif x === $T(32)
                return ComplexF32
            elseif x === $T(64)
                return Float64
            elseif x === $T(256)
                return Int8
            elseif x === $T(512)
                return UInt16
            elseif x === $T(768)
                return UInt32
            elseif x === $T(1024)
                return Int64
            elseif x === $T(1280)
                return UInt64
            elseif x === $T(1792)
                return ComplexF64
            else
                to_eltype_error(x)
            end
        end

        $(Symbol(:eltype_to_, symint))(::Type{Bool}) = $T(1)
        $(Symbol(:eltype_to_, symint))(::Type{UInt8}) = $T(2)
        $(Symbol(:eltype_to_, symint))(::Type{Int16}) = $T(4)
        $(Symbol(:eltype_to_, symint))(::Type{Int32}) = $T(8)
        $(Symbol(:eltype_to_, symint))(::Type{Float32}) = $T(16)
        $(Symbol(:eltype_to_, symint))(::Type{ComplexF32}) = $T(32)
        $(Symbol(:eltype_to_, symint))(::Type{Float64}) = $T(64)
        $(Symbol(:eltype_to_, symint))(::Type{Int8}) = $T(256)
        $(Symbol(:eltype_to_, symint))(::Type{UInt16}) = $T(512)
        $(Symbol(:eltype_to_, symint))(::Type{UInt32}) = $T(768)
        $(Symbol(:eltype_to_, symint))(::Type{Int64}) = $T(1024)
        $(Symbol(:eltype_to_, symint))(::Type{UInt64}) = $T(1280)
        $(Symbol(:eltype_to_, symint))(::Type{ComplexF64}) = $T(1792)
    end
end


to_dimensions(d::NTuple{8,Int16}) = ntuple(i ->Int(@inbounds(getfield(d, i + 1))), first(d))
to_dimensions(d::NTuple{8,Int}) = ntuple(i ->@inbounds(getfield(d, i + 1)), first(d))

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

