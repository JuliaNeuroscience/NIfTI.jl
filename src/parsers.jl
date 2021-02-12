

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

@inline function to_ecode(x::Int32)
    if x === Int32(2)
        return :DICOM
    elseif x === Int32(4)
        return :AFNI
    elseif x === Int32(6)
        return :Comment
    elseif x === Int32(8)
        return :XCEDE
    elseif x === Int32(10)
        return :JimDimInfo
    elseif x === Int32(12)
        return :WorkflowFWDS
    elseif x === Int32(14)
        return :Freesurfer
    elseif x === Int32(16)
        return :PyPickly
    elseif x === Int32(18)
        return :MindIdent
    elseif x === Int32(20)
        return :BValue
    elseif x === Int32(22)
        return :SphericalDirection
    elseif x === Int32(24)
        return :DTComponent
    elseif x === Int32(26)
        return :SHCDegreeOrder
    elseif x === Int32(28)
        return :Voxbo
    elseif x === Int32(30)
        return :Caret
    elseif x === Int32(32)
        return :CIfTI
    elseif x === Int32(34)
        return :VariableFrameTiming
    elseif x === Int32(38)
        return :Eval
    elseif x === Int32(40)
        return :Matlab
    else
        return :Ignore
    end
end

@inline function ecode_to_int(x::Symbol)
    if x === :DICOM
        return Int32(2)
    elseif x === :AFNI
        return Int32(4)
    elseif x === :Comment
        return Int32(6)
    elseif x === :XCEDE
        return Int32(8)
    elseif x === :JimDimInfo
        return Int32(10)
    elseif x === :WorkflowFWDS
        return Int32(12)
    elseif x === :Freesurfer
        return Int32(14)
    elseif x === :PyPickly
        return Int32(16)
    elseif x === :MindIdent
        return Int32(18)
    elseif x === :BValue
        return Int32(20)
    elseif x === :SphericalDirection
        return Int32(22)
    elseif x === :DTComponent
        return Int32(24)
    elseif x === :SHCDegreeOrder
        return Int32(26)
    elseif x === :Voxbo
        return Int32(28)
    elseif x === :Caret
        return Int32(30)
    elseif x === :CIfTI
        return Int32(32)
    elseif x === :VariableFrameTiming
        return Int32(34)
    elseif x === :Eval
        return Int32(38)
    elseif x === :Matlab
        return Int32(40)
    else
        return Int32(0)
    end
end

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

