const SIZEOF_HDR1 = Int32(348)
const SIZEOF_HDR2 = Int32(540)

primitive type Float128 <: AbstractFloat 128 end
const ComplexF128 = Complex{Float128}

BitTypes = Union{Integer,AbstractFloat,ComplexF128,ComplexF32}

const NiftiDatatypes = Dict{Int16,Type}([
    (Int16(2), UInt8),
    (Int16(4), Int16),
    (Int16(8), Int32),
    (Int16(16), Float32),
    (Int16(32), ComplexF32),
    (Int16(64), Float64),
    (Int16(128), RGB),  # 3 8 bit bytes
    (Int16(256), Int8),
    (Int16(512), UInt16),
    (Int16(768), UInt32),
    (Int16(1024), Int64),
    (Int16(1280), UInt64),
    (Int16(1536), Float128),
    (Int16(1792), ComplexF64),
    (Int16(2048), ComplexF128),
    (Int16(2304), RGBA)  # 4 8 bit bytes
])

const NiftiDatatypesReverse = Dict{Type,Int16}()
for (k, v) in NiftiDatatypes
    NiftiDatatypesReverse[v] = k
end

const NP1_MAGIC = (0x6e,0x2b,0x31,0x00)
const NI1_MAGIC = (0x6e,0x69,0x31,0x00)
const NP2_MAGIC = (0x6e,0x2b,0x32,0x00, 0x0d, 0x0a,0x1a,0x0a)

#@unit ppm "ppm" PartsPerMillion 1//1000000 false
#Unitful.register(ppm)

const NiftiUnits = Dict([
    (Int16(0), nothing),
    (Int16(1), u"m"),  # meter
    (Int16(2), u"mm"),  # millimeters
    (Int16(3), u"μm"),  # microns
    (Int16(8), u"s"),  # second
    (Int16(16), u"ms"),  # milliseconds
    (Int16(24), u"μs"), # microseconds
    (Int16(32), u"Hz"),  # hertz
    #(Int16(40), u"ppm"),  # parts per million
    (Int16(48), u"rad/s")  # radians per second
])

const NiftiUnitsReverse = Dict()
for (k, v) in NiftiUnits
    NiftiUnitsReverse[v] = k
end

const ANALYZE75_ORIENT = Dict{Symbol,Int16}([
    (:a75_transverse_unflipped, 0),
    (:a75_coronal_unflipped, 1),
    (:a75_sagittal_unflipped, 2),
    (:a75_transverse_flipped, 3),
    (:a75_coronal_flipped, 4),
    (:a75_sagittal_flipped, 5),
    (:a75_orient_unknown, 6)
])

const NIFTI_ORIENTATION = Dict{Int16, Symbol}([
    (Int16(1),  :L2R),
    (Int16(-1), :R2L),
    (Int16(2),  :P2A),
    (Int16(-2), :A2P),
    (Int16(3),  :I2S),
    (Int16(-3), :S2I)
])

ori2space = Dict{Symbol, Symbol}(
    :L2R => :L,
    :R2L => :R,
    :P2A => :P,
    :A2P => :A,
    :I2S => :I,
    :S2I => :S
)
#ori_lut = ((   1,   -1,    2,   -2,    3,   -3),
#           (:L2R, :R2L, :P2A, :A2P, :I2S, :S2I),
#           (  :L,   :R,   :P,   :A,   :I,   :S))

isradview(axnames::NTuple{3,Symbol}) = axnames == (:L, :A, :S)
isneuroview(axnames::NTuple{3,Symbol}) = axnames == (:R, :A, :S)

const NiftiSliceCodes = Dict{Int16,String}([
    (Int16(0), "Unkown"),
    (Int16(1), "Sequential+Increasing"),
    (Int16(2), "Sequential+Decreasing"),
    (Int16(3), "Alternating+Increasing"),
    (Int16(4), "Alternating+Decreasing"),
    (Int16(5), "Alternating+Increasing#2"),
    (Int16(6), "Alternating+Decreasing#2")
   ])

const NiftiSliceCodesReverse = Dict{String,Int16}()
for (k, v) in NiftiSliceCodes
    NiftiSliceCodesReverse[v] = k
end

const NIFTI_FTYPE = Dict{Int16, String}([
    (Int16(0), "Analyze"),
    (Int16(1), "NIfTI-1Single"),
    (Int16(2), "NIfTI-1Dual"),
    (Int16(3), "ASCII"),
    (Int16(4), "NIfTI-2Single"),
    (Int16(5), "NIfTI-2Dual")
])

const NiftiXForm = Dict{Int,Symbol}([
    (Int16(0), :Unkown),
    (Int16(1), :Scanner_anat),
    (Int16(2), :Aligned_anat),
    (Int16(3), :Talairach),
    (Int16(4), :MNI152)
])

const NiftiXFormReverse = Dict{Symbol,Int}()
for (k, v) in NiftiXForm
    NiftiXFormReverse[v] = k
end

const NIFTI_ECODE = Dict{Int16, Symbol}([
    (Int16(0), :ignore),
    (Int16(2), :dicom),
    (Int16(4), :afni),
    (Int16(6), :comment),
    (Int16(8), :xcede),
    (Int16(10), :jimdiminfo),
    (Int16(12), :workflow_fwds),
    (Int16(14), :freesurfer),
    (Int16(16), :pypickly),
    (Int16(18), :mind_ident),
    (Int16(20), :b_value),
    (Int16(22), :spherical_direction),
    (Int16(24), :dt_component),
    (Int16(26), :shc_degreeorder),
    (Int16(28), :voxbo),
    (Int16(30), :caret),
    (Int16(32), :cifti),
    (Int16(34), :variable_frame_timing),
    (Int16(38), :eval),
    (Int16(40), :matlab)
])


struct Correlation <: ContinuousUnivariateDistribution end
struct ZScore <: ContinuousUnivariateDistribution end
struct PValue <: ContinuousUnivariateDistribution end
struct LogPValue <: ContinuousUnivariateDistribution end
struct Log10PValue <: ContinuousUnivariateDistribution end
struct ParametricEstimate <: ContinuousUnivariateDistribution
    p1::Any
    p2::Any
    p3::Any
end

const StatP0 = Union{Correlation,ZScore,PValue,LogPValue,Log10PValue}
const StatP1 = Union{TDist,Chi,Chisq,Poisson}
const StatP2 = Union{FDist,Beta,Binomial,Gamma,Normal,NoncentralT,NoncentralChisq,Logistic,Uniform}
const StatP3 = Union{NoncentralF,GeneralizedExtremeValue,ParametricEstimate}

struct CustomLabels end
struct NeuroNames end

struct SymmetricMatrix end
struct DisplacementVector end
struct Dimensionless end
struct GiftiTimeSeries end
struct GiftiNodeIndex end
struct GiftiShape end
struct GiftiRGB end
struct GiftiRGBA end

const NiftiIntents = Dict{Int, Type}(
    # spm intent
    2    => Correlation,
    3    => TDist,
    4    => FDist,
    5    => ZScore,
    6    => Chisq,
    7    => Beta,
    8    => Binomial,
    9    => Gamma,
    10   => Poisson,
    11   => Normal,
    12   => NoncentralF,
    13   => NoncentralChisq,
    14   => Logistic,
    15   => Laplace,
    16   => Uniform,
    17   => NoncentralT,
    18   => Weibull,
    19   => Chi,
    20   => InverseGaussian,
    21   => GeneralizedExtremeValue,
    22   => PValue,
    23   => LogPValue,
    24   => Log10PValue,
    1001 => ParametricEstimate,
    1002 => CustomLabels,
    1003 => NeuroNames,
    1004 => SMatrix,
    1005 => SymmetricMatrix,
    1006 => DisplacementVector,
    1007 => SVector,
    1008 => Point,
    1009 => Triangle,
    1010 => Quat,
    1011 => Dimensionless,
    2001 => GiftiTimeSeries,
    2002 => GiftiNodeIndex,
    2003 => GiftiRGB,
    2004 => GiftiRGBA,
    2005 => GiftiShape)

const NiftiIntentsReverse = Dict{Type,Int16}()
for (k, v) in NiftiIntents
    NiftiIntentsReverse[v] = k
end

