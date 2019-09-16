
function num2intent(i::Integer)
    if i < 1
        return UnknownIntent
    elseif i < 1000
        statsintent(i)
    elseif i < 2000
        structintent(i)
    elseif i < 3000
        giftiintent(i)
    elseif i < 4000
        ciftiintent(i)
    else
        return UnknownIntent
    end
end

#intent(img::ImageMeta{T,N,A,ImageProperties{:NII}}) where {T,N,A} = intent(properties(img))

#intent(p::ImageProperties) = intent(header(p))
#intent(p::ImageProperties{:header}) = @get p "intent" NoIntent
intent(::Nothing) = NoIntent

intent(A::AbstractArray{T}) where {T<:BitTypes} = NoIntent

intent(A::AbstractArray{<:RGB}) = NoIntent
intent(A::AbstractArray{<:RGBA}) = NoIntent

intent(A::AbstractArray{<:Quat}) = Quat
intent(A::AbstractArray{<:SVector}) = SVector
intent(A::AbstractArray{<:SMatrix}) = SMatrix
intent(A::AbstractArray{T}) where {T<:Distribution} = T
intent(A::AbstractArray{<:Point}) = Point
intent(A::AbstractArray{<:Triangle}) = Triangle
intent(A::AbstractVector{<:RGB}) = GiftiRGB
intent(A::AbstractVector{<:RGBA}) = GiftiRGBA

# TODO
# CustomLabels
# NeuroNames
#intent(A::IndirectArray) = CustomLabels

# TODO
# Dimensionless
# GiftiTimeSeries
# GiftiShape)

# array dimensions are in row major form



intentaxis(s, pixdim, ::Type) = Axis{:dim5}(range(Float64(1), step=pixdim[5], length=size[5]))

function guessintent(ext::AbstractVector{<:AbstractString}, ::Type{A}) where {A<:AbstractArray}
    if ext[end] == "gz"  # CIfTI and GIfTI aren't supposed to be zipped
        return A
    elseif ext[end] == "nii"
        if ext[end-1] == "dconn"
            return ConnDense
        elseif ext[end-1] == "dtseries"
            return ConnDenseSeries
        elseif ext[end-1] == "pconn"
            return ConnParcels
        elseif ext[end-1] == "ptseries"
            return ConnParcelSries
        elseif ext[end-1] == "dscalar"
            return ConnDenseScalar
        elseif ext[end-1] == "dlabel"
            return ConnDenseLabel
        elseif ext[end-1] == "pscalar"
            return ConnParcelScalr
        elseif ext[end-1] == "pdconn"
            return ConnParcelDense
        elseif ext[end-1] == "dpconn"
            return ConnDenseParcel
        elseif ext[end-1] == "pconnscalar"
            return ConnPPSc
        elseif ext[end-1] == "dfan"
            return FiberFan
        elseif ext[end-1] == "dfansamp"
            return FanSamples
        else
            return A  # assume no meaningful extra extension was present
        end
    elseif ext[end] == "gii"
            if ext[end-1] == "coord"
            return GiftiCoordinate
        elseif ext[end-1] == "func"
            return GiftiFunctional
        elseif ext[end-1] == "label"
            return GiftiLabel
        elseif ext[end-1] == "rgba"
            return GiftiRGBA
        elseif ext[end-1] == "shape"
            return GiftiShape
        elseif ext[end-1] == "surf"
            return GiftiSurface
        elseif ext[end-1] == "tensor"
            return GiftiTensor
        elseif ext[end-1] == "time"
            return GiftiTime
        elseif ext[end-1] == "topo"
            return GiftiTopology
        elseif ext[end-1] == "vector"
            return GiftiVector
        elseif ext[end] == "gii"
            return GIfTI
        else
            return A
        end
    else
        return A
    end
end

const STRING_INTENTS = Dict(
    "NIFTI_INTENT_POINTSET" => Point,
    "NIFTI_INTENT_TRIANGLE" => Triangle,
    "NIFTI_INTENT_NODE_INDEX" => NodeIndex,
    "NIFTI_INTENT_NORMAL" => Normal,
    "NIFTI_INTENT_NONE" => UnknownIntent,
    "NIFTI_INTENT_SHAPE" => Polygon,
    "NIFTI_INTENT_LABEL" => NeuroLabel,
    "NIFTI_INTENT_TIME_SERIES" => "timeseries",
    "NIFTI_INTENT_RGB_VECTOR" => GiftiRGB,
    "NIFTI_INTENT_RGBA_VECTOR" => GiftiRGBA,
    "NIFTI_INTENT_GENMATRIX" => "genmatrix",

    "NIFTI_INTENT_CORREL" => CorrelationCoefficient,
    "NIFTI_INTENT_TTEST" => TDist,
    "NIFTI_INTENT_FTEST" => FDist,

    "NIFTI_INTENT_ZSCORE" => ZScore,
    "NIFTI_INTENT_CHISQ" => Chisq,
    "NIFTI_INTENT_BETA" => Beta,
    "NIFTI_INTENT_BINOM" => Binomial,
    "NIFTI_INTENT_GAMMA" => Gamma,
    "NIFTI_INTENT_POISSON" => Poisson,
    "NIFTI_INTENT_FTEST_NONC" => NoncentralF,
    "NIFTI_INTENT_CHISQ_NONC" => NoncentralChisq,
    "NIFTI_INTENT_LOGISTIC" => Logistic,
    "NIFTI_INTENT_LAPLACE" => Laplace,
    "NIFTI_INTENT_UNIFORM" => Uniform,
    "NIFTI_INTENT_TTEST_NONC" => NoncentralT,
    "NIFTI_INTENT_WEIBULL" => Weibull,
    "NIFTI_INTENT_CHI" => Chi,
    "NIFTI_INTENT_INVGAUSS" => InverseGaussian,
    "NIFTI_INTENT_EXTVAL" => "extval",
    "NIFTI_INTENT_PVAL" => PValue,
    "NIFTI_INTENT_LOGPVAL" => LogPValue,
    "NIFTI_INTENT_LOG10PVAL" => Log10PValue,
    "NIFTI_INTENT_ESTIMATE" => Estimate,


    "NIFTI_INTENT_NEURONAME" => NeuroName,

    "NIFTI_INTENT_SYMMATRIX" => Symmetric,
    "NIFTI_INTENT_DISPVECT" => DisplacementVector,
    "NIFTI_INTENT_VECTOR" => SVector,

    "NIFTI_INTENT_QUATERNION" => Quat,
    "NIFTI_INTENT_DIMLESS" => Dimensionless
  )


