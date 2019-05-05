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

struct NoIntent end

const NiftiIntents = Dict{Int, Type}(
    0    => NoIntent,
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

intent(img::ImageMeta{T,N,A,ImageProperties{:NII}}) where {T,N,A} = intent(properties(img))
intent(s::ImageStream) = intent(properties(s))

intent(p::ImageProperties) = intent(header(p))
intent(p::ImageProperties{:header}) = @get p "intent" NoIntent
intent(::Nothing) = NoIntent

intent(A::AbstractArray{T}) where {T<:BitTypes} = NoIntent

intent(A::AbstractArray{<:RGB}) = NoIntent
intent(A::AbstractArray{<:RGBA}) = NoIntent

intent(A::AbstractArray{<:Quat}) = Quat
intent(A::AbstractArray{<:SVector}) = SVector
intent(A::AbstractArray{<:SMatrix}) = SMatrix
intent(A::AbstractArray{T}) where {T<:Union{StatP0,StatP1,StatP2,StatP3}} = T
intent(A::AbstractArray{<:Point}) = Point
intent(A::AbstractArray{<:Triangle}) = Triangle
intent(A::AbstractVector{<:RGB}) = GiftiRGB
intent(A::AbstractVector{<:RGBA}) = GiftiRGBA

# TODO
# CustomLabels
# NeuroNames
#intent(A::IndirectArray) = CustomLabels

# TODO
# SymmetricMatrix
# Dimensionless
# GiftiTimeSeries
# GiftiShape)

"""
    intentname(img)
"""
intentname(img::ImageFormat{format"NII"}) = img["header"]["intentname"]


intentname(img::ImageMeta{T,N,A,ImageProperties{:NII}}) where {T,N,A} = intentname(properties(img))
intentname(s::ImageStream) = intentname(properties(s))

intentname(p::ImageProperties) = intentname(header(p))
intentname(p::ImageProperties{:header}) = @get p "intentname" intentname() 
intentname(::Nothing) = intentname()


intentname(A::AbstractArray) = intentname()
intentname() = String(fill(UInt8(0), 16))




"""
    intentparams(img)
"""
intentparams(img::ImageMeta{T,N,A,ImageProperties{:NII}}) where {T,N,A} = intentparams(properties(img))
intentparams(s::ImageStream) = intentparams(properties(s))
intentparams(A::AbstractArray) = ntuple(_->float(0), Val(3))::Tuple{Float64,Float64,Float64}

intentparams(p::ImageProperties) = intentparams(header(p))
intentparams(p::ImageProperties{:header}) = @get p "intentparams" Float64[0, 0, 0]
intentparams(::Nothing) = Float64[0, 0, 0]


# array dimensions are in row major form
intentparams(A::AbstractArray{<:AbstractMatrix}) = (Float64(size(A,2)), Float64(size(A,1)), Float64(0.0))


intentaxis(s::AbstractArray, ::Type) = Axis{:dim5}(range(Float64(1), step=s.pixdim[5], length=s.size[5]))
intentaxis(s::AbstractArray, ::Type{SMatrix}) = [Axis{:col}(range(1, length=intentparams(s)[2])), Axis{:row}(range(1, length=intentparams(s)[2]))]
intentaxis(s::AbstractArray, ::Type{Triangle}) = Axis{:triangledim}(1:3)
intentaxis(s::AbstractArray, ::Type{Point}) = Axis{:pointdim}(1:3)
intentaxis(s::AbstractArray, ::Type{GiftiNodeIndex}) = Axis{:pointdim}(1:3)
intentaxis(s::AbstractArray, ::Type{Quat}) = Axis{:quatdim}(1:4)
intentaxis(s::AbstractArray, ::Type{DisplacementVector}) = Axis{:vecdim}(range(1, length=s.size[5]))
intentaxis(s::AbstractArray, ::Type{SVector}) = Axis{:vecdim}(range(1, length=s.size[5]))
intentaxis(s::AbstractArray, ::Type{GiftiRGB}) = Axis{:colordim}(1:3)
intentaxis(s::AbstractArray, ::Type{GiftiRGBA}) = Axis{:colordim}(1:4)

function niaxes(sz::NTuple{N,Int}, xyzt_units::Integer, toffset::F,
                pixdim::NTuple{N}, s::IOMeta) where {N,F}
   # determine axisnames
    if s["header"]["sformcode"] != :Unkown
       ori = orientation(s["header"]["sform"])
    else
        ori = orientation(s["header"]["qform"])
    end

    su = get(NiftiUnits, xyzt_units & 0x07, 1)
    axs = Vector{Axis}(undef, min(4,N))
    axs[1] = Axis{ori[1]}(range(F(1), step=pixdim[1], length=sz[1])*su)
    if N > 1
        if pixdim[3] == 1
            axs[2] = Axis{ori[2]}(Base.OneTo(sz[2])*su)
        else
            axs[2] = Axis{ori[2]}(range(F(1), step=pixdim[2], length=sz[2])*su)
        end
    end
    if N > 2
        if pixdim[3] == 1
            axs[3] = Axis{ori[3]}(Base.OneTo(sz[3])*su)
        else
            axs[3] = Axis{ori[3]}(range(F(1), step=pixdim[3], length=sz[3])*su)
        end
    end
    if N > 3
        tu = get(NiftiUnits, xyzt_units & 0x38, 1)
        axs[4] = Axis{:time}(range(toffset, step=pixdim[4], length=sz[4])*tu)
    end
    if N > 4
        append!(axs, intentaxis(s, pixdim, s["header"]["intent"]))
    end
    if N > 5
        for i in 6:N
            append!(axs, Axis{Symbol("dim$i")}(range(F(1), step=pixdim[i], length=sz[i])))
        end
    end
    return (axs...,)
end

intentaxis(s, pixdim, ::Type) = Axis{:dim5}(range(Float64(1), step=s.pixdim[5], length=s.size[5]))
intentaxis(s, pixdim, ::Type{SMatrix}) = [Axis{:col}(Base.OneTo(intentparams(s)[2])), Axis{:row}(Base.OneTo(intentparams(s)[1]))]
intentaxis(s, pixdim, ::Type{Triangle}) = Axis{:triangledim}(1:3)
intentaxis(s, pixdim, ::Type{Point}) = Axis{:pointdim}(1:3)
intentaxis(s, pixdim, ::Type{GiftiNodeIndex}) = Axis{:pointdim}(1:3)
intentaxis(s, pixdim, ::Type{Quat}) = Axis{:quatdim}(1:4)
intentaxis(s, pixdim, ::Type{DisplacementVector}) = Axis{:vecdim}(range(1, length=s.size[5]))
intentaxis(s, pixdim, ::Type{SVector}) = Axis{:vecdim}(range(1, length=s.size[5]))
intentaxis(s, pixdim, ::Type{GiftiRGB}) = Axis{:colordim}(1:3)
intentaxis(s, pixdim, ::Type{GiftiRGBA}) = Axis{:colordim}(1:4)


