# DISCLAIMER!!! Everything in this file should be considered pre-alpha and
# subject to change.
#
# * This provides a tentative API for moving back and forth between a channelview
#   like space for stats ands storing the structs in a single dimension.
# * The reall work here was making sure that the appropriate parameters used in
#   the NIfTI standard matched up with those in Distribution. For the most part
#   they did, but there are some things that will need some creative handling
#   (such as p-values, z-scores, etc.)

"""
    CorrelationCoefficientIntent

    p1 = degrees of freedom
       R/sqrt(1-R*R) is t-distributed with p1 DOF. */
"""
struct CorrelationCoefficient end
numerictag(::Type{CorrelationCoefficient}) = 2
stringtag(::Type{CorrelationCoefficient}) = "NIFTI_INTENT_CORREL"


#=
    fielddims(::Type{TDist}, dimnames::Tuple)
=#
function fielddims(::Type{I}, dimnames::Tuple) where  {I<:TDist}
    dof1 = dim(dimnames, :ν1)
    if dof1 === 0
        dof1 = dim(dimnames, :dof1)
    end
    if dof1 === 0
        dof1 = dim(dimnames, :dof)
    end

    return dof1
end
numerictag(::Type{<:TDist}) = 3
stringtag(::Type{<:TDist}) = "NIFTI_INTENT_TTEST"



#=
    FScoreIntent

*! [C2, chap 27] Fisher F statistic (2 params):
       p1 = numerator DOF, p2 = denominator DOF. */
=#
function fielddims(::Type{I}, dimnames::Tuple) where  {I<:FDist}
    dof1 = dim(dimnames, :ν1)
    if dof1 === 0
        dof1 = dim(dimnames, :dof1)
    end
    if dof1 === 0
        dof1 = dim(dimnames, :dof)
    end

    dof2 = dim(dimnames, :ν2)
    if dof2 === 0
        dof2 = dim(dimnames, :dof2)
    end

    return dof1, dof2
end
numerictag(::Type{<:FDist}) = 4
stringtag(::Type{<:FDist}) = "NIFTI_INTENT_FTEST"

# TODO ZScore
struct ZScore end
numerictag(::Type{<:ZScore}) = 5
stringtag(::Type{<:ZScore}) = "NIFTI_INTENT_ZSCORE"

numerictag(::Type{<:Chisq}) = 6
stringtag(::Type{<:Chisq}) = "NIFTI_INTENT_CHISQ"

"""

 /*! [C2, chap 25] Beta distribution (2 params): p1=a, p2=b.
      Density(x) proportional to x^(a-1) * (1-x)^(b-1). */
"""
function fielddims(::Type{I}, dimnames::Tuple) where {I<:Beta}
    indb = dim(dimnames, :b)
    if ind === 0
        inda = dim(dimnames, :β)
    end
    inda = dim(dimnames, :a)
    if ind === 0
        inda = dim(dimnames, :α)
    end
    return inda, indb
end
numerictag(::Type{<:Beta}) = 7
stringtag(::Type{<:Beta}) = "NIFTI_INTENT_BETA"

"""
    fielddims(::Type{Binomial}, dimnames::Tuple)
"""
function fielddims(::Type{I}, dimnames::Tuple) where  {I<:Binomial}
    indn = dim(dimnames, :n)
    if indn === 0
        indn = dim(dimnames, :trials)
    end
    indp = dim(dimnames, :p)
    if indp === 0
        indp = dim(dimnames, :probability)
    end

    return indn, indp
end
numerictag(::Type{<:Binomial}) = 8
stringtag(::Type{<:Binomial}) = "NIFTI_INTENT_BINOM"


function fielddims(::Type{I}, dimnames::Tuple) where {I<:Gamma}
    indshape = dim(dimnames, :α)
    if indshape === 0
        indshape = dim(dimnames, :shape)
    end
    indscale = dim(dimnames, :θ)
    if indscale === 0
        indscale = dim(dimnames, :scale)
    end
    return indshape, indscale
end
numerictag(::Type{<:Gamma}) = 9
stringtag(::Type{<:Gamma}) = "NIFTI_INTENT_GAMMA"

function fielddims(::Type{I}, dimnames::Tuple) where {I<:Poisson}
    ind = dim(dimnames, :λ)
    if ind === 0
        ind = dim(dimnames, :mean)
    end
    return ind
end
numerictag(::Type{<:Poisson}) = 10
stringtag(::Type{<:Poisson}) = "NIFTI_INTENT_POISSON"


function fielddims(::Type{I}, dimnames::Tuple) where {I<:Normal}
    indmean = dim(dimnames, :μ)
    if indmean === 0
        indmean = dim(dimnames, :mean)
    end
    indsd = dim(dimnames, :σ)
    if indsd === 0
        indsd = dim(dimnames, :sd)
    end
    if indsd === 0
        indsd = dim(dimnames, :standard_deviation)
    end
 
    return indmean, indsd
end
numerictag(::Type{<:Normal}) = 11
stringtag(::Type{<:Normal}) = "NIFTI_INTENT_NORMAL"

function fielddims(::Type{I}, dimnames::Tuple) where  {I<:NoncentralF}
    dof1 = dim(dimnames, :ν1)
    if dof1 === 0
        dof1 = dim(dimnames, :dof1)
    end
    if dof1 === 0
        dof1 = dim(dimnames, :dof)
    end

    dof2 = dim(dimnames, :ν2)
    if dof2 === 0
        dof2 = dim(dimnames, :dof2)
    end

    noncentr = dim(dimnames, :λ)
    if noncentr === 0
        noncentr = dim(dimnames, :noncentrality)
    end
    return dof1, dof2, noncentrality
end
numerictag(::Type{<:NoncentralF}) = 12
stringtag(::Type{<:NoncentralF}) = "NIFTI_INTENT_FTEST_NONC"

function fielddims(::Type{I}, dimnames::Tuple) where  {I<:NoncentralChisq}
    dof1 = dim(dimnames, :ν)
    if dof1 === 0
        dof1 = dim(dimnames, :ν1)
    end
    if dof1 === 0
        dof1 = dim(dimnames, :dof1)
    end
    if dof1 === 0
        dof1 = dim(dimnames, :dof)
    end
    noncentr = dim(dimnames, :λ)
    if noncentr === 0
        noncentr = dim(dimnames, :noncentrality)
    end

    return dof1, noncentrality
end
numerictag(::Type{<:NoncentralChisq}) = 13
stringtag(::Type{<:NoncentralChisq}) = "NIFTI_INTENT_CHISQ_NONC"

"""
    Logistic
"""
function fielddims(::Type{I}, dimnames::Tuple) where  {I<:Logistic}
    indloc = dim(dimnames, :μ)
    if indloc === 0
        indloc = dim(dimnames, :location)
    end
    indscale = dim(dimnames, :θ)
    if indscale === 0
        indscale = dim(dimnames, :scale)
    end

    return indloc, indscale
end
numerictag(::Type{<:Logistic}) = 14
stringtag(::Type{<:Logistic}) = "NIFTI_INTENT_LOGISTIC"

"""
Laplace
"""
function fielddims(::Type{I}, dimnames::Tuple) where  {I<:Laplace}
    indloc = dim(dimnames, :μ)
    if indloc === 0
        indloc = dim(dimnames, :location)
    end
    indscale = dim(dimnames, :β)
    if indscale === 0
        indscale = dim(dimnames, :scale)
    end

    return indloc, indscale
end
numerictag(::Type{<:Laplace}) = 15
stringtag(::Type{<:Laplace}) = "NIFTI_INTENT_LAPLACE"

"""
    Uniform
"""
function fielddims(::Type{I}, dimnames::Tuple) where  {I<:Uniform}
    indlo = dim(dimnames, :a)
    if indlo === 0
        indlo = dim(dimnames, :low)
    end
    indhi = dim(dimnames, :b)
    if indhi === 0
        indhi = dim(dimnames, :high)
    end

    return indlo, indhi
end
numerictag(::Type{<:Uniform}) = 16
stringtag(::Type{<:Uniform}) = "NIFTI_INTENT_UNIFORM"

"""
    NoncentralT
"""
function fielddims(::Type{I}, dimnames::Tuple) where  {I<:NoncentralT}
    dof1 = dim(dimnames, :ν)
    if dof1 === 0
        dof1 = dim(dimnames, :ν1)
    end
    if dof1 === 0
        dof1 = dim(dimnames, :dof1)
    end
    if dof1 === 0
        dof1 = dim(dimnames, :dof)
    end
    noncentr = dim(dimnames, :λ)
    if noncentr === 0
        noncentr = dim(dimnames, :noncentrality)
    end

    return dof1, noncentrality
end
numerictag(::Type{<:NoncentralT}) = 17
stringtag(::Type{<:NoncentralT}) = "NIFTI_INTENT_TTEST_NONC"

"""
    Weibull
"""
function fielddims(::Type{I}, dimnames::Tuple) where {I<:Weibull}
    indshape = dim(dimnames, :α)
    if indshape === 0
        indshape = dim(dimnames, :shape)
    end
    indscale = dim(dimnames, :θ)
    if indscale === 0
        indscale = dim(dimnames, :scale)
    end
    return indshape, indscale
end
numerictag(::Type{<:Weibull}) = 18
stringtag(::Type{<:Weibull}) = "NIFTI_INTENT_WEIBULL"


"""
    Chi
"""
function fielddims(::Type{I}, dimnames::Tuple) where {I<:Chi}
    ind = dim(dimnames, :ν)
    if ind === 0
        ind = dim(dimnames, :dof)
    end
    if ind === 0
        ind = dim(dimnames, :dof1)
    end
    return ind
end
numerictag(::Type{<:Chi}) = 19
stringtag(::Type{<:Chi}) = "NIFTI_INTENT_CHI"

"""
    InverseGaussian
"""
function fielddims(::Type{I}, dimnames::Tuple) where {I<:InverseGaussian}
    indmean = dim(dimnames, :μ)
    if indmean === 0
        indmean = dim(dimnames, :mean)
    end
    indshape = dim(dimnames, :λ)
    if indshape === 0
        indshape = dim(dimnames, :shape)
    end
    return indmean, indshape
end
numerictag(::Type{<:InverseGaussian}) = 20
stringtag(::Type{<:InverseGaussian}) = "NIFTI_INTENT_INVGAUSS"

struct PValue end
numerictag(::Type{PValue}) = 22
stringtag(::Type{PValue}) = "NIFTI_INTENT_PVAL"


struct LogPValue end
numerictag(::Type{LogPValue}) = 23
stringtag(::Type{LogPValue}) = "NIFTI_INTENT_LOGPVAL"

struct Log10PValue end
numerictag(::Type{Log10PValue}) = 24
stringtag(::Type{Log10PValue}) = "NIFTI_INTENT_LOG10PVAL"
