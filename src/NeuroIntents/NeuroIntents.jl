using Distributions: TDist,
                     Chi,
                     Chisq,
                     Poisson,
                     FDist,
                     Beta,
                     Binomial,
                     Gamma,
                     Normal,
                     NoncentralT,
                     NoncentralChisq,
                     Logistic,
                     Uniform,
                     NoncentralF,
                     GeneralizedExtremeValue,
                     Weibull,
                     Distribution

struct UnknownIntent end

intent(x::Any) = getter(x, "intent", Type, UnknownIntent)

intent!(x::Any, val::Type) = setter!(x, "intent", val, Type)

include("cifti.jl")
include("gifti.jl")
include("stats.jl")
include("structures.jl")
include("intent.jl")
