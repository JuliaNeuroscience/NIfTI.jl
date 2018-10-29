# TODO
struct Correlation <: ContinuousUnivariateDistribution end
struct ZScore <: ContinuousUnivariateDistribution end
struct Pvalue <: ContinuousUnivariateDistribution end
struct logPvalue <: ContinuousUnivariateDistribution end
struct log10Pvalue <: ContinuousUnivariateDistribution end

# what to do about stats that have mixed many params and single params
# (eg., nsub is same but mean voxel is diff)
#
# params(img::AbstractArray{T<:Distribution})

function intent_stats(statistic::String, img::AbstractArray, header::Dict{String,Any})
    if ndims(img) < 5
        if intent_code == "Correlation"
            header["intent"] = Correlation
        elseif intent_code == "TTest"
            header["intent"] = TDist(hdr.intent_p1)
        elseif intent_code == "FTest"
            header["intent"] = FDist(hdr.intent_p1,
                                     hdr.intent_p2)
        elseif intent_code == "ZScore"
            header["intent"] = ZScore()
        elseif intent_code == "Chisq"
            header["intent"] = Chisq(hdr.intent_p1)
        elseif intent_code == "Beta"
            header["intent"] = Beta(hdr.intent_p1,
                                    hdr.intent_p2)
        elseif intent_code == "Binomial"
            header["intent"] = Binomial(hdr.intent_p1,
                                        hdr.intent_p2)
        elseif intent_code == "Gamma"
            header["intent"] = Gamma(hdr.intent_p1,
                                    hdr.intent_p2)
        elseif intent_code == "Poisson"
            header["intent"] = Poisson(hdr.intent_p1)
        elseif intent_code == "Normal"
            header["intent"] = Normal(hdr.intent_p1,
                                     hdr.intent_p2)
        elseif intent_code == "NoncentralFTest"
            header["intent"] = NoncentralF(hdr.intent_p1,
                                           hdr.intent_p2,
                                           hdr.intent_p3)
        elseif intent_code == "NoncentralChisq"
            header["intent"] = NoncentralChisq(hdr.intent_p1,
                                               hdr.intent_p2)
        elseif intent_code == "Logistic"
            header["intent"] = Logistic(hdr.intent_p1,
                                        hdr.intent_p2)
        elseif intent_code == "Laplace"
            header["intent"] = Laplace(hdr.intent_p1,
                                       hdr.intent_p2)
        elseif intent_code == "Uniform"
            header["intent"]  = Uniform(hdr.intent_p1,
                                        hdr.intent_p2)
        elseif intent_code == "NoncentralTTest"
            header["intent"] = NoncentralT(hdr.intent_p1,
                                           hdr.intent_p2)
        elseif intent_code == "Chi"
            header["intent"] = Chi(hdr.intent_p1)
        elseif intent_code == "Weibull"
            header["intent"] = Weibull(hdr.intent_p1,
                                       hdr.intent_p2)
        elseif intent_code == "InverseGaussian"
            header["intent"] = InverseGaussian(hdr.intent_p1,
                                               hdr.intent_p2)
        elseif intent_code == "ExtremeValue"
            # TODO double check this matches NIfTI
            header["intent"] = GeneralizedExtremeValue(hdr.intent_p1,
                                                       hdr.intent_p2,
                                                       hdr.intent_p3)
        elseif intent_code == "Pvalue"
            header["intent"] = Pvalue()
        elseif intent_code == "logPvalue"
            header["intent"] = logPvalue()
        elseif intent_code == "log10Pvalue"
            header["intent"] = log10Pvalue()
        end
    else
        # TODO 5th axis for parameters
        # probably do something like RGB but with Distribution parameters
    end
    ImageMeta(img, properties=Dict{String,Any}(
                                  "filetype" => ftype,
                                  "header" => header))
end

