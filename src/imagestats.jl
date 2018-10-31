# TODO
struct Correlation <: ContinuousUnivariateDistribution end
struct ZScore <: ContinuousUnivariateDistribution end
struct Pvalue <: ContinuousUnivariateDistribution end
struct logPvalue <: ContinuousUnivariateDistribution end
struct log10Pvalue <: ContinuousUnivariateDistribution end
struct ParametricEstimate <: ContinuousUnivariateDistribution end

# what to do about stats that have mixed many params and single params
# (eg., nsub is same but mean voxel is diff)
#
# params(img::AbstractArray{T<:Distribution})

function intent_stats(intent::Distribution,
                      img::AbstractArray,
                      header::Dict{String,Any})
    fp = fieldnames(intent)
    np = length(fp)
    if ndims(img) < 5
        if np == 0
            img.properties["header"]["intent"] = intent()
        elseif np == 1
            img.properties["header"]["intent"] = intent(hdr.intent_p1)
        elseif np == 2
            img.properties["header"]["intent"] = intent(hdr.intent_p1,
                                                        intent_p2)
        elseif np == 3
            img.properties["header"]["intent"] = intent(hdr.intent_p1,
                                                        hdr.intent_p2,
                                                        hdr.intent_p3)
        end
    else
        # TODO 5th axis for parameters
        # probably do something like RGB but with Distribution parameters
    end
    ImageMeta(img, properties=Dict{String,Any}(
                                  "filetype" => ftype,
                                  "header" => header))
end

function setintent_stats!(hdr::NiftiHeader, img::ImageMeta)
    hdr.intent_code = NIFTI_INTENT_REVERSE[img.properties["header"]["intent"]]
    fp = fieldnames(NIFTI_INTENT_REVERSE[img.properties["header"]["intent"]])
    np = length(fp)
    if ndims(img) < 5
        if np == 0
            img.properties["header"]["intent"] = intent()
        elseif np == 1
            img.properties["header"]["intent"] = intent(hdr.intent_p1)
        elseif np == 2
            img.properties["header"]["intent"] = intent(hdr.intent_p1,
                                                        intent_p2)
        elseif np == 3
            img.properties["header"]["intent"] = intent(hdr.intent_p1,
                                                        hdr.intent_p2,
                                                        hdr.intent_p3)
        end
    else
        # TODO 5th axis for parameters
        # probably do something like RGB but with Distribution parameters
    end
end
