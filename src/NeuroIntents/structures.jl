
# TODO This is a terrible name "Estimate"
struct Estimate{T}
    p1::T
    p2::T
    p3::T
end
numerictag(::Type{<:Estimate}) = 1001
stringtag(::Type{<:Estimate}) = "NIFTI_INTENT_ESTIMATE"
intentaxis(::Type{<:Estimate}) = (:estimate,)

# TODO NeuroLabel
struct NeuroLabel <: AbstractLabel end
numerictag(::Type{<:NeuroLabel}) = 1002
stringtag(::Type{<:NeuroLabel}) = "NIFTI_INTENT_LABEL"
intentaxis(::Type{<:NeuroLabel}) = (:neurolabel,)

struct NeuroName <: AbstractLabel end
numerictag(::Type{<:NeuroName}) = 1003
stringtag(::Type{<:NeuroName}) = "NIFTI_INTENT_NEURONAME"
intentaxis(::Type{<:NeuroName}) = (:neuroname,)

numerictag(::Type{<:SMatrix}) = 1004
stringtag(::Type{<:SMatrix}) = "NIFTI_INTENT_GENMATRIX"
intentaxis(::Type{<:SMatrix}) = (:genmatrix,)

numerictag(::Type{<:Symmetric}) = 1005
stringtag(::Type{<:Symmetric}) = "NIFTI_INTENT_SYMMATRIX"
intentaxis(::Type{<:Symmetric}) = (:symmatrix,)

struct DisplacementVector end
numerictag(::Type{<:DisplacementVector}) = 1006
stringtag(::Type{<:DisplacementVector}) = "NIFTI_INTENT_DISPVECT"
intentaxis(::Type{<:DisplacementVector}) = (:displacement,)

numerictag(::Type{<:SVector}) = 1007
stringtag(::Type{<:SVector}) = "NIFTI_INTENT_VECTOR"
intentaxis(::Type{<:SVector}) = (:vector,)

numerictag(::Type{<:Point}) = 1008
stringtag(::Type{<:Point}) = "NIFTI_INTENT_POINTSET"
intentaxis(::Type{<:Point}) = (:point,)

numerictag(::Type{<:Triangle}) = 1009
stringtag(::Type{<:Triangle}) = "NIFTI_INTENT_TRIANGLE"
intentaxis(::Type{<:Triangle}) = (:triangle,)

numerictag(::Type{<:Quat}) = 1010
stringtag(::Type{<:Quat}) = "NIFTI_INTENT_QUATERNION"
intentaxis(::Type{<:Quat}) = (:quaternion,)

struct Dimensionless end
numerictag(::Type{<:Dimensionless}) = 1010
stringtag(::Type{<:Dimensionless}) = "NIFTI_INTENT_DIMLESS"
intentaxis(::Type{<:Dimensionless}) = (:dimensionaless,)

function structintent(i::Integer)
    if i == 1001
        return Estimate
    elseif i == 1002
        return NeuroLabel
    elseif i == 1003
        return NeuroNames
    elseif i == 1004
        return SMatrix
    elseif i == 1005
        return Symmetric
    elseif i == 1006
        return DisplacementVector
    elseif i == 1007
        return SVector
    elseif i == 1008
        return Point
    elseif i == 1009
        return Triangle
    elseif i == 1010
        return Quat
    elseif i == 1011
        return Dimensionless
    else
        return UnknownIntent
    end
end
