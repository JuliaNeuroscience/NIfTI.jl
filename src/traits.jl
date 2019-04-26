const NiftiSliceCodes = Dict{Int,String}(
    0 => "Unkown",
    1 => "Sequential+Increasing",
    2 => "Sequential+Decreasing",
    3 => "Alternating+Increasing",
    4 => "Alternating+Decreasing",
    5 => "Alternating+Increasing#2",
    6 => "Alternating+Decreasing#2")

const NiftiSliceCodesReverse = Dict{String,Int16}()
for (k, v) in NiftiSliceCodes
    NiftiSliceCodesReverse[v] = k
end

"""
    scaleslope
"""
scaleslope(img::ImageFormat{format"NII"}) = img["header"]["scaleslope"]
scaleslope(A::AbstractArray) = 1

"""
    scaleintercept
"""
scaleintercept(img::ImageFormat{format"NII"}) = img["header"]["scaleintercept"]
scaleintercept(A::AbstractArray) = 0

diminfo(img::ImageFormat{format"NII"}) = img["header"]["diminfo"]
diminfo(A::AbstractArray) = Int8(0)

"""
    frequencydim
"""
frequencydim(s::AbstractArray) = (diminfo(s) >> 4) & Int8(3)

"""
    phasedim
"""
phasedim(s::AbstractArray) = (diminfo(s) >> 4) & Int8(3)

"""
    slicedim
"""
slicedim(s::AbstractArray) = (diminfo(s) >> 4) & Int8(3)

"""
    slicecode
"""
slicecode(img::ImageFormat{format"NII"}) = img["header"]["slicecode"]
slicecode(A::AbstractArray) = "Unkown"

"""
    sliceduration
"""
sliceduration(img::ImageFormat{format"NII"}) = img["header"]["sliceduration"]
sliceduration(A::AbstractArray) = 0

"""
    slicestart
"""
slicestart(img::ImageFormat{format"NII"}) = img["header"]["slicestart"]
slicestart(A::AbstractArray) = 0

"""
    sliceend
"""
sliceend(img::ImageFormat{format"NII"}) = img["header"]["sliceend"]
sliceend(A::AbstractArray) = 0

"""
    qform(img)
"""
function qform(img::ImageMeta)
    if haskey(properties(img), "nifti")
        return get(img["nifti"], "qform", qform(data(img)))
    else
        return qform(data(img))
    end
end

function qform(A::AbstractArray)
    SMatrix{4,4,Float64,16}(1.0, 0.0, 0.0, 0.0,
                            0.0, 1.0, 0.0, 0.0,
                            0.0, 0.0, 1.0, 0.0,
                            0.0, 0.0, 0.0, 1.0)
end

"""
    sform(A)
"""
function sform(img::ImageMeta)
    if haskey(properties(img), "nifti")
        return get(img["nifti"], "sform", sform(data(img)))
    else
        return sform(data(img))
    end
end

function sform(A::AbstractArray)
    SMatrix{4,4,Float64,16}(1.0, 0.0, 0.0, 0.0,
                            0.0, 1.0, 0.0, 0.0,
                            0.0, 0.0, 1.0, 0.0,
                            0.0, 0.0, 0.0, 1.0)
end

"""
qformcode
"""
qformcode(img::ImageFormat{format"NII"}) = img["header"]["qformcode"]
qformcode(A::AbstractArray) = :Unkown

"""
sformcode
"""
sformcode(img::ImageFormat{format"NII"}) = img["header"]["sformcode"]
sformcode(A::AbstractArray) = :Unkown
