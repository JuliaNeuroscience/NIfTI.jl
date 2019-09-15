"""
    intentname(x)

Returns the name for the the image intent specified in the NIfTI header.
"""
intentname(x::Any) = getter(x, "intentname", String, String(fill(UInt8(0), 16)))

intentname!(x::Any, val::AbstractString) = setter!(x, "intentname", String(val), String)

"""
    scaleslope(x) -> Float64

The values stored in each voxel can be scaled, allowing storage of voxels as
smaller datatypes (`scaleslope * stored_value + scaleintercept -> actual_value`).
These values are ignored for RGB(A) data types.
"""
scaleslope(x::Any) = getter(x, "scaleslope", Float64, one(Float64))
scaleslope!(x::Any, val::Real) = setter!(x, "scaleslope", Float64(val), Float64)

"""
    scaleintercept(x) -> Float64

The values stored in each voxel can be scaled, allowing storage of voxels as
smaller datatypes (`scaleslope * stored_value + scaleintercept -> actual_value`).
These values are ignored for RGB(A) data types.
"""
scaleintercept(x::Any) = getter(x, "scaleintercept", Float64, one(Float64))
scaleintercept!(x::Any, val::Real) = setter!(x, "scaleintercept", Float64(val), Float64)

"""
    quaternb(x) -> Float64

Returns the quaternion b parameter.
"""
quaternb(x::Any) = getter(x, "quaternb", Float64, zero(Float64))
quaternb!(x::Any, val::Real) = setter!(x, "quaternb", Float64(val), Float64)

"""
    quaternc(x) -> Float64

Returns the quaternion c parameter.
"""
quaternc(x::Any) = getter(x, "quaternc", Float64, zero(Float64))
quaternc!(x::Any, val::Real) = setter!(x, "quaternc", Float64(val), Float64)

"""
    quaternd(x) -> Float64

Returns the quaternion d parameter.
"""
quaternd(x::Any) = getter(x, "quaternd", Float64, zero(Float64))
quaternd!(x::Any, val::Real) = setter!(x, "quaternd", Float64(val), Float64)

"""
    intentparams(x) -> NTuple{3,Float64}

Returns the parameters associated with the NIfTI intent.
"""
intentparams(x::Any) = getter(x, "intentparams", NTuple{3,Float64}, (zero(Float64),zero(Float64),zero(Float64)))

"""
    intentparams!(x, val)

Set the `intentparams` property.
"""
intentparams!(x::Any, val::NTuple{3,Float64}) = setter!(x, "intentparams", val, NTuple{3,Float64})
intentparams!(x::Any, val::NTuple{3,Real}) = setter!(x, "intentparams", Float64.(val), NTuple{3,Float64})



# NIfTI doesn't deal with 2D images so we convert to 3D. This assumes that the top-left
# 2x2 is the affine and and the far right is the linear transformation. The linear
# transformation has to be moved to the 4th column.

# TODO document slicecode
"""
    slicecode(x) -> String
"""
slicecode(x::Any) = getter(x, "slicecode", String, "Unkown")

# TODO document slicecode!
"""
    slicecode!(x, val)
"""
slicecode!(x::Any, val::AbstractString) = setter!(x, "slicecode", val, String)

function numeric2slicecode(i::Integer)
    if i == 1
        return "Sequential+Increasing"
    elseif i == 2
        return "Sequential+Decreasing"
    elseif i == 3
        return "Alternating+Increasing"
    elseif i == 4
        return "Alternating+Decreasing"
    elseif i == 5
        return "Alternating+Increasing#2"
    elseif i == 6
        return "Alternating+Decreasing#2"
    else
        return "Unkown"
    end
end

function slicecode2numeric(i::AbstractString)
    if i == "Sequential+Increasing"
        return 1
    elseif i == "Sequential+Decreasing"
        return 2
    elseif i == "Alternating+Increasing"
        return 3
    elseif i == "Alternating+Decreasing"
        return 4
    elseif i == "Alternating+Increasing#2"
        return 5
    elseif i == "Alternating+Decreasing#2"
        return 6
    else
        return 0
    end
end

diminfo(x::Any) = to_diminfo(freqdim(x), phasedim(x), slicedim(x))
function to_diminfo(f::Int, p::Int, s::Int)
    return ((((Int8(f - 1)) & 0x03)     ) |
            (((Int8(p - 1)) & 0x03) << 2) |
            (((Int8(s - 1)) & 0x03) << 4))
end


