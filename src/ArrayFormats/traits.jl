"""
    dataoffset(x) -> Int

Returns the IO stream offset to data given an type instance. Defaults to 0.
"""
dataoffset(x::Any) = getter(x, "dataoffset", Int, 1)

"""
    dataoffset!(x, val)

Set the the `data_offset` property.
"""
dataoffset!(x::Any, i::Integer) = setter!(x, "dataoffset", i, Int)

"""

    auxfiles(x) -> Vector{String}

Retrieves string for auxiliary file associated with the image.
"""
auxfiles(x::Any) = getter(x, "auxfiles", Vector{String}, [""])

"""
    auxfiles!(x, val)

Sets the `auxfiles` property. `val` should be a `String` or `Vector{String}`.
"""
auxfiles!(x::Any, i::Vector{String}) = setter!(x, "auxfiles", i, Vector{String})

"""
    srcfile(x) -> String

Retrieves the file name that the image comes from.
"""
srcfile(x::Any) = getter(x, "srcfile", String, "")

"""
    srcfile!(x, f::String)

Change `srcfile` property.
"""
srcfile!(x::Any, val::AbstractString) = setter!(x, "srcfile", val, String)

"""
    description(x) -> String

Retrieves description field that may say whatever you like.
"""
description(x::Any) = getter(x, "description", String, "")

"""
    description!(x, descrip::String)

Change description defined in an properties type.
"""
description!(x::Any, val::AbstractString) = setter!(x, "description", val, String)


