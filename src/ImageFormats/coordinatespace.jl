struct CoordinateSpace{S} end

"""
    UnkownSpace
"""
const UnkownSpace = CoordinateSpace{:unkown}()


"""
    ScannerSpace

In scanner space.
"""
const ScannerSpace = CoordinateSpace{:scanner}()


"""
    AnatomicalSpace

equivalent of 'aligned' space in NIfTI standard.
"""
const AnatomicalSpace = CoordinateSpace{:anatomical}()

"""
    TailarachSpace

Tailarach space
"""
const TailarachSpace = CoordinateSpace{:tailarach}()

"""
    MNI152Space

MNI152 space
"""
const MNI152Space = CoordinateSpace{:MNI152}()

"""
    coordinatespace(x)

Return the coordinate space that `x` is in.
"""
function coordinatespace(x::Any)
    getter(x, "coordinatespace", CoordinateSpace, UnkownSpace)
end

"""
    coordinatespace!(x, val)

Set the coordinate space for `x` to `val`.
"""
function coordinatespace!(x::Any, val::CoordinateSpace)
    setter!(x, "coordinatespace", val, CoordinateSpace)
end


