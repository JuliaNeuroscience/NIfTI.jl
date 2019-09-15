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
    orientationspace(x)

"""
orientationspace(x::T) where {T} = orientationspace(HasProperties(T), x)
orientationspace(::HasProperties{true}, x) = orientationspace(properties(x))
orientationspace(::HasProperties{false}, x) = UnkownSpace
