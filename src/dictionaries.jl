const SIZEOF_HDR1 = Int32(348)
const SIZEOF_HDR2 = Int32(540)

const NiftiDatatypes = Dict{Int16,Type}([
    (Int16(2), UInt8),
    (Int16(4), Int16),
    (Int16(8), Int32),
    (Int16(16), Float32),
    (Int16(32), ComplexF32),
    (Int16(64), Float64),
    (Int16(128), RGB),  # 3 8 bit bytes
    (Int16(256), Int8),
    (Int16(512), UInt16),
    (Int16(768), UInt32),
    (Int16(1024), Int64),
    (Int16(1280), UInt64),
    (Int16(1536), Float128),
    (Int16(1792), ComplexF64),
    (Int16(2048), ComplexF128),
    (Int16(2304), RGBA)  # 4 8 bit bytes
])

const NiftiDatatypesReverse = Dict{Type,Int16}()
for (k, v) in NiftiDatatypes
    NiftiDatatypesReverse[v] = k
end

const NP1_MAGIC = [0x6e,0x2b,0x31,0x00]
const NI1_MAGIC = [0x6e,0x69,0x31,0x00]
const NP2_MAGIC = [0x6e,0x2b,0x32,0x00, 0x0d, 0x0a,0x1a,0x0a]
const NI2_MAGIC = [0x6e,0x69,0x32,0x00, 0x0d, 0x0a,0x1a,0x0a]


#@unit ppm "ppm" PartsPerMillion 1//1000000 false
#Unitful.register(ppm)

const NiftiUnits = Dict([
    (Int16(1), u"m"),  # meter
    (Int16(2), u"mm"),  # millimeters
    (Int16(3), u"μm"),  # microns
    (Int16(8), u"s"),  # second
    (Int16(16), u"ms"),  # milliseconds
    (Int16(24), u"μs"), # microseconds
    (Int16(32), u"Hz"),  # hertz
    #(Int16(40), u"ppm"),  # parts per million
    (Int16(48), u"rad/s")  # radians per second
])

const NiftiUnitsReverse = Dict()
for (k, v) in NiftiUnits
    NiftiUnitsReverse[v] = k
end

const ANALYZE75_ORIENT = Dict{Symbol,Int16}([
    (:a75_transverse_unflipped, 0),
    (:a75_coronal_unflipped, 1),
    (:a75_sagittal_unflipped, 2),
    (:a75_transverse_flipped, 3),
    (:a75_coronal_flipped, 4),
    (:a75_sagittal_flipped, 5),
    (:a75_orient_unknown, 6)
])


"""
    xform(i) -> CoordinateSpace
"""
function xform(i::Union{Int16,Int32})
    if i == 1
        return ScannerSpace
    elseif i == 2
        return AnatomicalSpace
    elseif i == 3
        return TalairachSpace
    elseif i == 4
        return MNI152Space
    else
        return UnkownSpace
    end
end

numerictag(::Type{CoordinateSpace{Sym}}) where {Sym} = Int16(0)
numerictag(::Type{CoordinateSpace{:scanner}}) = Int16(1)
numerictag(::Type{CoordinateSpace{:anatomical}}) = Int16(2)
numerictag(::Type{CoordinateSpace{:tailarach}}) = Int16(3)
numerictag(::Type{CoordinateSpace{:MNI152}}) = Int16(4)
