const SIZEOF_HDR = Int32(348)

const NIfTI_DT_BITSTYPES = Dict{Int16,Type}([
    (Int16(2), UInt8),
    (Int16(4), Int16),
    (Int16(8), Int32),
    (Int16(16), Float32),
    (Int16(32), ComplexF32),
    (Int16(64), Float64),
    (Int16(256), Int8),
    (Int16(512), UInt16),
    (Int16(768), UInt32),
    (Int16(1024), Int64),
    (Int16(1280), UInt64),
    (Int16(1792), ComplexF64)
])
const NIfTI_DT_BITSTYPES_REVERSE = Dict{Type,Int16}()
for (k, v) in NIfTI_DT_BITSTYPES
    NIfTI_DT_BITSTYPES_REVERSE[v] = k
end

const NP1_MAGIC = (0x6e,0x2b,0x31,0x00)
const NI1_MAGIC = (0x6e,0x69,0x31,0x00)

# Conversion factors to mm/ms
# http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/xyzt_units.html
const SPATIAL_UNIT_MULTIPLIERS = [
    1000,   # 1 => NIfTI_UNITS_METER
    1,      # 2 => NIfTI_UNITS_MM
    0.001   # 3 => NIfTI_UNITS_MICRON
]
const TIME_UNIT_MULTIPLIERS = [
    1000,   # NIfTI_UNITS_SEC
    1,      # NIfTI_UNITS_MSEC
    0.001,  # NIfTI_UNITS_USEC
    1,      # NIfTI_UNITS_HZ
    1,      # NIfTI_UNITS_PPM
    1       # NIfTI_UNITS_RADS
   ]
