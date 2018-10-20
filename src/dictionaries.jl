const SIZEOF_HDR1 = Int32(348)
const SIZEOF_HDR2 = Int32(540)

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
const NP2_MAGIC = (0x6e,0x2b,0x32,0x00, 0x0d, 0x0a,0x1a,0x0a)

const NIFTI_UNITS = Dict([
    (Int16(0), nothing),
    (Int16(1), u"m"),  # meter
    (Int16(2), u"mm"),  # millimeters
    (Int16(3), u"μm"),  # microns
    (Int16(8), u"s"),  # second
    (Int16(16), u"ms"),  # milliseconds
    (Int16(24), u"μs"), # microseconds
    (Int16(32), u"Hz"),  # hertz
    (Int16(40), u"ppm"),  # parts per million
    (Int16(48), u"rad/s")  # radians per second
])

const ANALYZE75_ORIENT = Dict{Symbol,Int16}([
    (:a75_transverse_unflipped, 0),
    (:a75_coronal_unflipped, 1),
    (:a75_sagittal_unflipped, 2),
    (:a75_transverse_flipped, 3),
    (:a75_coronal_flipped, 4),
    (:a75_sagittal_flipped, 5),
    (:a75_orient_unknown, 6)
])

const NIFTI_ORIENTATION = Dict{Int16, Symbol}([
    (Int16(1), :L2R),
    (Int16(-1), :R2L),
    (Int16(2), :P2A),
    (Int16(-2), :A2P),
    (Int16(3), :I2S),
    (Int16(-3), :S2I)
])

const NIFTI_SLICE = Dict{Int16, Symbol}([
    (Int16(0), :unkown),
    (Int16(1), :seq_inc),
    (Int16(2), :seq_dec),
    (Int16(3), :alt_inc),
    (Int16(4), :alt_dec),
    (Int16(5), :alt_inc2),
    (Int16(6), :alt_dec2)
   ])

const NIFTI_FTYPE = Dict{NTuple{4,UInt8}
, Symbol}([
    (Int16(0), :Analyze),
    (Int16(1), :NIfTI1_1),
    (Int16(2), :NIfTI1_2),
    (Int16(3), :ASCII),
    (Int16(4), :NIfTI2_1),
    (Int16(5), :NIfTI2_2)
])

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

const NIFTI_INTENT = Dict{Int16,Symbol}([
    (Int16(0), :none),
    (Int16(2), :correl),
    (Int16(3), :ttest),
    (Int16(4), :ftest),
    (Int16(5), :zscore),
    (Int16(6), :chisq),
    (Int16(7), :beta),
    (Int16(8), :binom),
    (Int16(9), :gamma),
    (Int16(10), :poisson),
    (Int16(11), :normal),
    (Int16(12), :ftest_nonc),
    (Int16(13), :chisq_nonc),
    (Int16(14), :logistic),
    (Int16(15), :laplace),
    (Int16(16), :uniform),
    (Int16(17), :ttest_nonc),
    (Int16(18), :weibull),
    (Int16(19), :chi),
    (Int16(20), :invgauss),
    (Int16(21), :extval),
    (Int16(22), :pval),
    (Int16(23), :logpval),
    (Int16(24), :log10pval),

    # not stats codes
    (Int16(1001), :estimate),
    (Int16(1002), :label),
    (Int16(1003), :neuroname),
    (Int16(1004), :genmatrix),
    (Int16(1005), :symmatrix),
    (Int16(1006), :dispvect),
    (Int16(1007), :vector),
    (Int16(1008), :pointset),
    (Int16(1009), :triangle),
    (Int16(1010), :quaternion),
    (Int16(1011), :dimless),

    # gifti datasets
    (Int16(2001), :time_series),
    (Int16(2002), :node_index),
    (Int16(2003), :rgb_vector),
    (Int16(2004), :rgba_vector),
    (Int16(2005), :shape)
])

const NIFTI_XFORM = Dict{Int16,Symbol}([
    (Int16(0), :unkown),
    (Int16(1), :scanner_anat),
    (Int16(2), :aligned_anat),
    (Int16(3), :talairach),
    (Int16(4), :mni152)
])

const NIFTI_ECODE = Dict{Int16, Symbol}([
    (Int16(0), :ignore),
    (Int16(2), :dicom),
    (Int16(4), :afni),
    (Int16(6), :comment),
    (Int16(8), :xcede),
    (Int16(10), :jimdiminfo),
    (Int16(12), :workflow_fwds),
    (Int16(14), :freesurfer),
    (Int16(16), :pypickly),
    (Int16(18), :mind_ident),
    (Int16(20), :b_value),
    (Int16(22), :spherical_direction),
    (Int16(24), :dt_component),
    (Int16(26), :shc_degreeorder),
    (Int16(28), :voxbo),
    (Int16(30), :caret),
    (Int16(32), :cifti),
    (Int16(34), :variable_frame_timing),
    (Int16(38), :eval),
    (Int16(40), :matlab)
])

