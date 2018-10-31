const SIZEOF_HDR1 = Int32(348)
const SIZEOF_HDR2 = Int32(540)

const NIFTI_DT_BITSTYPES = Dict{Int16,Type}([
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

const NIFTI_DT_BITSTYPES_REVERSE = Dict{Type,Int16}()
for (k, v) in NIFTI_DT_BITSTYPES
    NIFTI_DT_BITSTYPES_REVERSE[v] = k
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
    # (Int16(40), u"ppm"),  # parts per million TODO
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
    (Int16(1),  :L2R),
    (Int16(-1), :R2L),
    (Int16(2),  :P2A),
    (Int16(-2), :A2P),
    (Int16(3),  :I2S),
    (Int16(-3), :S2I)
])

ori2space = Dict{Symbol, Symbol}(
    :L2R => :L,
    :R2L => :R,
    :P2A => :P,
    :A2P => :A,
    :I2S => :I,
    :S2I => :S
)
ori_lut = [[   1   -1    2   -2    3   -3]
           [:L2R :R2L :P2A :A2P :I2S :S2I]
           [  :L   :R   :P   :A   :I   :S]]

function isradview(axnames::NTuple{3,Symbol})
    axnames == (:L, :A, :S)
end

function isneuroview(axnames::NTuple{3,Symbol})
    axnames == (:R, :A, :S)
end

const NIFTI_SLICE = Dict{Int16, Symbol}([
    (Int16(0), :unkown),
    (Int16(1), :seq_inc),
    (Int16(2), :seq_dec),
    (Int16(3), :alt_inc),
    (Int16(4), :alt_dec),
    (Int16(5), :alt_inc2),
    (Int16(6), :alt_dec2)
   ])

const NIFTI_FTYPE = Dict{Int16, Symbol}([
    (Int16(0), :Analyze),
    (Int16(1), :NIfTI1_1),
    (Int16(2), :NIfTI1_2),
    (Int16(3), :ASCII),
    (Int16(4), :NIfTI2_1),
    (Int16(5), :NIfTI2_2)
])

const NIFTI_XFORM = Dict{Int16,Symbol}([
    (Int16(0), :unkown),
    (Int16(1), :scanner_anat),
    (Int16(2), :aligned_anat),
    (Int16(3), :talairach),
    (Int16(4), :mni152)
])

const NIFTI_XFORM_REVERSE = Dict{Symbol,Int16}()

for (k, v) in NIFTI_XFORM
    NIFTI_XFORM_REVERSE[v] = k
end

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
