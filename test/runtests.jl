# TODO
# - test io on example
# - test orientation on LR/RL image
# - test NIfTI-2 io
# - test endian values
# - test intent
#   - statistics
using NIfTI, ImageMetadata, FileIO
using Test

using FileIO
# FileIO stuff

function detectnii(io::IO)
    isgzipped = isgz(io)
    hdrio = isgzipped ? gzdopen(io) : io
    version, swap = checkfile(hdrio)
    if version == 1
        magic = seek(hdrio, number_to_magic)
        return magic == NP1_MAGIC ? true : false
    elseif version == 2
        magic = read(hdrio, Array{UInt8}(undef,8))
        return magic == NP1_MAGIC ? true : false
    else
        return false
    end
end

detecthdr(io::IO) = detectnii(io)

function detectimg(io::IO)
    img_file = filename(io)
    hdr_file = gethdr(img_file)
    open(hdr_file, "r") do io
        detectnii(io)
    end
end

add_format(format"NII", detectnii, [".nii"], [:NIfTI])
# HDR is already FileIO registered
#add_format(format"HDR", detecthdr, [".hdr"], [:NIfTI])
add_format(format"IMG", detectimg, [".img"], [:NIfTI])

#file = "/Volumes/SD/projects/NIfTI.jl/test/data/example4d.nii"

file = "/Volumes/SD/projects/NIfTI.jl/test/data/avg152T1_LR_nifti.nii.gz"
hdr = load(file; returntype=NiftiHeader)
sform, qform = NIfTI.getaffine(hdr)
@test NIfTI.mat2ori(qform) == (:L2R, :P2A, :I2S)
@test NIfTI.mat2ori(sform) == (:R2L, :P2A, :I2S)

file = "/Volumes/SD/projects/NIfTI.jl/test/data/avg152T1_RL_nifti.nii.gz"
hdr = load(file; returntype=NiftiHeader)
sform, qform = NIfTI.getaffine(hdr)
@test NIfTI.mat2ori(qform) == (:L2R, :P2A, :I2S)
@test NIfTI.mat2ori(sform) == (:L2R, :P2A, :I2S)

img = load(file)

function extractto(gzname, out)
    open(out, "w") do io
        gzopen(gzname) do gz
            write(io, read(gz))
        end
    end
end

img = NIfTI.load(file; mode="r", mmap=false)

# from
# https://nifti.nimh.nih.gov/nifti-1/data/avg152T1_LR_nifti_ntool.txt
function test_avgLR(img::ImageMeta)
    # check intent
    intent = getintent(img)
    intent.intent_p1 == 0.0
    intent.intent_p2 == 0.0
    intent.intent_p3 == 0.0
    intent.intent_code == 0
    intent.intent_name == 0  # this was left blank from the link

    # check dims
    size(img) == (91, 109, 91)
    # check pixdim (in actual hdr length(pixdim) == 8)
    pixelspacing(img) == (2.000000, 2.000000, 2.000000) * u""
    xyzt_units     = 10

    # check slicing
    img.properties["header"]["slice_duration"] == 0.000000
    img.properties["header"]["toffset"] == 0.000000
    img.properties["header"]["slice_end"] == 0
    img.properties["header"]["slice_code"] == 0
    img.properties["header"]["slice_start"] == 0

    # check affine
    img.properties["header"]["qform_code"] == 0
    img.properties["header"]["sform_code"] == 4

    # TODO figure out how this relates to spacedirections
    img["spacedirections"] == ((-2.000000 0.000000 0.000000 90.000000),
                               (0.000000 2.000000 0.000000 -126.000000),
                               (0.000000 0.000000 2.000000 -72.000000))

    img.properties["header"]["cal_max"] == 255.000000
    img.properties["header"]["cal_min"] == 0.000000

    NIfTI.NIFTI_DT_BITSTYPES_REVERSE[eltype(img)] == 2
    img.properties["header"]["bitpix"] == 8

    img.properties["descrip"] == "FSL3.2beta"
    get(img.properties["header"], "aux_file", nothing) == nothing

    # check extension TODO
    img.properties["header"]["extension"]
end

function test_example4d(img::ImageMeta)
    size(img) == (128, 96, 24)
    datatype == Int16
    pixelspacing(img) == (2.00, 2.00, 2.20)
    @test img[64, 48, 12, :] == [265, 266]
    @test img[69, 56, 13, :] == [502, 521]
end

# single file storage
const GZIPPED_NII = joinpath(dirname(@__FILE__), "data/example4d.nii.gz")
const NII = "$(tempname()).nii"
extractto(GZIPPED_NII, NII)

# dual file storage
const GZIPPED_HDR = joinpath(dirname(@__FILE__), "data/example4d.hdr.gz")
hdr_stem = tempname()
const HDR = "$hdr_stem.hdr"
const IMG = "$hdr_stem.img"
extractto(GZIPPED_HDR, HDR)
extractto(joinpath(dirname(@__FILE__), "data/example4d.img.gz"), IMG)

for (fname, mmap) in ((NII, false), (NII, true), (HDR, false), (HDR, true),
                      (GZIPPED_NII, false), (GZIPPED_HDR, false))

    img = load(fname; mmap=mmap)

    # Header
    @test eltype(img) == Int16
    @test pixelspacing(file.header) ≈ Float32[2.0, 2.0, 2.2] * u"mm"
    @test sdim(img) == (128, 96, 24)
    @test nimages(img) == 2
    @test timedim(img) == (2000000) * u"ms"

    # Content
    @test img[64, 48, 12, :] == [265, 266]
    @test img[69, 56, 13, :] == [502, 521]
    @test spacedirections(img) == ((-2.00, 0.00,  0.00  117.86),
                                   ( 0.00, 1.97, -0.36, -35.72),
                                   ( 0.00, 0.32,  2.17,  -7.25))
    # TODO
    @test axisnames(img) == ()

    @assert maximum(file) == maximum(file.raw)
end

# Test NIfTI-2
# data from: https://nifti.nimh.nih.gov/pub/dist/data/nifti2/

@test_throws ErrorException load(GZIPPED_NII; mmap=true)
@test_throws ErrorException load(GZIPPED_HDR; mmap=true)

# Test writing
const TEMP_FILE = "$(tempname()).nii"
vol = NiftiImage()
write(TEMP_FILE, vol)
read(TEMP_FILE)

# Site is currently down TODO: reintroduce this test when site is up
# Big endian
# const BE = "$(tempname()).nii"
# download("https://nifti.nimh.nih.gov/nifti-1/data/avg152T1_LR_nifti.nii.gz", BE)
img = niread("data/avg152T1_LR_nifti.nii.gz")
@test size(img) == (91,109,91)
axisnames(img) == (:L,:r,:f)

img = niread("data/avg152T1_RL_nifti.nii.gz")
axisnames(img) == (:,:,:)


# Clean up
rm(NII)
rm(HDR)
rm(IMG)
rm(TEMP_FILE)
# rm(BE)
