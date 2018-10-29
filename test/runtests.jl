using NIfTI, GZip, ImageMetadata
using Test

include("src/NIfTI.jl")
using FileIO
file = "/Volumes/SD/projects/NIfTI.jl/test/data/example4d.nii"
img = load(file)

function extractto(gzname, out)
    open(out, "w") do io
        gzopen(gzname) do gz
            write(io, read(gz))
        end
    end
end

# from
# https://nifti.nimh.nih.gov/nifti-1/data/avg152T1_LR_nifti_ntool.txt
function test_avgLR(img::ImageMeta)
    # check intent
    intent = NIfTI.getintent(img)
    intent.intent_p1 == 0.0
    intent.intent_p2 == 0.0
    intent.intent_p3 == 0.0
    intent.intent_code == 0
    intent.intent_name == nothing # this was left blank from the link

    # check dims
    size(img) == (91, 109, 91)
    # check pixdim (in actual hdr length(pixdim) == 8)
    pixelspacing(img) == (2.000000, 2.000000, 2.000000)
    xyzt_units     => 10

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
             srow_x         => -2.000000 0.000000 0.000000 90.000000,
             srow_y         => 0.000000 2.000000 0.000000 -126.000000,
             srow_z         => 0.000000 0.000000 2.000000 -72.000000,

    img.properties["header"]["cal_max"] == 255.000000
    img.properties["header"]["cal_min"] == 0.000000

    NIfTI.NIFTI_DT_BITSTYPES_REVERSE[eltype(img)] == 2
    img.properties["header"]["bitpix"] == 8

    img.properties["header"]["scl_slope"] == 0.000000
    img.properties["header"]["scl_inter"] == 0.000000
    img.properties["descrip"] == "FSL3.2beta"
    img.properties["header"]["aux_file"] == nothing

    # check extension TODO
    img.properties["header"]["extension"]
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

# test load(.nii)
file = load(NII{format"NII"}; mode=false)

# test load(.hdr)

# test load(.img)

# test load(.nii.gz)

# test load(.hdr.gz)

# test save()

# test orientation
for (fname, mmap) in ((NII, false),
                      (NII, true),
                      (HDR, false),
                      (HDR, true),
                      (GZIPPED_NII, false),
                      (GZIPPED_HDR, false))
    load(fname{format"NII"}
    file = niread(fname, mmap=mmap)

    # Header
    @test time_step(file.header) == 2000000 # Actually an error in the file AFAIK
    @test voxel_size(file.header) ≈ Float32[2.0, 2.0, 2.2]
    @test size(file) == (128, 96, 24, 2)

    # Content
    @test file.raw[65, 49, 13, :][:] == [265, 266]
    @test vox(file, 64, 48, 12, :)[:] == [265, 266]
    @test vox(file, 69, 56, 13, :)[:] == [502, 521]

    @assert maximum(file) == maximum(file.raw)
end



@test_throws ErrorException niread(GZIPPED_NII; mmap=true)
@test_throws ErrorException niread(GZIPPED_HDR; mmap=true)

# Test writing
const TEMP_FILE = "$(tempname()).nii"
vol = NIVolume()
niwrite(TEMP_FILE, vol)
niread(TEMP_FILE)

# Site is currently down TODO: reintroduce this test when site is up
# Big endian
# const BE = "$(tempname()).nii"
# download("https://nifti.nimh.nih.gov/nifti-1/data/avg152T1_LR_nifti.nii.gz", BE)
img = niread("data/avg152T1_LR_nifti.nii.gz")
@test size(img) == (91,109,91)

# Clean up
rm(NII)
rm(HDR)
rm(IMG)
rm(TEMP_FILE)
# rm(BE)
