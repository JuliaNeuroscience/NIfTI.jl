
using CodecZlib
using NIfTI
using Test
using TranscodingStreams

function extractto(gzname, out)
    open(out, "w") do io_out
        open(gzname, "r") do io_in
            gz = GzipDecompressorStream(io_in)
            write(io_out, read(gz))
        end
    end
end

const TEMP_DIR_NAME = mktempdir()

# empty file test
const EMPTY_NII = joinpath(dirname(@__FILE__), "data/empty.nii.gz")
@test_throws EOFError niread(EMPTY_NII)

# single file storage
const GZIPPED_NII = joinpath(dirname(@__FILE__), "data/example4d.nii.gz")
const NII = joinpath(TEMP_DIR_NAME, "$(tempname()).nii")
extractto(GZIPPED_NII, NII)

# dual file storage
const GZIPPED_HDR = joinpath(dirname(@__FILE__), "data/example4d.hdr.gz")
hdr_stem = joinpath(TEMP_DIR_NAME, tempname())
const HDR = "$hdr_stem.hdr"
const IMG = "$hdr_stem.img"
extractto(GZIPPED_HDR, HDR)
extractto(joinpath(dirname(@__FILE__), "data/example4d.img.gz"), IMG)

function image_tests(fname, mmap)
    img = niread(fname, mmap=mmap)

    # Header
    @test time_step(img.header) == 2000000 # Actually an error in the file AFAIK
    # @test all(isapprox.(voxel_size(img.header), (2.0, 2.0, 2.2)))
    vs1, vs2, vs3 = voxel_size(img.header)
    @test isapprox(vs1, Float32(2.0))
    @test isapprox(vs2, Float32(2.0))
    @test isapprox(vs3, Float32(2.2))
    @test size(img) == (128, 96, 24, 2)

    # Content
    @test img.raw[65, 49, 13, :][:] == [265, 266]
    @test vox(img, 64, 48, 12, :)[:] == [265, 266]
    @test vox(img, 69, 56, 13, :)[:] == [502, 521]

    @assert maximum(img) == maximum(img.raw)

    @test getaffine(img) ≈ [
        -2.0 6.714715653593746e-19 9.081024511081715e-18 117.8551025390625
        6.714715653593746e-19 1.9737114906311035 -0.35552823543548584 -35.72294235229492
        8.25548088896093e-18 0.3232076168060303 2.171081781387329 -7.248798370361328
        0.0 0.0 0.0 1.0
    ]

    @test NIfTI.get_qform(img) ≈ [
        -2.0 7.75482f-26 -6.93824f-27 117.855
        7.75482f-26 1.97371 -0.355528 -35.7229
        6.30749f-27 0.323208 2.17108 -7.2488
        0.0 0.0 0.0 1.0
    ]
    @test NIfTI.orientation(img) == (:right, :posterior, :inferior)
end

image_tests(NII, false)
image_tests(NII, true)
image_tests(HDR, false)
image_tests(HDR, true)
image_tests(GZIPPED_NII, false)
image_tests(GZIPPED_HDR, false)

@testset "Header Field Accessors" begin
    img = niread(GZIPPED_NII)
    @test NIfTI.freqdim(img) == 1
    @test NIfTI.phasedim(img) == 2
    @test NIfTI.slicedim(img) == 3
    @test NIfTI.slice_start(img) == 1
    @test NIfTI.slice_end(img) == 24
    @test NIfTI.slice_duration(img) == 0
end

@test_throws ErrorException niread(GZIPPED_NII; mmap=true)
@test_throws ErrorException niread(GZIPPED_HDR; mmap=true)

# Test writing
const TEMP_FILE = joinpath(TEMP_DIR_NAME, "$(tempname()).nii")
vol = NIVolume()
niwrite(TEMP_FILE, vol)
niread(TEMP_FILE)

const TEMP_GZIPPED_FILE = joinpath(TEMP_DIR_NAME, "$(tempname()).nii.gz")
niwrite(TEMP_GZIPPED_FILE, vol)
niread(TEMP_GZIPPED_FILE)

# Write and read DT_BINARY
const BOOL_WRITE = joinpath(TEMP_DIR_NAME, "$(tempname()).nii")
const BIT_WRITE = joinpath(TEMP_DIR_NAME, "$(tempname()).nii")
mask = rand(Bool, 3, 5, 7) # Array{Bool}
mask_bitarray = BitArray(mask) # BitArray
niwrite(BOOL_WRITE, NIVolume(mask))
niwrite(BIT_WRITE, NIVolume(mask_bitarray))
@test niread(BOOL_WRITE).raw == mask
@test niread(BIT_WRITE).raw == mask_bitarray

# Write and read INT16 volume
const INT16_WRITE = joinpath(TEMP_DIR_NAME, "$(tempname()).nii")
vol_INT16 = rand(Int16, 3, 5, 7) # Array{Int16}
niwrite(INT16_WRITE, NIVolume(vol_INT16))
@test niread(INT16_WRITE).raw == vol_INT16

# Open mmaped file for reading and writing
const WRITE = joinpath(TEMP_DIR_NAME, "$(tempname()).nii")
const VERIFY_WRITE = joinpath(TEMP_DIR_NAME, "$(tempname()).nii")
cp(NII, WRITE)
img = niread(WRITE; mmap=true, mode="r+")
img.raw[1, 1, 1, 1] = 5
img.raw[:, 2, 1, 1] = ones(size(img)[1])
cp(WRITE, VERIFY_WRITE)
@test niread(VERIFY_WRITE)[1, 1, 1, 1] == 5
@test niread(VERIFY_WRITE)[:, 2, 1, 1] == ones(size(img)[1])
# Site is currently down TODO: reintroduce this test when site is up
# Big endian
# const BE = "$(tempname()).nii"
# download("https://nifti.nimh.nih.gov/nifti-1/data/avg152T1_LR_nifti.nii.gz", BE)
img = niread(joinpath(dirname(@__FILE__), "data/avg152T1_LR_nifti.nii.gz"))
@test size(img) == (91, 109, 91)

GC.gc() # closes mmapped files

@test NIfTI._dir2ori(-1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0) == (:right, :posterior, :inferior)

@test NIfTI._dir2ori(1.0, 0.0, 0.0,
    0.0, -1.0, 0.0,
    0.0, 0.0, 1.0) == (:left, :anterior, :inferior)


@test NIfTI._dir2ori(1.0, 0.0, 0.0,
    0.0, -1.0, 0.0,
    0.0, 0.0, -1.0) == (:left, :anterior, :superior)


