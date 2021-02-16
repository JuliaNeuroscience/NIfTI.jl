
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

function image_tests(img)
    # Header
    @test time_step(img.header) == 2000000 # Actually an error in the file AFAIK
    @test voxel_size(img.header) ≈ Float32[2.0, 2.0, 2.2]
    @test size(img) == (128, 96, 24, 2)

    # Content
    @test img.raw[65, 49, 13, :][:] == [265, 266]
    @test vox(img, 64, 48, 12, :)[:] == [265, 266]
    @test vox(img, 69, 56, 13, :)[:] == [502, 521]

    @assert maximum(img) == maximum(img.raw)
end

image_tests(niread(NII, mmap=false))
image_tests(niread(NII, mmap=true))

file = joinpath(TEMP_DIR_NAME, "$(tempname()).nii")
niwrite(file, niread(NII))
image_tests(niread(NII, mmap=false))


image_tests(niread(HDR, mmap=false))
image_tests(niread(HDR, mmap=true))

image_tests(niread(GZIPPED_NII, mmap=false))
image_tests(niread(GZIPPED_HDR, mmap=false))

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
mask = [true false true; false true false; false true true]
mask_bitarray = BitArray(mask) # BitArray

niwrite(BOOL_WRITE, NIVolume(mask))
@test niread(BOOL_WRITE).raw == mask

niwrite(BIT_WRITE, NIVolume(mask_bitarray))
@test niread(BIT_WRITE).raw == mask_bitarray

# Open mmaped file for reading and writing
const WRITE = joinpath(TEMP_DIR_NAME, "$(tempname()).nii")
const VERIFY_WRITE = joinpath(TEMP_DIR_NAME, "$(tempname()).nii")
cp(NII, WRITE)
img = niread(WRITE; mmap=true, mode="r+")
img.raw[1,1,1,1] = 5
img.raw[:,2,1,1] = ones(size(img)[1])
cp(WRITE, VERIFY_WRITE)
@test niread(VERIFY_WRITE)[1,1,1,1] == 5
@test niread(VERIFY_WRITE)[:,2,1,1] == ones(size(img)[1])
# Site is currently down TODO: reintroduce this test when site is up
# Big endian
# const BE = "$(tempname()).nii"
# download("https://nifti.nimh.nih.gov/nifti-1/data/avg152T1_LR_nifti.nii.gz", BE)
img = niread(joinpath(dirname(@__FILE__), "data/avg152T1_LR_nifti.nii.gz"))
@test size(img) == (91,109,91)

GC.gc() # closes mmapped files

