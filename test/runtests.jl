# TODO
# - test io on example
# - test orientation on LR/RL image
# - test NIfTI-2 io
# - test endian values
# - test intent
#   - statistics
using NIfTI, ImageMetadata, GZip
using Test

function extractto(gzname, out)
    open(out, "w") do io
        gzopen(gzname) do gz
            write(io, read(gz))
        end
    end
end

# single file storage
const GZIPPED_NII = joinpath(dirname(@__FILE__), "data/example4d.nii.gz")
const NII = "$(tempname()).nii"
extractto(GZIPPED_NII, NII)

const TEMPNII_FILE = "$(tempname()).nii"

# dual file storage
const GZIPPED_HDR = joinpath(dirname(@__FILE__), "data/example4d.hdr.gz")
hdr_stem = tempname()
const HDR = "$hdr_stem.hdr"
const IMG = "$hdr_stem.img"
extractto(GZIPPED_HDR, HDR)
extractto(joinpath(dirname(@__FILE__), "data/example4d.img.gz"), IMG)

@testset "NII1" begin
    include("readnii1.jl")
end
@testset "NII1GZ" begin
    include("readnii1gz.jl")
end

@testset "HDR1" begin
    include("readhdr1.jl")
end

@testset "HDR1GZ" begin
    include("readhdr1gz.jl")
end

@testset "NII2" begin
    include("readnii2.jl")
end
@testset "NII2GZ" begin
    include("readnii2gz.jl")
end

# check header and ImageMeta
# all NiftiReader API has to work for this to work

function test2()
    @test true == true
    @test true == false
end

dim_info = 9
# writing
@test niunits(img)
@test nidim(img)
@test nipixdim(img)
@test nibitpix(img) == 16
@test voxoffset(img) == 416

f = "./data/avg152T1_LR_nifti.nii.gz"
s = load(file; returntype=NiftiStream)
@test size(s) == (91,109,91)
@test spatialorder(s) == (:R2L, :P2A, :I2S)
close(s)

f = "./data/avg152T1_RL_nifti.nii.gz"
s = load(f; returntype=NiftiStream)
@test size(s) == (91,109,91)
@test axisnames(s) == (:L2R, :P2A, :I2S)
close(s)







# Test NIfTI-2
# data from: https://nifti.nimh.nih.gov/pub/dist/data/nifti2/

@test_throws ErrorException load(GZIPPED_NII; mmap=true)
@test_throws ErrorException load(GZIPPED_HDR; mmap=true)

# Test writing
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
