# NIfTI-1 not gzipped

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


for file_name in (NII, GZIPPED_NII, HDR, GZIPPED_HDR)
    img = niread(file_name)

    @testset "NIfTI Interface" begin
        @test slicecode(img) == "Unkown"
        @test sliceduration(img) == 0
        @test slicestart(img) == 0
        @test sliceend(img) == 23
        @test frequencydim(img) == 1
        @test phasedim(img) == 2
        @test slicedim(img) == 3
        @test qformcode(img)  == :Scanner_anat
        @test qform(img)  # TODO
        @test sform(img) == ((-2.00, 0.00,  0.00  117.86),
                             ( 0.00, 1.97, -0.36, -35.72),
                             ( 0.00, 0.32,  2.17,  -7.25))
    end

    @testset "ImageFormat Interface" begin
        @test spatunits(img) == u"mm"
        @test timeunits(img) == u"ms"
        @test description(img) == "FSL3.3\0 v2.25 NIfTI-1 Single file format\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
        @test auxfile(img) == "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
        @test data_offset(img) == 416
    end

    @testset "AbstractArray Interface" begin
        @test ndims(img) == 4
        @test size(img) == (128, 96, 24, 2)
        @test length(img) == prod((128, 96, 24, 2))
        @test eltype(img) == Int16
        # @test img[64, 48, 13, 1] == 496 FIXME
        @test img[64, 48, 12, :] == [265, 266]
        @test img[69, 56, 13, :] == [502, 521]
    end

    @testset "Images Interface" begin
        @test ustrip(pixelspacing(img)[1]) == 2.0
        @test ustrip(pixelspacing(img)[2]) == 2.0
        # TODO
        @test ustrip(pixelspacing(img)[3]) == 2.1999990940093994
        @test sdims(img) == 3
        @test nimages(img) == 2

        @test ustrip(timeaxis(img)[end]) == 2000000
        @test ustrip(timeaxis(img)[1]) == 0
        @test timedim(img) == 4
        #@test indices_spatial(img) == ((1.0:2.0:255.0)*u"mm",
        #                               (1.0:2.0:191.0)*u"mm",
        #                               (1.0:2.1999990940093994:51.59997916221619)*u"mm")
        @test size_spatial(img) == (128, 96, 24)
        @test spatialorder(img) == (:R2L, :P2A, :I2S)
    end

    if file_name == NII
        niwrite()
    end
end

# big endian with orientation
ALR = joinpath(dirname(@__FILE__), "data/avg152T1_LR_nifti.nii.gz")
ARL = joinpath(dirname(@__FILE__), "data/avg152T1_RL_nifti.nii.gz")

# LR
s = load(file; returntype=NiftiStream)
@test size(s) == (91,109,91)
@test spatialorder(s) == (:R2L, :P2A, :I2S)
close(s)

# RL
s = load(f; returntype=NiftiStream)
@test size(s) == (91,109,91)
@test axisnames(s) == (:L2R, :P2A, :I2S)
close(s)

