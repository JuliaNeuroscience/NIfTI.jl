# NIfTI-1 not gzipped

#=
if epoch_idx - last_improvement >= 5 && opt.eta > 1e-6
    opt.eta /= 10.0
    @warn(" -> Haven't improved in a while, dropping learning rate to $(opt.eta)!")

    # After dropping learning rate, give it a few epochs to improve
    last_improvement = epoch_idx
end
=#

# single file storage

function nii1_tests(img::AbstractArray)
    @testset "NIfTI Interface" begin
        @test slicecode(img) == "Unkown"
        @test sliceduration(img) == 0
        @test slicestart(img) == 1

        # dim_info = 57
        @test sliceend(img) == 24
        @test frequencydim(img) == 2
        @test phasedim(img) == 3
        @test slicedim(img) == 4
        @test qformcode(img)  == :Scanner_anat
        # TODO Something is messed up with these values
        # The NIfTI format tries to normalize the qform to orientation. Therefore, it's not
        # always (rarely) used by other libraries, but we need it for when we interact with
        # DICOMs or some versions of FSL (I think).
        #NIfTI.orientation(qform(img)) == AxisArrays.axisnames(img)[1:3]

        @test sform(img) ≈ [ -2.0                    6.714715653593746e-19  9.081024511081715e-18  117.8551025390625
                              6.714715653593746e-19  1.9737114906311035    -0.35552823543548584   -35.72294235229492
                              8.25548088896093e-18   0.3232076168060303     2.171081781387329     -7.248798370361328
                              0.0                    0.0                    0.0                   -1.0]
    end

    @testset "ImageFormat Interface" begin
        @test all(spatunits(img) .== u"mm")
        @test timeunits(img) == u"s"
        @test description(img) == "FSL3.3\0 v2.25 NIfTI-1 Single file format\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
        @test auxfiles(img) == ["\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"]
        # TODO this should be tested on a stream, not an image
        #@test data_offset(img) == 416
    end

    @testset "AbstractArray Interface" begin
        @test ndims(img) == 4
        @test size(img) == (128, 96, 24, 2)
        @test length(img) == prod((128, 96, 24, 2))
        @test eltype(img) == Int16
        # @test img[64, 48, 13, 1] == 496 FIXME
        @test img[65, 49, 13, 1] == 265
        @test img[65, 49, 13, 2] == 266
        @test img[70, 57, 14, 1] == 502
        @test img[70, 57, 14, 2] == 521
    end

    @testset "Images Interface" begin
        @test ustrip(pixelspacing(img)[1]) == 2.0
        @test ustrip(pixelspacing(img)[2]) == 2.0
        # TODO
        @test ustrip(pixelspacing(img)[3]) == 2.1999990940093994
        @test sdims(img) == 3
        @test nimages(img) == 2

        @test ustrip(timeaxis(img)[end]) == 2000
        @test ustrip(timeaxis(img)[1]) == 0
        @test timedim(img) == 4
        #@test indices_spatial(img) == ((1.0:2.0:255.0)*u"mm",
        #                               (1.0:2.0:191.0)*u"mm",
        #                               (1.0:2.1999990940093994:51.59997916221619)*u"mm")
        @test size_spatial(img) == (128, 96, 24)
        @test spatialorder(img) == (:R2L, :P2A, :I2S)
    end
end
function hdr1_tests(img::AbstractArray)
    @testset "NIfTI Interface" begin
        @test slicecode(img) == "Unkown"
        @test sliceduration(img) == 0
        @test slicestart(img) == 1
        # These fields aren't written in the HDR example data so we test defaults
        @test sliceend(img) == 1
        @test frequencydim(img) == 1
        @test phasedim(img) == 1
        @test slicedim(img) == 1
        @test qformcode(img)  == :Scanner_anat
        # TODO Something is messed up with these values
        # The NIfTI format tries to normalize the qform to orientation. Therefore, it's not
        # always (rarely) used by other libraries, but we need it for when we interact with
        # DICOMs or some versions of FSL (I think).
        #NIfTI.orientation(qform(img)) == AxisArrays.axisnames(img)[1:3]

        @test sform(img) ≈ [ -2.0                    6.714715653593746e-19  9.081024511081715e-18  117.8551025390625
                              6.714715653593746e-19  1.9737114906311035    -0.35552823543548584   -35.72294235229492
                              8.25548088896093e-18   0.3232076168060303     2.171081781387329     -7.248798370361328
                              0.0                    0.0                    0.0                   -1.0]
    end

    @testset "ImageFormat Interface" begin
        @test all(spatunits(img) .== u"mm")
        @test timeunits(img) == u"s"
        @test description(img) == "FSL5.0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
        @test auxfiles(img) == ["\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"]
        # TODO this should be tested on a stream, not an image
        #@test data_offset(img) == 416
    end

    @testset "AbstractArray Interface" begin
        @test ndims(img) == 4
        @test size(img) == (128, 96, 24, 2)
        @test length(img) == prod((128, 96, 24, 2))
        @test eltype(img) == Int16
        # @test img[64, 48, 13, 1] == 496 FIXME
        @test img[65, 49, 13, 1] == 265
        @test img[65, 49, 13, 2] == 266
        @test img[70, 57, 14, 1] == 502
        @test img[70, 57, 14, 2] == 521
    end

    @testset "Images Interface" begin
        @test ustrip(pixelspacing(img)[1]) == 2.0
        @test ustrip(pixelspacing(img)[2]) == 2.0
        # TODO
        @test ustrip(pixelspacing(img)[3]) == 2.1999990940093994
        @test sdims(img) == 3
        @test nimages(img) == 2

        @test ustrip(timeaxis(img)[end]) == 2000
        @test ustrip(timeaxis(img)[1]) == 0
        @test timedim(img) == 4
        #@test indices_spatial(img) == ((1.0:2.0:255.0)*u"mm",
        #                               (1.0:2.0:191.0)*u"mm",
        #                               (1.0:2.1999990940093994:51.59997916221619)*u"mm")
        @test size_spatial(img) == (128, 96, 24)
        @test spatialorder(img) == (:R2L, :P2A, :I2S)
    end
end



@testset ".nii+mmap" begin
    img = load(NII, mmap=true)
    nii1_tests(img)
end

@testset ".nii.gz" begin
    img = load(File(format"NII", GZIPPED_NII))
    nii1_tests(img)
end

@testset "save" begin
    img = load(File(format"NII", GZIPPED_NII))
    save(TEMPNII_FILE, img)
    img = load(File(format"NII", TEMPNII_FILE))
    nii1_tests(img)
end

@testset ".hdr" begin
    img = load(File(format"NII", HDR), mmap=true)
    hdr1_tests(img)
end

@testset ".hdr.gz" begin
    img = load(File(format"NII", GZIPPED_HDR))
    hdr1_tests(img)
end

@testset "orientation" begin
    # big endian with orientation
    ALR = joinpath(dirname(@__FILE__), "data/avg152T1_LR_nifti.nii.gz")
    ARL = joinpath(dirname(@__FILE__), "data/avg152T1_RL_nifti.nii.gz")

    # LR
    img = load(File(format"NII", ALR))
    @test size(img) == (91,109,91)
    @test spatialorder(img) == (:R2L, :P2A, :I2S)

    # RL
    img = load(File(format"NII", ARL))
    @test size(img) == (91,109,91)
    @test spatialorder(img) == (:L2R, :P2A, :I2S)
end

#img2 = permutedims(img, (2,1,3,4))
