
function nii2_tests(img)
    @testset "NIfTI Interface" begin
        @test slicecode(img) == "Unkown"
        @test sliceduration(img) == 0
        @test slicestart(img) == 1  # zero in header
        @test sliceend(img) == 1  # zero in header

        @test frequencydim(img) == 1
        @test phasedim(img) == 1
        @test slicedim(img) == 1

        #=
        intent
        intentname
        intentparams
        =#
        @test scaleintercept(img) == 0.0
        @test scaleslope(img) == 1.0

        @test intentparams(img) == (0.0, 0.0, 0.0)
        @test qformcode(img)  == :MNI152
    #    @test qform(img) # TODO
        @test sform(img) ==  [-1.0  0.0  0.0    90.0
                               0.0  1.0  0.0  -126.0
                               0.0  0.0  1.0   -72.0
                               0.0  0.0  0.0    -1.0]
    end

    @testset "AbstractArray Interface" begin
        @test ndims(img) == 3
        @test size(img) == (182, 218, 182)
        @test length(img) == 7221032
        @test eltype(img) == Float32

        # TODO need to test actual values in image
    end


    @testset "ImageFormat Interface" begin
        @test all(spatunits(img) .== u"mm")
        @test description(img)[1:6] == "FSL3.3"

        @test auxfiles(img) == ["\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"]
        @test data_offset(img) == 544
        @test calmin(img) == 3000
        @test calmax(img) == 8000
    end

    @testset "Images Interface" begin
        @test ustrip(pixelspacing(img)[1]) == 1
        @test ustrip(pixelspacing(img)[2]) == 1
        @test ustrip(pixelspacing(img)[3]) == 1
        @test sdims(img) == 3
        @test nimages(img) == 1
        @test timedim(img) == 0
        @test size_spatial(img) == (182, 218, 182)
        @test spatialorder(img) == (:R2L, :P2A, :I2S)
    end
end

@testset "NIfTI-2 read" begin
    nii2_tests(load(File(format"NII", "data/MNI152_T1_1mm_nifti2.nii.gz")))
end

@testset "NIfTI-2 save" begin
    img = load(File(format"NII", "data/MNI152_T1_1mm_nifti2.nii.gz"))
    save(TEMPNII_FILE, img, version=2)
    img = load(File(format"NII", TEMPNII_FILE))
    nii2_tests(img)
end

@testset "" begin
    nii2_tests(metadata(File(format"NII", "data/MNI152_T1_1mm_nifti2.nii.gz")))
end


#=
Time Axis Shift = 0.0
Quaternion Parameters:  b = 0.0  c = 1.0  d = 0.0
Quaternion Offsets:  x = 90.0  y = -126.0  z = -72.0

Intent Code = 0

=#
