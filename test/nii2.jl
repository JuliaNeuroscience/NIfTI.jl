"data/MNI152_T1_1mm_nifti2.nii.gz"


@testset "NIfTI Interface" begin
    @test slicecode(img) == "Unkown"
    @test sliceduration(img) == 0
    @test slicestart(img) == 1  # zero in header
    @test sliceend(img) == 1  # zero in header

    @test frequencydim(img) == 1
    @test phasedim(img) == 1
    @test slicedim(img) == 1

    intent
    intentname
    intentparams


    @test sliceintercept(img) = 0.0
    @test sliceslope(img) = 1.0

    @test intentparams(img) == (0.0, 0.0, 0.0)
    @test qformcode(img)  == :MNI152
    @test qform(img) # TODO
        @test sform(img) == ((-1.0, 0.0, 0.0    90.0),
                             ( 0.0, 1.0, 0.0, -126.0),
                             ( 0.0, 0.0, 1.0,  -72.0))
                         #TODO add fourth row
end


@testset "AbstractArray Interface" begin
    @test ndims(img) == 3
    @test size(img) == (182, 218, 182)
    @test length(img) == 7221032
    @test eltype(img) == Int16

    # TODO
    # @test img[64, 48, 13, 1] == 496 FIXME
    #@test img[64, 48, 12, :] == [265, 266]
    #@test img[69, 56, 13, :] == [502, 521]
end


@testset "ImageFormat Interface" begin
    @test spatunits(img) == u"mm"
    @test timeunits(img) == u"s"
    @test description(img) == "FSL3.3\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0
                               \0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0
                               \0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
    @test auxfile(img) == "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
    @test data_offset(img) == 544
    @test calmin(img) = 3000.0
    @test calmin(img) = 8000.0
end

@testset "Images Interface" begin
    @test ustrip(pixelspacing(img)[1]) == 1
    @test ustrip(pixelspacing(img)[2]) == 1
    @test ustrip(pixelspacing(img)[3]) == 1
    @test sdims(img) == 3
    @test nimages(img) == 1
    @test timedim(img) == 0
    @test indices_spatial(img) == (Base.OneTo(182), Base.OneTo(218), Base.OneTo(182))
    @test size_spatial(img) == (182, 218, 182)
    @test spatialorder(img) == (:R2L, :P2A, :I2S)
end

#=
Time Axis Shift = 0.0
Quaternion Parameters:  b = 0.0  c = 1.0  d = 0.0
Quaternion Offsets:  x = 90.0  y = -126.0  z = -72.0

Intent Code = 0

=#
