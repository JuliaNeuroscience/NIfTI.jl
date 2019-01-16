# NIfTI-1 not gzipped
io = load(NII)
s = NiftiStream(io)

f = "test/data/example4d.nii"
#f = "test/data/example4d.nii.gz" # MethodError: no method matching unread(::TranscodingStreams.TranscodingStream{CodecZlib.GzipDecompressor,IOStream}, ::Int32)
img = niread(f)

s = NoopStream(open(f))
hdr = niread(s, NIfTI.Nifti1Header, false)
ext = niread(s, hdr, false, NiftiExtension)
props = properties(hdr, ext)


io = open(f)
img = niread(io, Array)

# until FileIO implementation
if isa(query(f),File{DataFormat{:UNKNOWN}})
    img = niread(NII)
else
    img = load(NII)
end

@testset "NiftiVolume" begin
    vol = NiftiVolume(s)
    @test vol[64, 48, 12, :] == [265, 266]
    @test vol[69, 56, 13, :] == [502, 521]
end

@testset "ImageMeta" begin
    @test ndims(img) == 4
    @test size(img) == (128, 96, 24, 2)
    @test length(img) == prod((128, 96, 24, 2))
    @test eltype(img) == Int16

    # nifti traits (implemented only in this package)
    @test spatunits(img)
    @test timeunits(img)
    @test frequencydim(img) == 1
    @test phasedim(img) == 2
    @test slicedim(img) == 3
    # TODO
    @test sliceinfo(img) == (slice_code = "Unkown",
                             slice_duration = 0.0f0,
                             slice_start = 0,
                             slice_end = 23)
    @test description(img) == "FSL3.3\0 v2.25 NIfTI-1 Single file format\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
    @test auxfile(img) == "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
    # TODO
    @test xform(img)  == :Scanner_anat

    # Image traits
    su = u"mm"
    @test spatunits(img) = (su, su, su,)
    @test ustrip(pixelspacing(img)[1]) == 2.0
    @test ustrip(pixelspacing(img)[2]) == 2.0
    # TODO
    @test ustrip(pixelspacing(img)[3]) == 2.1999990940093994
    @test sdims(img) == (128, 96, 24)
    @test nimages(img) == 2

    # TODO: getting (2000.0)
    @test ustrip(timeaxis(img)[end]) == 2000000
    @test ustrip(timeaxis(img)[1]) == 0
    @test timedim(img) == 4
    #@test indices_spatial(img) == ((1.0:2.0:255.0)*u"mm",
    #                               (1.0:2.0:191.0)*u"mm",
    #                               (1.0:2.1999990940093994:51.59997916221619)*u"mm")
    @test size_spatial(img) == (128, 96, 24)
    @test spacedirections(img) == ((-2.00, 0.00,  0.00  117.86),
                                   ( 0.00, 1.97, -0.36, -35.72),
                                   ( 0.00, 0.32,  2.17,  -7.25))
    @test spatialorder(img) == (:R2L, :P2A, :I2S)
    #@test voxoffset(img) == 416
    # @test img[64, 48, 13, 1] == 496 FIXME
    @test img[64, 48, 12, :] == [265, 266]
    @test img[69, 56, 13, :] == [502, 521]
end

@testset "Write ImageMeta as NIfTI" begin
    save(TEMPNII_FILE, img)
end

@testset "ImageMeta Mapped" begin
    img = load(file; mode="r", mmap=true, returntype=ImageMeta)
    img_traits_test(img)
end
