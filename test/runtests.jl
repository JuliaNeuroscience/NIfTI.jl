# TODO
# - test io on example
# - test orientation on LR/RL image
# - test NIfTI-2 io
# - test endian values
# - test intent
#   - statistics
using FileIO, NIfTI, ImageMetadata, ImageCore, ImageAxes, GZip, Test, Unitful
using NIfTI.ImageFormats
import AxisArrays

import NIfTI: slicecode, sliceduration, slicestart, sliceend, frequencydim, phasedim,
              slicedim, qformcode, qform, sform, scaleslope, scaleintercept, intentparams

function extractto(gzname, out)
    open(out, "w") do io
        gzopen(gzname) do gz
            write(io, read(gz))
        end
    end
end

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


#
if FileIO.unknown((query(NII)))
    add_format(format"NII", (), [".nii", ".img"], [:NIfTI])
    #add_format(format"ANZ", detectnii, [".hdr"], [:NIfTI])
end



@testset "NIfTI-1" begin
    include("nii1.jl")
end

# data from: https://nifti.nimh.nih.gov/pub/dist/data/nifti2/
@testset "NIfTI-2" begin
    include("nii2.jl")
end


#@test_throws ErrorException load(GZIPPED_NII; mmap=true)
#@test_throws ErrorException load(GZIPPED_HDR; mmap=true)

# Test writing
#vol = NiftiImage()
#write(TEMP_FILE, vol)
#read(TEMP_FILE)

# Clean up
rm(NII)
rm(HDR)
rm(IMG)
