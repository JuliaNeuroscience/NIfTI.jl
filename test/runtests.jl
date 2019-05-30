# TODO
# - test io on example
# - test orientation on LR/RL image
# - test NIfTI-2 io
# - test endian values
# - test intent
#   - statistics
using NIfTI, ImageMetadata, ImageCore, ImageAxes, GZip, Test, Unitful
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
