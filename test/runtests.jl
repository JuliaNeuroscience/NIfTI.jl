using NIfTI, GZip
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

# dual file storage
const GZIPPED_HDR = joinpath(dirname(@__FILE__), "data/example4d.hdr.gz")
hdr_stem = tempname()
const HDR = "$hdr_stem.hdr"
const IMG = "$hdr_stem.img"
extractto(GZIPPED_HDR, HDR)
extractto(joinpath(dirname(@__FILE__), "data/example4d.img.gz"), IMG)

for (fname, mmap) in ((NII, false), (NII, true), (HDR, false), (HDR, true),
	                  (GZIPPED_NII, false), (GZIPPED_HDR, false))
	  file = niread(fname, mmap=mmap)

	# Header
	@test time_step(file.header) == 2000000 # Actually an error in the file AFAIK
  @test voxel_size(file.header) ≈ Float32[2.0, 2.0, 2.2]
	@test size(file) == (128, 96, 24, 2)

	# Content
	@test file.raw[65, 49, 13, :][:] == [265, 266]
	@test vox(file, 64, 48, 12, :)[:] == [265, 266]
	@test vox(file, 69, 56, 13, :)[:] == [502, 521]

	@assert maximum(file) == maximum(file.raw)
end

@test_throws ErrorException niread(GZIPPED_NII; mmap=true)
@test_throws ErrorException niread(GZIPPED_HDR; mmap=true)

# Test writing
const TEMP_FILE = "$(tempname()).nii"
vol = NIVolume()
niwrite(TEMP_FILE, vol)
niread(TEMP_FILE)

# Site is currently down TODO: reintroduce this test when site is up
# Big endian
# const BE = "$(tempname()).nii"
# download("https://nifti.nimh.nih.gov/nifti-1/data/avg152T1_LR_nifti.nii.gz", BE)
img = niread("data/avg152T1_LR_nifti.nii.gz")
@test size(img) == (91,109,91)

# Clean up
rm(NII)
rm(HDR)
rm(IMG)
rm(TEMP_FILE)
# rm(BE)
