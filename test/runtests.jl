using NIfTI, GZip, Base.Test

function extractto(gzname, out)
	open(out, "w") do io
		gzopen(gzname) do gz
			write(io, readall(gz))
		end
	end
end

# single file storage
const GZIPPED_NII = joinpath(dirname(@__FILE__), "example4d.nii.gz")
const NII = "$(tempname()).nii"
extractto(GZIPPED_NII, NII)

# dual file storage
const GZIPPED_HDR = joinpath(dirname(@__FILE__), "example4d.hdr.gz")
hdr_stem = tempname()
const HDR = "$hdr_stem.hdr"
const IMG = "$hdr_stem.img"
extractto(GZIPPED_HDR, HDR)
extractto(joinpath(dirname(@__FILE__), "example4d.img.gz"), IMG)

for (fname, mmap) in ((NII, false), (NII, true), (HDR, false), (HDR, true),
	                  (GZIPPED_NII, false), (GZIPPED_HDR, false))
	file = niread(fname, mmap=mmap)

	# Header
	@test time_step(file.header) == 2000000 # Actually an error in the file AFAIK
	@test_approx_eq voxel_size(file.header) Float32[2.0, 2.0, 2.2]
	@test size(file) == (128, 96, 24, 2)

	# Content
	@test file.raw[65, 49, 13, :][:] == [265, 266]
	@test vox(file, 64, 48, 12, :)[:] == [265, 266]
	@test vox(file, 69, 56, 13, :)[:] == [502, 521]
end

@test_throws ErrorException niread(GZIPPED_NII; mmap=true)
@test_throws ErrorException niread(GZIPPED_HDR; mmap=true)

# Test writing
const TEMP_FILE = "$(tempname()).nii"
vol = NIVolume()
niwrite(TEMP_FILE, vol)
niread(TEMP_FILE)

# Clean up
rm(NII)
rm(HDR)
rm(IMG)
rm(TEMP_FILE)
