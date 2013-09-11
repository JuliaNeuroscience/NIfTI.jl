using NIfTI, Base.Test

file = niread("test/example4d.nii")

# Header
@test time_step(file.header) == 2000000 # Actually an error in the file AFAIK
@test_approx_eq voxel_size(file.header) Float32[2.0, 2.0, 2.2]
@test size(file) == (128, 96, 24, 2)

# Content
@test file.raw[65, 49, 13, :][:] == [265, 266]
@test vox(file, 64, 48, 12, :)[:] == [265, 266]
@test vox(file, 69, 56, 13, :)[:] == [502, 521]

# Test writing
const TEMP_FILE = "/tmp/my.nii"
vol = NIVolume()
niwrite(TEMP_FILE, vol)
niread(TEMP_FILE)
rm(TEMP_FILE)