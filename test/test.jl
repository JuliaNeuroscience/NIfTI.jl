using NIFTI, Base.Test

file = niftiread("test/example4d.nii")

# Header
@test tr(file.header) == 2000 # Actually an error in the file AFAIK
@test_approx_eq voxelsize(file.header) Float32[2.0, 2.0, 2.2]
@test size(file) == (128, 96, 24, 2)

# Content
@test file.raw[65, 49, 13, :][:] == [265, 266]
@test vox(file, 64, 48, 12, :)[:] == [265, 266]
@test vox(file, 69, 56, 13, :)[:] == [502, 521]