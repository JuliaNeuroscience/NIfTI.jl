## Usage

To read a NIFTI file:

```julia
nii = niftiread("my.nii")
```

The header is in `nii.header`; NIFTI extensions are in `nii.extensions`; the raw volume is in `nii.volume.raw`.

To mmap the NIFTI file:

```julia
nii = niftiread("my.nii", true)
```