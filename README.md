NIfTI.jl
=======

[![Build Status](https://travis-ci.org/JuliaIO/NIfTI.jl.svg?branch=master)](https://travis-ci.org/JuliaIO/NIfTI.jl?branch=master)

## Usage

To read a NIfTI file:

```julia
using NIfTI
ni = niread("my.nii")
ni = niread("my.nii.gz") # gzipped NIfTI files are detected automatically
```

The header is in `nii.header`; NIfTI extensions are in `nii.extensions`; the raw
volume is in `nii.raw`.

To mmap the NIfTI file:

```julia
ni = niread("my.nii", mmap=true)
```

To get the TR and voxel size:
```julia
vsize = voxel_size(ni.header)    # In mm
tr = time_step(ni.header)        # In ms
```

To get the value of the volume along a given dimension:
```julia
d = vox(ni, x, y, z, t)       # Scaled by slope and intercept given in header, zero-based indexes
d = ni[x, y, z, t]            # Scaled by slope and intercept given in header, one-based indexes
d = ni.raw[x, y, z, t]        # Unscaled, one-based indexes
```
Colons works everywhere, even with `vox`

To write a volume:
```julia
niwrite("my.nii", ni)
```

It is also possible to construct a new volume from scratch; see the
`NIVolume` constructor in the source for documentation.
