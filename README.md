NIfTI.jl
=======

## Usage

To read a NIfTI file:

```julia
using NIfTI
ni = niftiread("my.nii")
```

The header is in `nii.header`; NIfTI extensions are in `nii.extensions`; the raw
volume is in `nii.raw`.

To mmap the NIfTI file:

```julia
ni = niftiread("my.nii", true)
```

To get the TR and voxel size:
```julia
vsize = voxel_size(ni.header)     # In mm
tr = time_step(ni.header)        # In ms
```

To get the value of the volume along a given dimension:
```
d = vox(ni, x, y, z, t)       # Scaled by slope and intercept given in header,
zero-based indexes
d = ni[x, y, z, t]            # Scaled by slope and intercept given in header,
one-based indexes
d = ni.raw[x, y, z, t]        # Unscaled, one-based indexes
```
Colons works everywhere, even with `vox`

To write a volume:
```julia
niftiwrite("my.nii", ni)
```

It is also possible to construct a new volume from scratch; see the
`NIfTIVolume` constructor in the source for documentation.