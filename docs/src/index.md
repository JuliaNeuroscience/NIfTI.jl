
```@meta
CurrentModule = NIfTI
```

# NIfTI

Documentation for [NIfTI](https://github.com/JuliaNeuroscience/NIfTI.jl).

```@docs
NIfTI.freqdim
Accepts a NIfTI.NIfTI1Header
Extracts frequency dimension from dim_info field. (dim_info & 0x03)

NIfTI.phasedim
Accepts a NIfTI.NIfTI1Header
Extracts frequency dimension from dim_info field. (dim_info >> 0x02) & 0x03

NIfTI.slicedim
Accepts a NIfTI.NIfTI1Header
Extracts frequency dimension from dim_info field. (dim_info >> 0x04)

NIfTI.slice_start
NIfTI.slice_end
NIfTI.slice_duration
NIfTI.sdims
NIfTI.voxel_size
```
