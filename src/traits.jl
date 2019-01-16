Base.ndims(hdr::NiftiHeader) = hdr.dim[1]
Base.size(hdr::NiftiHeader, i::Int) = Int(hdr.dim[i+1])
Base.size(hdr::NiftiHeader) = convert(Tuple{Vararg{Int}}, hdr.dim[2:hdr.dim[1]+1])

function Base.eltype(hdr::NiftiHeader)
    t = get(NiftiDatatypes, hdr.datatype, nothing)
    if t == nothing
        @error "Unsupported data type $(hdr.datatype)"
    end
    t
end

voxoffset(hdr::NiftiHeader) = Int16(hdr.vox_offset)
"""
    spatunits(ImageMeta)
    spatunits(NiftiHeader)
    spatunits(AbstractArray)

Retrieves tuple of Unitful units for spatial dimensions. If Unitful units are
not set it returns `Unitful.FreeUnits{(),Unitful.Dimensions{()}`, which
produces a blank line in REPL (similar to `nothing`).
"""
spatunits(hdr::NiftiHeader) = get(NiftiUnits, hdr.xyzt_units & 0x07, 1)
function spatunits(img::AbstractArray)
    axv = axisvalues(img)
    return map(i -> unit(axv[i][1]), coords_spatial(img))
end

"""
    timeunits(ImageMeta)
    timeunits(NiftiHeader)
    timeunits(AbstractArray)

Retrieves Unitful units for time dimension. If time dimension is not present or
Unitful units are not set it returns `nothing`
"""
timeunits(hdr::NiftiHeader) = get(NiftiUnits, hdr.xyzt_units & 0x38, 1)
function timeunits(img::AbstractArray)
    ta = timeaxis(img)
    if ta == nothing
        return nothing
    else
        return unit(ta[1])
    end
end
"""
    description(ImageMeta)
    description(NiftiHeader)
    description(AbstractArray)

Retrieves description field that may say whatever you like.
"""
description(hdr::NiftiHeader) = String([hdr.descrip...])
description(img::ImageMeta) = @get img "description" ""
description(A::AbstractArray) = ""

"""
    auxfile(ImageMeta)
    auxfile(NiftiHeader)
    auxfile(AbstractArray)

Retrives string for auxiliary file associated with the image.
"""
auxfile(hdr::NiftiHeader) = String([hdr.aux_file...])
auxfile(s::ImageMeta) = @get s "auxfile" ""
auxfile(A::AbstractArray) = ""

qformcode(hdr::NiftiHeader) = get(NiftiXForm, hdr.qform_code, :Unkown)
sformcode(hdr::NiftiHeader) = get(NiftiXForm, hdr.sform_code, :Unkown)

xform(img::ImageMeta; form::Symbol=:sform_code) = @get img string(form) :Unkown
xform(a::AbstractArray) = :Unknown

"""
    qform(ImageMeta)
    qform(NiftiHeader)
    qform(AbstractArray)
"""
function qform(hdr::NiftiHeader)
    qform(hdr.qform_code, hdr.quatern_b, hdr.quatern_c, hdr.quatern_d,
          hdr.qoffset_x,  hdr.qoffset_y, hdr.qoffset_z, hdr.pixdim[2],
          hdr.pixdim[3],  hdr.pixdim[4])
end
function qform(qform_code::I, qb::T, qc::T, qd::T, qx::T, qy::T, qz::T, dx::T,
               dy::T, dz::T) where {I<:Integer,T<:AbstractFloat}
    if qform_code <= 0
        # TODO double check this:s
        # if not nifti or qform_code <= 0, use grid spacing for qto_xyz
        return ((T[dx,  0,  0, 0]...,),
                (T[ 0, dy,  0, 0]...,),
                (T[ 0,  0, dz, 0]...,),
                (T[ 0,  0,  0, 1]...,))
    else
        # use the quaternion-specified transformation
        qfac = dx < T(0.0) ? T(-1.0) : T(1.0)
        return quat2affine(qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac)
    end
end
function qform(A::ImageMeta)
    hdr = get(A.properties, "header", nothing)
    if isa(hdr, Nothing)
        return qform(data(A))
    else
        return get(hdr, "qform", sform(data(A)))
    end
end
qform(A::AbstractArray) = (spacedirections(A)..., (Int[0, 0, 0, 1]...))

"""
    sform(ImageMeta)
    sform(NiftiHeader)
    sform(AbstractArray)
"""
function sform(hdr::NiftiHeader)
    T = eltype(hdr.srow_x)
    if hdr.sform_code > 0
        # set the sto transformation from srow_*[]
        return (hdr.srow_x, hdr.srow_y, hdr.srow_z, (T[0,0,0,1]...,))
    else
        # not nifti or sform_code <= 0, then no sto transformation
        return 0
    end
end
function sform(A::ImageMeta)
    hdr = get(A.properties, "header", nothing)
    if isa(hdr, Nothing)
        return sform(data(A))
    else
        return get(hdr, "sform", sform(data(A)))
    end
end
sform(A::AbstractArray) = (spacedirections(A)..., (Int[0, 0, 0, 1]...))


"""
    frequencydim(ImageMeta)
    frequencydim(NiftiHeader)
    frequencydim(AbstractArray)
"""
frequencydim(hdr::NiftiHeader) = (diminfo(hdr) >> 2) & int8(3)
function code(A::ImageMeta)
    hdr = get(A.properties, "header", nothing)
    if isa(hdr, Nothing)
        return frequencydim(data(A))
    else
        return get(hdr, "frequency", frequencydim(data(A)))
    end
end
frequencydim(A::AbstractArray) = 0

"""
    phasedim(ImageMeta)
    phasedim(NiftiHeader)
    phasedim(AbstractArray)
"""
phasedim(hdr::NiftiHeader) = (diminfo(hdr) >> 2) & int8(3)
function phasecode(A::ImageMeta)
    hdr = get(A.properties, "header", nothing)
    if isa(hdr, Nothing)
        return phasedim(data(A))
    else
        return get(hdr, "phasedim", phasedim(data(A)))
    end
end
phasedim(A::AbstractArray) = 0

"""
    slicedim(ImageMeta)
    slicedim(NiftiHeader)
    slicedim(AbstractArray)
"""
slicedim(hdr::NiftiHeader) = (diminfo(hdr) >> 4) & int8(3)
function slicedim(A::ImageMeta)
    hdr = get(A.properties, "header", nothing)
    if isa(hdr, Nothing)
        return slicedim(data(A))
    else
        return get(hdr, "slicedim", slicedim(data(A)))
    end
end
slicedim(A::AbstractArray) = 0

diminfo(hdr::NiftiHeader) = header.dim_info
function diminfo(A::AbstractArray)
    fdim = frequencydim(A)
    pdim = phasedim(A)
    sdim = slicedime(A)
    if fdim > 3 || fdim < 0
        error("Invalid frequency dimension $(fdim)")
    elseif pdim > 3 || pdim < 0
        error("Invalid phase dimension $(dim_info[2])")
    elseif sdim > 3 || sdim < 0
        error("Invalid slice dimension $(dim_info[3])")
    end
    dim_info = Int8(fdim | (pdim << 2) | (sdim << 4))
end

"""
    slicecode(ImageMeta)
    slicecode(NiftiHeader)
    slicecode(AbstractArray)
"""
slicecode(hdr::NiftiHeader) = get(NiftiSliceCodes, hdr.slice_code, "Unkown")
function slicecode(A::ImageMeta)
    hdr = get(A.properties, "header", nothing)
    if isa(hdr, Nothing)
        return slicecode(data(A))
    else
        return get(hdr, "slicecode", slicecode(data(A)))
    end
end
slicecode(A::AbstractArray) = "Unkown"

"""
    sliceduration(ImageMeta)
    sliceduration(NiftiHeader)
    sliceduration(AbstractArray)
"""
sliceduration(hdr::NiftiHeader) = hdr.slice_duration
function sliceduration(A::ImageMeta)
    hdr = get(A.properties, "header", nothing)
    if isa(hdr, Nothing)
        return sliceduration(data(A))
    else
        return get(hdr, "sliceduration", sliceduration(data(A)))
    end
end
sliceduration(A::AbstractArray) = 0.0

"""
    slicestart(ImageMeta)
    slicestart(NiftiHeader)
    slicestart(AbstractArray)
"""
slicestart(hdr::NiftiHeader) = hdr.slice_start
function slicestart(A::ImageMeta)
    hdr = get(A.properties, "header", nothing)
    if isa(hdr, Nothing)
        return slicestart(data(A))
    else
        return get(hdr, "slicestart", slicestart(data(A)))
    end
end
slicestart(A::AbstractArray) = 0.0

"""
    sliceend(ImageMeta)
    sliceend(NiftiHeader)
    sliceend(AbstractArray)
"""
sliceend(hdr::NiftiHeader) = hdr.slice_end
function sliceend(A::ImageMeta)
    hdr = get(A.properties, "header", nothing)
    if isa(hdr, Nothing)
        return sliceend(data(A))
    else
        return get(hdr, "sliceend", sliceend(data(A)))
    end
end
