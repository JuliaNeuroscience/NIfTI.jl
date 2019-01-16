function scale_data(A::AbstractArray)
    slope = 0
    intercept = 0
    return slope, intercept
end

# AbstractArray → NIfTI
function NiftiHeader(a::AbstractArray{T,N}, ext::NiftiExtension; v::Int=1,
                     cal_max::Union{AbstractFloat,Nothing}=nothing,
                     cal_min::Union{AbstractFloat,Nothing}=nothing,
                     scl_slope::Float64=Float64(0),
                     scl_inter::Float64=Float64(0)) where {T,N}
    nitimeunits = unit(timeaxis(a)[1])
    nispatunits = ([unit(i[1]) for i in indices_spatial(a)]...,)

    to_diminfo(a) 
    # Gets spatial and time units for header and converts them to type `Int` for use
    # in a NIfTI file. (currently assumes all spatial units are same).
    ss = get(NiftiUnitsReverse, spatunits(A)[1], 0)
    tt = get(NiftiUnitsReverse, timeunits(A), 0)
    niunits = (((char)(ss)) & 0x07) | (((char)(tt)) & 0x38)


    # get nifti formatted dims
    dims = ones(Int16, 8)
    dims[1] = N
    dims[2:dims[1]+1] = [size(a)...]
    dims = (dim...,)

    # Retrieves the element type of an image and converts it for encoding the element
    # type stored in a NIfTI file.
    niftiT = get(NIFTI_DT_BITSTYPES_REVERSE, T, nothing)
    if niftiT == nothing
        @error "Unsupported data type $(T)"
    end

    # nibitpix
    nibitipix = Int16(sizeof(T)*8)

    # pixdim
    pixdim = ones(Int16, 8)
    pixdim[1] = Int16(N)
    pixdim[2:dim[1]+1] = [ustrip.(pixelspacing(a))...]
    pixdim = (pixdim...,)

    si = sliceinfo(a)

    affinetuple = spacedirections(a)

    if v == 1
        sizeof_hdr = SIZEOF_HDR1
        offsettype = Int32
    elseif v == 2
        sizeof_hdr = SIZEOF_HDR2
        offsettype = Int64
    else
        @error "Version $v is not a supported NIfTI version"
    end

    if isempty(ext)
        voxoffset = offsettype(sizeof_hdr)
    else
        voxoffset = offsettype(mapreduce(esize, +, ext) + sizeof_hdr)
    end

    # time offset
    ta = timeaxis(a)
    if ta == nothing
        timeoffset = 1
    else
        timeoffset = ta[1]
    end

    niintent = NiftiIntentsReverse(imageintent(a))

    NiftiHeader(dim_info(a), dims, niintent.p1, niintent.p2, niintent.p3,
                niintent.code, niftiT, nibitpix, si.slice_start,
                pixdim, voxoffset, scl_slope, scl_inter, si.slice_end,
                get(NIFTI_SLICE_REVERSE, si.slice_code, 0), niunits,
                ifelse(cal_max == nothing, maximum(a), cal_max),
                ifelse(cal_min == nothing, minimum(a), cal_min),
                si.slice_duration, timeoffset, string_tuple(description(a), 80),
                string_tuple(auxfile(a), 24), 0,  # default to blank qform
                get(NiftiXFormReverse, xform(a), 0),  # xform_code
                zero(Float32), zero(Float32), zero(Float32), zero(Float32),
                zero(Float32), zero(Float32),
                affinetuple[1], affinetuple[2], affinetuple[3],
                niintent.name; v=v)
end



niwrite_volume(io::IO, hdr::NiftiHeader, a::A) where {A<:ImageMeta} = niwrite_volume(io, hdr, data(a))
niwrite_volume(io::IO, hdr::NiftiHeader, a::A) where {A<:AxisArray} = niwrite_volume(io, hdr, data(a))

function niwrite_volume(io::IO, hdr::NiftiHeader, a::A) where {A<:AbstractArray{T,N}} where {T<:RGBA,N}
    # 4 x dim1 x dim2 x ... → 4 * dim1 x dim2 x ....
    rawsize = (size(hdr, 1) * 4, [hdr.dim[i] for i in 3:ndims(hdr)]...,)
    a = reinterpret(Float32, reshape(channelview(a), rawsize))
    write(io, a)
end

function niwrite_volume(io::IO, hdr::NiftiHeader, a::A) where {A<:AbstractArray{T,N}} where {T<:RGB,N}
    # 3 x dim1 x dim2 x ... → 3 * dim1 x dim2 x ....
    rawsize = (size(hdr, 1) * 3, [hdr.dim[i] for i in 3:ndims(hdr)]...,)
    a = reinterpret(Float32, reshape(channelview(a), rawsize))
    write(io, a)
end

function niwrite_volume(io::IO, hdr::NiftiHeader, a::A) where {A<:AbstractArray{T,N}} where {T<:Gray,N}
    write(io, channelview(a))
end

function niwrite_volume(io::IO, hdr::NiftiHeader, a::A) where {A<:AbstractArray{T,N}} where {T<:Number,N}
    write(io, a)
end

function niwrite(io::IO, a::AbstractArray; v::Int=1, scale_data::Bool=false,
                 cal_max::Union{AbstractFloat,Nothing}=nothing,
                 cal_min::Union{AbstractFloat,Nothing}=nothing)
    ext = NiftiExtension(a)

    # TODO need to convert array to same spatial axis
    # TODO implement this
    if scale_data
        scl_slope, scl_inter = scale_data(a)
    else
        scl_slope = 0.0
        scl_inter = 0.0
    end

    hdr = NiftiHeader(a, ext; v=v, cal_max=cal_max, cal_min=cal_min,
                      scl_slope=scl_sclope, scl_inter=scl_inter)

    write(io, typeof(hdr))
    write(io, ext)

    niwrite_volume(io, hdr, a)
end

function niwrite(io, A::AbstractArray{T,N}; kwarg...) where {T,N}
    niwrite(io, NiftiExtension(A; kwarg...), NiftiExtension(A))
end
function niwrite(io::IO, hdr::NiftiExtension, ext::NiftiExtension, A::Array)
    niwrite(io, hdr)
    niwrite(io, ext)
    niwrite(io, A)
end
