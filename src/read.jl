function checkfile(io::TranscodingStream)
    ret = read(io, 4)
    hdrsz = reinterpret(Int32, ret)
    if hdrsz[1] == Int32(348)
        TranscodingStreams.unread(io, ret)
        return 1, false
    elseif hdrsz[1] == Int32(540)
        TranscodingStreams.unread(io, ret)
        return 2, false
    elseif hdrsz[1] == ntoh(Int32(348))
        TranscodingStreams.unread(io, ret)
        return 1, true
    elseif hdrsz[1] == ntoh(Int32(540))
        TranscodingStreams.unread(io, ret)
        return 2, true
    else
        TranscodingStreams.unread(io, ret)
        return 0, false
    end
end

function isgz(io::IO)
    gzbits = read(io, 2)
    ret = gzbits == [0x1F,0x8B]
    seek(io, 0)
    ret
end

function getimg(f::AbstractString)
    img_file = replace(f, ".hdr" => ".img")
    if isfile(out)
        return img_file
    else
        error("NIfTI file is dual file storage, but $img_file does not exist")
    end
end

function gethdr(f::AbstractString)
    hdr_file = replace(f, ".img" => ".hdr")
    if isfile(out)
        return hdr_file
    else
        error("NIfTI file is dual file storage, but $hdr_file does not exist")
    end
end

"""
NiftiSchema

* qformcode/sformcode: orientation relative to qform affine and sform affine respectively
* sformcode: orientation relative to sform affine
    - Unkown
    - Scanner_anat
    - Aligned_anat
    - Talairach)
    - MNI152
* qform
* sform:
* spacedirections: affine orientation see @ref(`spacedirections`) for more information
* `description`: see [`description`!](@ref)
* `auxfile`: see [`auxfile`!](@ref)
* "calibration": `cal_max` and `cal_min` refer to the maximum and minimu display intensity.
  Values above `cal_max` and below `cal_min` should be binarized to a single
  color (e.g., all values above `cal_max` are white and all values below `cal_min`
  are black).

"""
struct NiftiSchema{S,T,intent}
    axes::Tuple
    props::Dict{String,Any}
    needswap::Bool
end

function NiftiSchema(io::IO)
    version, needswap = checkfile(io)
    if version == 1
        hdr = niread(io, Nifti1Header, needswap)
    elseif version == 2
        hdr = niread(io, Nifti2Header, needswap)
    else
        @error "File does not have NIfTI v1 or v2 header formatting."
    end

    N = ndims(hdr)
    sz = size(hdr)
    su = spatunits(hdr)
    props = Dict{String,Any}(
        "description" => String([hdr.descrip...]),
        "auxfile" => String([hdr.aux_file...]),
        "header" => Dict{String,Any}(
            "intentname" => String([hdr.intent_name...]),
            "slice_duration" => hdr.slice_duration,
            "extension" => niread(io, hdr, needswap, NiftiExtension),
            "qformcode" => sformcode(hdr),
            "qform" => qform(hdr),
            "sformcode" => sformcode(hdr),
            "sform" => sform(hdr),
            "scale" => (slope = hdr.scl_slope, intercept = hdr.scl_inter),
            "calibration" => (min = hdr.cal_min, max = hdr.cal_max)))

    if hdr.sform_code > 0
        props["spacedirections"] = ntuple(i -> props["header"]["sform"][i], 3)
        ori = spatialorder(props["header"]["sform"])
    else
        props["spacedirections"] = ntuple(i -> props["header"]["qform"][i], 3)
        ori = spatialorder(props["header"]["qform"])
    end
    intent = get(NiftiIntents, hdr.intent_code, eltype(hdr))

    if N < 5
        if intent <: Distribution
            if intent <: StatP0
                props["header"]["intent"] = intent()
            elseif props["header"]["intent"] <: StatP1
                props["header"]["intent"] = intent(hdr.intent_p1)
            elseif props["header"]["intent"] <: StatP2
                props["header"]["intent"] = intent(hdr.intent_p1,intent_p2)
            elseif props["header"]["intent"] <: StatP3
                props["header"]["intent"] = intent(hdr.intent_p1, hdr.intent_p2, hdr.intent_p3)
            end
        else
            props["header"]["intent"] = intent
        end
    end

    axs = Axis{ori[1]}(range(1, step=hdr.pixdim[2], length=hdr.dim[2])*su)
    if N > 1
        axs = (axs, Axis{ori[2]}(range(1, step=hdr.pixdim[3], length=size(hdr, 2))*su))
    end
    if N > 2
        axs = (axs..., Axis{ori[3]}(range(1, step=hdr.pixdim[4], length=size(hdr,3))*su))
    end
    if N > 3
        tu = timeunits(hdr)
        axs = (axs..., Axis{:time}(range(hdr.toffset, step=hdr.pixdim[5], length=size(hdr,4))*tu))
    end
    if N > 4
        for i in 5:N
            axs = (axs..., Axis{Symbol("dim$i")}(range(1, step=hdr.pixdim[i+1], length=size(hdr,i))))
        end
    end
    T = eltype(hdr)
    if T == RGB
        T = Float32
        intent = RGB
        sz = (3, sz...)
    elseif T == RGBA
        T = Float32
        intent = RGBA
        sz = (4, sz...)
    end

    return NiftiSchema{Tuple{sz...},T,intent}(axs, props, needswap)
end

const MetaAxisArray{T,N} = ImageMeta{T,N,AxisArray{T,N,Array{T,N}}}
# NIfTI
function niread(f::AbstractString, sink::Type{A}=MetaAxisArray; mmap::Bool=false) where {A}
    open(f) do s
        if isgz(s)
            return niread(GzipDecompressorStream(s), sink; mmap=mmap)
        else
            return niread(NoopStream(s), sink; mmap=mmap)
        end
    end
end

function niread(s::IO, sink::Type{A}; mmap::Bool=false) where {A}
    schema = NiftiSchema(s)
    niread(s, schema, sink; mmap=mmap)
end

# Analyze
function niread(hdrf::AbstractString, imgf::AbstractString, ::Type{A}=Array; mmap::Bool=false) where {A}
    open(hdrf) do hdrs
        open(imgf) do imgs
            if isgz(hdrs)
                niread(GzipDecompressorStream(hdrs), GzipDecompressorStream(imgs), sink; mmap=mmap)
            else
                niread(NoopStream(hdrs), NoopStream(imgs), sink; mmap=mmap)
            end
        end
    end
end

function niread(hdrs::IO, imgs::IO, sink::Type{A}; mmap::Bool=false) where {A}
    schema = NiftiSchema(hdrs)
    niread(imgs, schema, sink; mmap=mmap)
end

niread(io::IO, schema::NiftiSchema, sink::Type{<:NiftiSchema}; mmap::Bool=false) = schema
niread(io::IO, schema::NiftiSchema, sink::Type{<:AxisArray}; mmap::Bool=false) =
    AxisArray(niread(io, schema, fieldtype(sink, :data); mmap=mmap), schema.axes)
niread(io::IO, schema::NiftiSchema, sink::Type{<:ImageMeta}; mmap::Bool=false) =
    ImageMeta(niread(io, schema, fieldtype(sink, :data); mmap=mmap), schema.props)

# Array
function niread(io::IO, schema::NiftiSchema{S,T}, sink::Type{<:Array};
                mmap::Bool=false) where {S,T}
    if mmap
        _niread(schema, Mmap.mmap(io, Array{T}, (S.parameters...,)))
    else
        _niread(schema, read!(io, Array{T}(undef, (S.parameters...,))))
    end
end

# StaticArray
function niread(io::IO, schema::NiftiSchema{S,T}, sink::Type{<:StaticArray};
                mmap::Bool=false) where {S,T}
    dims = (S.parameters...,)
    if mmap
        _niread(schema, sink{S,T,length(dims),prod(dims)}(Mmap.mmap(io, Array{T}, sims)))
    else
        _niread(schema, read(io, sink{S,T,length(dims),prod(dims)}))
    end
end

function _niread(schema::NiftiSchema{S,T,intent}, A::AbstractArray) where {S,T,intent}
    if schema.needswap
        mappedarray((ntoh, hton), A)
    else
        A
    end
end
