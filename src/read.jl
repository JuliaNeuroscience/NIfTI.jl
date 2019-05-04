const MetaAxisArray{T,N} = ImageMeta{T,N,AxisArray{T,N,Array{T,N}}}
const MetaArray{T,N} = ImageMeta{T,N,Array{T,N}}

function readhdr(io::IO)
    ret = read(io, Int32)
    if ret == Int32(348)
        readhdr1(IOMeta(io, ImageProperties{:NII}()))
    elseif ret == Int32(540)
        readhdr2(IOMeta(io, ImageProperties{:NII}()))
    elseif ret == ntoh(Int32(348))
        readhdr1(IOMeta(SwapStream(io), ImageProperties{:NII}()))
    elseif ret == ntoh(Int32(540))
        readhdr2(IOMeta(SwapStream(io, needswap=true), ImageProperties{:NII}()))
    else
        error("Not a supported NIfTI format")
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
    if isfile(img_file)
        return img_file
    else
        error("NIfTI file is dual file storage, but $img_file does not exist")
    end
end

function gethdr(f::AbstractString)
    hdr_file = replace(f, ".img" => ".hdr")
    if isfile(hdr_file)
        return hdr_file
    else
        error("NIfTI file is dual file storage, but $hdr_file does not exist")
    end
end

niread(f::String, sink::Type{<:AbstractArray}=MetaAxisArray, args...; mode="r", mmap::Bool=false) = open(f, mode) do io
    read(nistreaming(io, f), sink; mmap=mmap)
end

function niread(io::IO, sink::Type{<:AbstractArray}=MetaAxisArray, args...; mode="r", mmap::Bool=false)
    read(nistreaming(io, ), sink; mmap=mmap)
end

nistreaming(f::String; mode::String="r") = open(f, mode) do io
    nistreaming(io, f)
end

function nitype(s::ImageStream)
    if s["header"]["magic"] == NP1_MAGIC
        return "NIfTI-1Single"
    elseif s["header"]["magic"] == NI1_MAGIC
        return "NIfTI-1Double"
    elseif s["header"]["magic"] == NP2_MAGIC
        return "NIfTI-2Single"
    elseif s["header"]["magic"] == NI1_MAGIC
        return "NIfTI-2Double"
    else
        error("Unsupported file type")
    end
end

function nistreaming(io::IO, f::String)
    if isgz(io)
        gzs = gzdopen(io)
        s = readhdr(gzs)
        s["filename"] = [f]
        if nitype(s) == "NIfTI-1Double" || nitype(s) == "NIfTI-2Double"
            append!(s["filenames"], getimg(s["filenames"]))
            close(s.io)
            s.io = gzdopen(s["filenames"][2])
        end
    else
        s = readhdr(io)
        s["filename"] = [f]
        if nitype(s) == "NIfTI-1Double" || nitype(s) == "NIfTI-2Double"
            append!(s["filenames"], getimg(s["filenames"]))
            close(s.io)
            s.io = open(s["filenames"][2])
        end
    end
    return s
end

function readhdr1(s::IOMeta)
    # Uncecessary fields
    skip(s, 35)

    s["header"] = Dict{String,Any}()
    s["header"]["diminfo"] = read(s, Int8)
    N = Int(read(s, Int16))
    sz = ([Int(read(s, Int16)) for i in 1:N]...,)
    skip(s, (7-N)*2)  # skip filler dims

    # intent parameters
    s["header"]["intentparams"] = (float.(read!(s, Vector{Int32}(undef, 3)))...,)
    s["header"]["intent"] = get(NiftiIntents, read(s, Int16), NoIntent)
    T = get(NiftiDatatypes, read(s, Int16), UInt8)

    # skip bitpix
    skip(s, 2)
    s["header"]["slicestart"] = read(s, Int16)

    qfac = read(s, Float32)
    pixdim = (read!(s, Vector{Float32}(undef, N))...,)

    skip(s, (7-N)*4)  # skip filler dims

    s["data_offset"] = Int(read(s, Float32))
    s["header"]["scaleslope"] = read(s, Float32)
    s["header"]["scaleintercept"] = read(s, Float32)
    s["header"]["sliceend"] = read(s, Int16)
    s["header"]["slicecode"] = get(NiftiSliceCodes, read(s, Int8), "Unkown")

    xyzt_units = Int32(read(s, Int8))

    s["calmax"] = read(s, Float32)
    s["calmin"] = read(s, Float32)
    s["header"]["sliceduration"] = read(s, Float32)
    toffset = read(s, Float32)



    skip(s, 8)
    s["description"] = String(read(s, 80))
    s["auxfile"] = String(read(s, 24))

    s["header"]["qformcode"] = get(NiftiXForm, read(s, Int16), :Unkown)
    s["header"]["sformcode"] = get(NiftiXForm, read(s, Int16), :Unkown)
    s["header"]["qform"] = qform(s["header"]["qformcode"],
                                 read!(s, Vector{Float32}(undef, 6))...,
                                 pixdim[1], pixdim[2], pixdim[3], qfac)
    # FIXME
    s["header"]["sform"] = permutedims(SArray{Tuple{4,4},Float32,2,16}(
                                        (read!(s, Vector{Float32}(undef, 12))...,
                                         Float32[0, 0, 0, 1]...)), (2,1))
    s["header"]["intentname"] = String(read(s, 16))
    s["header"]["magic"] = (read(s, 4)...,)

    if s["header"]["sformcode"] ==  :Unkown
        s["spacedirections"] = (Tuple(s["header"]["qform"][1,1:3]),
                                Tuple(s["header"]["qform"][2,1:3]),
                                Tuple(s["header"]["qform"][3,1:3]))
    else
        s["spacedirections"] = (Tuple(s["header"]["sform"][1,1:3]),
                                Tuple(s["header"]["sform"][2,1:3]),
                                Tuple(s["header"]["sform"][3,1:3]))
    end

    s["header"]["extension"] = read(s, NiftiExtension)


    return ImageStream{T}(s, niaxes(sz, xyzt_units, toffset, pixdim, s))
end

function readhdr2(s::IOMeta)
    s["header"] = Dict{String,Any}()
    s["header"]["magic"] = (read(s, 8)...,)
    T = get(NiftiDatatypes, read(s, Int16), UInt8)
    skip(s, 2)  # skip bitpix

    N = read(s, Int64)
    sz = read(s, Vector{Int64}(undef, N))
    skip(s, (7-N)*8)  # skip filler dims

    s["header"]["intentparams"] = (read!(s, Vector{Float64}(undef, 3))...,)

    qfac = read(s, Float64)  # FIXME: qfac is first dimension of otherwise unused pixdim[1]
    pixdim = (read!(s, Vector{Float64}(undef, N))...,)
    skip(s, (7-N)*8)  # skip filler dims

    voxoffset = read(s, Int64)
    s["header"]["scaleslope"] = read(s, Float64)
    s["header"]["scaleintercept"] = read(s, Float64)
    s["calmax"] = read(s, Float64)
    s["calmin"] = read(s, Float64)

    s["header"]["sliceduration"] = read(s, Float64)
    toffset = read(s, Float64)

    s["header"]["slicestart"] = read(s, Int64)
    s["header"]["sliceend"] = read(s, Int64)
    s["description"] = String(read(s, 80))
    s["auxfile"] = String(read(s, 24))
    s["header"]["qformcode"] = get(NiftiXForm, read(s, Int32), :Unkown)
    s["header"]["sformcode"] = get(NiftiXForm, read(s, Int32), :Unkown)

    s["header"]["qform"] = qform(s["header"]["qformcode"], read(s, Vector{Float64}(undef, 6))..., pixdim[1], pixdim[2], pixdim[3], qfac)
    s["header"]["sform"] = transpose(
                                   SMatrix{4,4,Float64,16}(
                                        (read!(s, Vector{Float64}(undef, 12))..., 0.0, 0.0, 0.0, 1.0)))

    s["header"]["slicecode"] = get(NiftiSliceCodes, read(s, Int32), "Unkown")

    xyzt_units = read(s, Int32)
    s["header"]["intent"] = get(NiftiIntents, read(s, Int32), NoIntent)
    s["header"]["intentname"] = String(read(s, 16))
    s["header"]["diminfo"] = read(s, Int8)
    skip(s, 15)

    if s["header"]["sformcode"] ==  :Unkown
        s["spacedirections"] = (Tuple(s["header"]["qform"][1,1:3]),
                                Tuple(s["header"]["qform"][2,1:3]),
                                Tuple(s["header"]["qform"][3,1:3]))
    else
        s["spacedirections"] = (Tuple(s["header"]["sform"][1,1:3]),
                                Tuple(s["header"]["sform"][2,1:3]),
                                Tuple(s["header"]["sform"][3,1:3]))
    end
    s["header"]["extension"] = niread(s, NiftiExtension)

    return ImageStream(s, niaxes(sz, xyzt_units, toffset, pixdim, s))
end
