const MetaAxisArray{T,N} = ImageMeta{T,N,AxisArray{T,N,Array{T,N}}}
const MetaArray{T,N} = ImageMeta{T,N,Array{T,N}}

function isgz(io::IO)
    gzbits = read(io, 2)
    ret = gzbits == [0x1F,0x8B]
    seek(io, 0)
    ret
end

function isgz(f::AbstractString)
    split(f, '.')[end] == "gz"
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

function load(f::File{DataFormat{:NII}}, sink::Type{<:AbstractArray}=MetaAxisArray;
              mode::String="r", mmap::Bool=false, grow::Bool=true)
    open(f, mode) do s
        load(s, sink, mmap=mmap, grow=grow)
    end
end

function load(s::Stream{DataFormat{:NII}}, sink::Type{<:AbstractArray}=MetaAxisArray;
              mmap::Bool=false, grow::Bool=true)
    read(_metadata(stream(s), filename(s)), sink, mmap=mmap, grow=grow)
end

function metadata(s::File{DataFormat{:NII}})
    open(f) do s
        ret = metadata(s)
    end
    return ret
end

metadata(s::Stream{DataFormat{:NII}}) =  getinfo(_metadata(stream(s), filename(s)))

function _metadata(io::IO, f::String)
    if isgz(io)
        gzs = gzdopen(io)
        s = readhdr(gzs)
        s["filename"] = f
        if nitype(s) == "NIfTI-1Double" || nitype(s) == "NIfTI-2Double"
            close(s.io)
            s.io = gzdopen(open(getimg(s["filename"])))
        end
    else
        s = readhdr(io)
        s["filename"] = f
        if nitype(s) == "NIfTI-1Double" || nitype(s) == "NIfTI-2Double"
            close(s.io)
            s.io = open(getimg(s["filename"]))
        end
    end
    return s
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

function readhdr(io::IO)
    ret = read(io, Int32)
    if ret == Int32(348)
        readhdr1(io, ImageProperties{format"NII"}())
    elseif ret == Int32(540)
        readhdr2(io, ImageProperties{format"NII"}())
    elseif ret == ntoh(Int32(348))
        readhdr1(SwapStream(io, needswap=true), ImageProperties{format"NII"}())
    elseif ret == ntoh(Int32(540))
        readhdr2(SwapStream(io, needswap=true), ImageProperties{format"NII"}())
    else
        error("Not a supported NIfTI format")
    end
end

function readhdr1(s::IO, p::ImageProperties)
    # Uncecessary fields
    skip(s, 35)

    p["header"] = ImageProperties{:header}()
    p["header"]["diminfo"] = read(s, Int8)
    N = Int(read(s, Int16))
    sz = ([Int(read(s, Int16)) for i in 1:N]...,)
    skip(s, (7-N)*2)  # skip filler dims

    # intent parameters
    p["header"]["intentparams"] = (float.(read!(s, Vector{Int32}(undef, 3)))...,)
    p["header"]["intent"] = get(NiftiIntents, read(s, Int16), NoIntent)
    T = get(NiftiDatatypes, read(s, Int16), UInt8)

    # skip bitpix
    skip(s, 2)
    p["header"]["slicestart"] = Int(read(s, Int16)) + 1  # to 1 based indexing

    qfac = read(s, Float32)
    pixdim = (read!(s, Vector{Float32}(undef, N))...,)

    skip(s, (7-N)*4)  # skip filler dims

    p["data_offset"] = Int(read(s, Float32))
    p["header"]["scaleslope"] = Float64(read(s, Float32))
    p["header"]["scaleintercept"] = Float64(read(s, Float32))
    p["header"]["sliceend"] = Int(read(s, Int16) + 1)  # to 1 based indexing
    p["header"]["slicecode"] = get(NiftiSliceCodes, read(s, Int8), "Unkown")

    xyzt_units = Int32(read(s, Int8))

    p["calmax"] = read(s, Float32)
    p["calmin"] = read(s, Float32)
    p["header"]["sliceduration"] = read(s, Float32)
    toffset = read(s, Float32)

    skip(s, 8)
    p["description"] = String(read(s, 80))
    p["auxfiles"] = [String(read(s, 24))]

    p["header"]["qformcode"] = get(NiftiXForm, read(s, Int16), :Unkown)
    p["header"]["sformcode"] = get(NiftiXForm, read(s, Int16), :Unkown)
    qb = read(s, Float32)
    qc = read(s, Float32)
    qd = read(s, Float32)
    qx = read(s, Float32)
    qy = read(s, Float32)
    qz = read(s, Float32)

    # FIXME
    p["header"]["sform"] = permutedims(SArray{Tuple{4,4},Float32,2,16}(
                                        (read!(s, Vector{Float32}(undef, 12))...,
                                         Float32[0, 0, 0, 1]...)), (2,1))
    p["header"]["intentname"] = String(read(s, 16))
    p["header"]["magic"] = (read(s, 4)...,)

    if p["header"]["sformcode"] == :Unkown
        qf = quat2mat(qb, qc, qd, qx, qy, qz, pixdim[1:min(N,3)]..., zeros(Float64, 3-min(N,3))..., qfac)
        p["spacedirections"] = (Tuple(qf[1,1:3]), Tuple(qf[2,1:3]), Tuple(qf[3,1:3]))
    else
        p["spacedirections"] = (Tuple(p["header"]["sform"][1,1:3]),
                                Tuple(p["header"]["sform"][2,1:3]),
                                Tuple(p["header"]["sform"][3,1:3]))
    end

    p["header"]["extension"] = read(s, p, NiftiExtension)

    return ImageStream{T}(s, niaxes(sz, xyzt_units, toffset, qx, qy, qz, pixdim, p), p)
end

function readhdr2(s::IO, p::ImageProperties)
    p["header"] = ImageProperties{:header}()
    p["header"]["magic"] = (read(s, 8)...,)
    T = get(NiftiDatatypes, read(s, Int16), UInt8)
    skip(s, 2)  # skip bitpix

    N = read(s, Int64)
    sz = Tuple(read!(s, Vector{Int64}(undef, N)))
    skip(s, (7-N)*8)  # skip filler dims

    p["header"]["intentparams"] = (read!(s, Vector{Float64}(undef, 3))...,)

    qfac = read(s, Float64)
    pixdim = (read!(s, Vector{Float64}(undef, N))...,)
    skip(s, (7-N)*8)  # skip filler dims

    p["data_offset"] = read(s, Int64)
    p["header"]["scaleslope"] = read(s, Float64)
    p["header"]["scaleintercept"] = read(s, Float64)
    p["calmax"] = read(s, Float64)
    p["calmin"] = read(s, Float64)

    p["header"]["sliceduration"] = read(s, Float64)
    toffset = read(s, Float64)

    p["header"]["slicestart"] = read(s, Int64) + 1  # to 1 based indexing
    p["header"]["sliceend"] = read(s, Int64) + 1
    p["description"] = String(read(s, 80))
    p["auxfiles"] = [String(read(s, 24))]
    p["header"]["qformcode"] = get(NiftiXForm, read(s, Int32), :Unkown)
    p["header"]["sformcode"] = get(NiftiXForm, read(s, Int32), :Unkown)

    qb = read(s, Float64)
    qc = read(s, Float64)
    qd = read(s, Float64)
    qx = read(s, Float64)
    qy = read(s, Float64)
    qz = read(s, Float64)

    p["header"]["sform"] = transpose(
                                   SMatrix{4,4,Float64,16}(
                                        (read!(s, Vector{Float64}(undef, 12))..., 0.0, 0.0, 0.0, 1.0)))

    p["header"]["slicecode"] = get(NiftiSliceCodes, read(s, Int32), "Unkown")

    xyzt_units = read(s, Int32)
    p["header"]["intent"] = get(NiftiIntents, read(s, Int32), NoIntent)
    p["header"]["intentname"] = String(read(s, 16))
    p["header"]["diminfo"] = read(s, Int8)
    skip(s, 15)

    if p["header"]["sformcode"] ==  :Unkown
        qf = quat2mat(qb, qc, qd, qx, qy, qz, pixdim[1:min(N,3)]..., zeros(Float64, 3-min(N,3))..., qfac)
        p["spacedirections"] = (Tuple(qf[1,1:3]), Tuple(qf[2,1:3]), Tuple(qf[3,1:3]))
    else
        p["spacedirections"] = (Tuple(p["header"]["sform"][1,1:3]),
                                Tuple(p["header"]["sform"][2,1:3]),
                                Tuple(p["header"]["sform"][3,1:3]))
    end
    p["header"]["extension"] = read(s, p, NiftiExtension)

    return ImageStream{T}(s, niaxes(sz, xyzt_units, toffset, qx, qy, qz, pixdim, p), p)
end
