const MetaAxisArray{T,N} = ImageMeta{T,N,AxisArray{T,N,Array{T,N}}}
const MetaArray{T,N} = ImageMeta{T,N,Array{T,N}}

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

function load(f::Formatted{format"NII"}, sink=MetaAxisArray;
              mode::String="r", mmap::Bool=false, grow::Bool=true)
        read(loadstreaming(f, mode=mode), sink; mmap=mmap, grow=grow)
end

#=
function load(f::File{DataFormat{:NII}}, sink::Type{<:AbstractArray}=MetaAxisArray;
              mode::String="r", mmap::Bool=false, grow::Bool=true)
    loadstreaming() do s
    open(f, mode) do s
        load(s, sink, mmap=mmap, grow=grow)
    end
end
=#

loadstreaming(s::File{DataFormat{:NII}}; mode::String="r") = loadstreaming(open(s, mode))

function loadstreaming(s::Stream{DataFormat{:NII},GZip.GZipStream})
    m = loadmeta(s)
    if nitype(m) == "NIfTI-1Double" || nitype(m) == "NIfTI-2Double"
        ret = ImageStream(gzdopen(open(getimg(filename(s)))), m)
        close(s)
    else
        ret = ImageStream(stream(s), m)
    end
    return ret
end

function loadstreaming(s::Stream{DataFormat{:NII},IOType}) where IOType
    if file_extension(s) == ".gz"
        return loadstreaming(Stream(DataFormat{:NII}, gzdopen(stream(s)), filename(s)))
    else
        m = loadmeta(s)
        if nitype(m) == "NIfTI-1Double" || nitype(m) == "NIfTI-2Double"
            close(s)
            return ImageStream(open(getimg(filename(s))), m)
        else
            return ImageStream(stream(s), m)
        end
    end
end

function metadata(s::File{DataFormat{:NII}})
    open(s) do s
        metadata(s)
    end
end

function nitype(s::ImageInfo)
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

function metadata(s::Stream{DataFormat{:NII},IOType}) where IOType
    if file_extension(s) == ".gz"
        return loadmeta(Stream(DataFormat{:NII}, gzdopen(stream(s)), filename(s)))
    else
        return loadmeta(s)
    end
end

metadata(s::Stream{DataFormat{:NII},GZip.GZipStream}) = loadmeta(s)

function loadmeta(s::Stream{DataFormat{:NII}})
    m = loadmeta(stream(s))
    m["filename"] = filename(s)
    return m
end

function loadmeta(io::IO)
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
    p["header"]["intentparams"] = Tuple(float.(read!(s, Vector{Int32}(undef, 3))))::NTuple{3,Float64}
    p["header"]["intent"] = get(NiftiIntents, read(s, Int16), NoIntent)
    T = get(NiftiDatatypes, read(s, Int16), UInt8)

    # skip bitpix
    skip(s, 2)
    p["header"]["slicestart"] = Int(read(s, Int16)) + 1  # to 1 based indexing

    p["header"]["qfac"] = read(s, Float32)
    pixdim = (read!(s, Vector{Float32}(undef, N))...,)

    skip(s, (7-N)*4)  # skip filler dims

    p["data_offset"] = Int(read(s, Float32))
    p["header"]["scaleslope"] = Float64(read(s, Float32))
    p["header"]["scaleintercept"] = Float64(read(s, Float32))
    p["header"]["sliceend"] = Int(read(s, Int16) + 1)  # to 1 based indexing
    p["header"]["slicecode"] = get(NiftiSliceCodes, read(s, Int8), "Unkown")::String

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

    p["header"]["quaternb"] = read(s, Float32)
    p["header"]["quaternc"] = read(s, Float32)
    p["header"]["quaternd"] = read(s, Float32)
    p["header"]["qoffsetx"] = read(s, Float32)
    p["header"]["qoffsety"] = read(s, Float32)
    p["header"]["qoffsetz"] = read(s, Float32)
    # FIXME
    p["header"]["sform"] = permutedims(SArray{Tuple{4,4},Float32,2,16}(
                                        (read!(s, Vector{Float32}(undef, 12))...,
                                         Float32[0, 0, 0, 1]...)), (2,1))
    p["header"]["intentname"] = String(read(s, 16))
    p["header"]["magic"] = Tuple(read(s, 4))::NTuple{4,UInt8}

    if p["header"]["sformcode"] == :Unkown
        qf = quat2mat(p["header"]["quaternb"],
                      p["header"]["quaternc"],
                      p["header"]["quaternc"],
                      p["header"]["qoffsetx"],
                      p["header"]["qoffsety"],
                      p["header"]["qoffsetz"], pixdim[1:min(N,3)]..., zeros(Float64, 3-min(N,3))..., p["header"]["qfac"])

        p["spacedirections"] = (Tuple(qf[1,1:3]), Tuple(qf[2,1:3]), Tuple(qf[3,1:3]))
    else
        p["spacedirections"] = (Tuple(p["header"]["sform"][1,1:3]),
                                Tuple(p["header"]["sform"][2,1:3]),
                                Tuple(p["header"]["sform"][3,1:3]))
    end

    p["header"]["extension"] = read(s, p, NiftiExtension)

    return ImageInfo{T}(niaxes(sz, xyzt_units, toffset, p["header"]["qoffsetx"], p["header"]["qoffsety"], p["header"]["qoffsetz"], pixdim, p), p)
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

    p["header"]["qfac"] = read(s, Float64)
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

    p["header"]["quaternb"] = read(s, Float64)
    p["header"]["quaternc"] = read(s, Float64)
    p["header"]["quaternd"] = read(s, Float64)
    p["header"]["qoffsetx"] = read(s, Float64)
    p["header"]["qoffsety"] = read(s, Float64)
    p["header"]["qoffsetz"] = read(s, Float64)

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
        qf = quat2mat(p["header"]["quaternb"],
                      p["header"]["quaternc"],
                      p["header"]["quaternc"],
                      p["header"]["qoffsetx"],
                      p["header"]["qoffsety"],
                      p["header"]["qoffsetz"], pixdim[1:min(N,3)]..., zeros(Float64, 3-min(N,3))..., p["header"]["qfac"])
        p["spacedirections"] = (Tuple(qf[1,1:3]), Tuple(qf[2,1:3]), Tuple(qf[3,1:3]))
    else
        p["spacedirections"] = (Tuple(p["header"]["sform"][1,1:3]),
                                Tuple(p["header"]["sform"][2,1:3]),
                                Tuple(p["header"]["sform"][3,1:3]))
    end
    p["header"]["extension"] = read(s, p, NiftiExtension)

    return ImageInfo{T}(niaxes(sz, xyzt_units, toffset, p["header"]["qoffsetx"], p["header"]["qoffsety"], p["header"]["qoffsetz"], pixdim, p), p)
end
