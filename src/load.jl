"""
    load(::DataFormat{:NII}, args...)

# Arguments

* `arraytype::Type`
* `mode::String`
* `mmap::Bool`
* `grow::Bool`
"""
function load(
    f::Formatted{format"NII"},
    ::Type{A}=Array,
    mode::String="r",
    mmap::Bool=false,
    grow::Bool=true) where {A<: AbstractArray}
    read(loadstreaming(f, guessintent(split(filename(f), '.'), A), mode); mmap=mmap, grow=grow)
end


"""
    loadstreaming(::DataFormat{:NII}, args...)
"""
function loadstreaming(s::File{DataFormat{:NII}}, ::Type{A}, mode::String="r") where {A<:AbstractArray}
    loadstreaming(open(s, mode))
end

function loadstreaming(s::Stream{DataFormat{:NII},GZip.GZipStream})
    m = loadmeta(s)
    if magicbytes(m) == NI1_MAGIC || magicbytes(m) == NI2_MAGIC
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
        if magicbytes(m) == NI1_MAGIC || magicbytes(m) == NI2_MAGIC
            close(s)
            return ImageStream(open(getimg(filename(s))), m)
        else
            return ImageStream(stream(s), m)
        end
    end
end


"""
    metadata(::DataFormat{:NII}) -> ImageInfo
"""
function metadata(s::File{DataFormat{:NII}})
    open(s) do s
        metadata(s)
    end
end

function metadata(s::Stream{DataFormat{:NII},IOType}) where IOType
    if split(filname(s), '.')[end] == ".gz"
        return loadmeta(Stream(DataFormat{:NII}, gzdopen(stream(s)), filename(s)))
    else
        return loadmeta(s)
    end
end

metadata(s::Stream{DataFormat{:NII},GZip.GZipStream}) = loadmeta(s)

function loadmeta(s::Stream{DataFormat{:NII}})
    m = loadmeta(stream(s))
    srcfile!(m, filename(s))
    return m
end

function loadmeta(io::IO)
    ret = read(io, Int32)
    if ret == Int32(348)
        readhdr1(io, Dict{String,Any}())
    elseif ret == Int32(540)
        readhdr2(io, Dict{String,Any}())
    elseif ret == ntoh(Int32(348))
        readhdr1(SwapStream(io, needswap=true), Dict{String,Any}())
    elseif ret == ntoh(Int32(540))
        readhdr2(SwapStream(io, needswap=true), Dict{String,Any}())
    else
        error("Not a supported NIfTI format")
    end
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


