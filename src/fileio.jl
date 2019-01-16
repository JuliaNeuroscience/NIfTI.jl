function detectnii(io::IO)
    isgzipped = isgz(io)
    hdrio = isgzipped ? gzdopen(io) : io
    version, swap = checkfile(hdrio)
    if version == 1
        magic = seek(hdrio, number_to_magic)
        return magic == NP1_MAGIC ? true : false
    elseif version == 2
        magic = read(hdrio, Array{UInt8}(undef,8))
        return magic == NP1_MAGIC ? true : false
    else
        return false
    end
end

detecthdr(io::IO) = detectnii(io)

function detectimg(io::IO)
    img_file = filename(io)
    hdr_file = gethdr(img_file)
    open(hdr_file, "r") do io
        detectnii(io)
    end
end

add_format(format"NII", detectnii, [".nii"], [:NIfTI])
add_format(format"NIIGZ", (), [".nii.gz"], [:NIfTI])

# LoadError: format FileIO.DataFormat{:HDR} is already registered
#add_format(format"HDR", detecthdr, [".hdr"], [:NIfTI])
#add_format(format"HDRGZ", (), [".hdr.gz"], [:NIfTI])

add_format(format"IMG", detectimg, [".img"], [:NIfTI])
add_format(format"IMGGZ", (), [".img.gz"], [:NIfTI])

"""
using NIfTI, AxisArrays, ImageMetadata

const MetaAxisArray{T,N} = ImageMeta{T,N,AxisArray{T,N,Array{T,N}}}
f = "test/data/example4d.nii"
io = open(f)

tmp = load(io, MetaAxisArray)
"""
function load(f::File{format"NII"}, sinkType=MetaAxisArray; mmap::Bool=false)
    open(f) do s
        load(stream(s), sink; mmap=mmap)
    end
end

function load(f::File{format"NIIGZ"}, sink::Type=MetaAxisArray)
    GzipDecompressorStream(f) do s
        load(stream(s), sink; mmap=mmap)
    end
end

function load(hdrf::File{format"HDR"}, sink::Type=MetaAxisArray; mmap::Bool=false)
    imgf = getimg(hdrf)
    open(hdrf) do hdrio
        open(imgf) do imgio
            load(stream(hdrio), stream(imgio), sink; mmap=mmap)
        end
    end
end

function load(hdrf::File{format"HDRGZ"}, sink=NiftiStream)
    imgf = getimg(hdrf)
    GzipDecompressorStream(hdrf) do hdrio
        GzipDecompressorStream(imgf) do imgiio
            load(stream(hdrio), stream(imgio), sink; mmap=mmap)
        end
    end
end

function load(imgff::File{format"IMG"}, sink=NiftiStream; mmap::Bool=false)
    hdrf = gethdr(imgf)
    open(hdrf) do hdrio
        open(imgf) do imgio
            load(stream(hdrio), stream(imgio), sink; mmap=mmap)
        end
    end
end

function load(imgf::File{format"IMGGZ"}, sink)
    hdrf = gethdr(imgf)
    GzipDecompressorStream(hdrf) do hdrio
        GzipDecompressorStream(imgf) do imgio
            load(stream(hdrio), stream(imgio), sink; mmap=mmap)
        end
    end
end
