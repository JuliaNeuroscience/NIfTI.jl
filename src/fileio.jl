function load(f::File{format"NII"}, sink::Type{<:AbstractArray}=MetaAxisArray,
    args...; mode="r", mmap::Bool=false)
    open(f, mode) do s
        load(s, sink, args...; mmap=mmap)
    end
end

function load(s::Stream{format"NII"}, sink::Type{<:AbstractArray}=MetaAxisArray, args...; mmap::Bool=false)
    read(loadstreaming(s, args...), sink; mmap=mmap)
end

function loadstreaming(f::File{format"NII"}, args...; mode="r")
    open(f, mode) do io
        loadstreaming(io, filename(f))
    end
end

loadstreaming(s::Stream{format"NII"}, args...) = nistreaming(stream(s), filename(s))

### This ↓ needs to go in the FileIO registry
function detectnii(io::IO)
    ret = read(io, Int32)
    if ret == Int32(348)
        seek(io, 344)
        return read(io, 4) == UInt8[0x6e,0x2b,0x31,0x00] |
                read(io, 4) == UInt8[0x6e,0x69,0x31,0x00]
    elseif ret == Int32(540)
        return read(io, 8) == UInt8[0x6e,0x2b,0x32,0x00,0x0d,0x0a,0x1a,0x0a] |
                read(io, 8) == UInt8[0x6e,0x69,0x32,0x00,0x0d,0x0a,0x1a,0x0a]
    elseif ret == ntoh(Int32(348))
        return read(io, 4) == UInt8[0x6e,0x2b,0x31,0x00] |
                read(io, 4) == UInt8[0x6e,0x69,0x31,0x00]
    elseif ret == ntoh(Int32(540))
        return read(io, 8) == UInt8[0x6e,0x2b,0x32,0x00,0x0d,0x0a,0x1a,0x0a] |
                read(io, 8) == UInt8[0x6e,0x69,0x32,0x00,0x0d,0x0a,0x1a,0x0a]
    else
        false
    end
end

function detectimg(io::IO)
    img_file = filename(io)
    hdr_file = gethdr(img_file)
    open(hdr_file, "r") do io
        detectnii(io)
    end
end

add_format(format"NII", detectnii, [".nii", ".nii.gz"], [:NIfTI])
# Analyze format
add_format(format"ANZ", detectimg, [".img"], [:NIfTI])
#add_format(format"ANZ", detectnii, [".hdr"], [:NIfTI])


"""
using NIfTI, AxisArrays, ImageMetadata

const MetaAxisArray{T,N} = ImageMeta{T,N,AxisArray{T,N,Array{T,N}}}
f = "test/data/example4d.nii"
io = open(f)

tmp = load(io, MetaAxisArray)
"""

