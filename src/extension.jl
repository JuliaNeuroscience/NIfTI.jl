
mutable struct NiftiExtension
    ecode::Int32
    edata::Vector{UInt8}
end

"""
    setoffset!(hdr::NiftiHeader, ext::NiftiExtension)

Set the offset to volume data accounting for extension size
"""
setoffset!(hdr::NiftiHeader, ext::NiftiExtension) =
    hdr.vox_offset = esize(ext) + hdr.sizeof_hdr

# Calculates the size of a NIfTI extension
esize(ex::NiftiExtension) = 8 + ceil(Int, length(ex.edata)/16)*16

# TODO figure this garbage out
# https://www.nitrc.org/forum/forum.php?thread_id=4380&forum_id=1955
# https://www.nitrc.org/forum/attachment.php?attachid=341&group_id=454&forum_id=1955
function read_extensions(io::IO, hdr::NiftiHeader)
    if eof(io)
        return NiftiExtension[]
    end

    extension = read!(io, Array{UInt8}(undef, 4))
    if extension[1] != 1
        return NiftiExtension[]
    end

    extensions = NiftiExtension[]
    should_bswap = hdr.dim[1] > 7
    while hdr.magic == NP1_MAGIC ? position(io) < hdr.vox_offset : !eof(io)
        esize = read(io, Int32)
        ecode = read(io, Int32)

        if should_bswap
            esize = bswap(esize)
            ecode = bswap(ecode)
        end

        push!(extensions, NiftiExtension(ecode, read!(io, Array{UInt8}(undef, esize-8))))
    end
    extensions
end

# Write a NIfTI file
function write_extension(io::IO, ext::NiftiExtension)
    if isempty(ext)
        write(io, Int32(0))
    else
        for ex in ext
            sz = esize(ex)
            write(io, Int32(sz))
            write(io, Int32(ex.ecode))
            write(io, ex.edata)
            write(io, zeros(UInt8, sz - length(ex.edata)))
        end
    end
end
