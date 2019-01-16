mutable struct NiftiExtension
    ecode::Int32
    edata::Vector{UInt8}
end

# Calculates the size of a NIfTI extension
esize(ex::NiftiExtension) = 8 + ceil(Int, length(ex.edata)/16)*16

# TODO figure this garbage out
# https://www.nitrc.org/forum/forum.php?thread_id=4380&forum_id=1955
# https://www.nitrc.org/forum/attachment.php?attachid=341&group_id=454&forum_id=1955

# makes empty extension
NiftiExtension() = NiftiExtension(Int32(0), zeros(Int8,4))
NiftiExtension(img::ImageMeta) = @get img "extension" NiftiExtension()
NiftiExtension(img::AbstractArray) = NiftiExtension

function niread(io::IO, hdr::NiftiHeader, needswap::Bool, ::Type{NiftiExtension})
    if eof(io)
        return NiftiExtension[]
    end

    extension = read!(io, Array{UInt8}(undef, 4))
    if extension[1] != 1
        return NiftiExtension[]
    else

        esize = read(io, Int32)
        ecode = read(io, Int32)

        if needswap
            esize = bswap(esize)
            ecode = bswap(ecode)
        end

        return NiftiExtension(ecode, read!(io, Array{UInt8}(undef, esize-8)))
    end
end

# TODO
# io is at write position for writing extension. this is only an issue for
# cases of writing dynamically writing extension or volume data
function write(io::IO, ext::NiftiExtension)
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
