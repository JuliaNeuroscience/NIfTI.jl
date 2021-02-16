
function read_volume(io::TranscodingStream{GzipDecompressor}, hdr, mmap::Bool)
    T = to_eltype(hdr.datatype)
    dims = to_dimensions(hdr.dim)
    if T === Bool
        A = BitArray{length(dims)}
    else
        A = Array{T,length(dims)}
    end
    if mmap
        close(io)
        error("cannot mmap a gzipped NIfTI file")
    else
        return read!(io, A(undef, dims))
    end
end

function read_volume(io, hdr, mmap::Bool)
    T = to_eltype(hdr.datatype)
    dims = to_dimensions(hdr.dim)
    if T === Bool
        A = BitArray{length(dims)}
    else
        A = Array{T,length(dims)}
    end
    if mmap
        return Mmap.mmap(io, A, dims)
    else
        return read!(io, A(undef, dims))
    end
end

function read_volume(file::AbstractString, mode::AbstractString, hdr, mmap::Bool)
    io = open(file, mode)
    b1 = read(io, UInt8)
    b2 = read(io, UInt8)
    if b1 === 0x1F && b2 === 0x8B
        seek(io, 0)
        return read_volume(GzipDecompressorStream(io), hdr, mmap)
    else
        seek(io, 0)
        return read_volume(io, hdr, mmap)
    end
end


write_volume(io, x::AbstractArray{T}) where {T} = write(io, x)
write_volume(io, x::AbstractArray{Bool}) = write_volume(io, BitArray(x))
write_volume(io, x::BitArray) = write(io, x)

