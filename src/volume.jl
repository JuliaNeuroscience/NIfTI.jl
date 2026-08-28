
function read_volume(
    io::TranscodingStream{GzipDecompressor},
    ::Type{T},
    dim::Tuple{Vararg{Any,N}},
    map::Bool,
) where {T,N}

    if T === Bool
        A = BitArray{N}
    else
        A = Array{T,N}
    end
    if map
        close(io)
        error("cannot mmap a gzipped NIfTI file")
    else
        return read!(io, A(undef, dim))
    end
end

function read_volume(io, ::Type{T}, dim::Tuple{Vararg{Any,N}}, map::Bool) where {T,N}
    if T === Bool
        A = BitArray{N}
    else
        A = Array{T,N}
    end
    if map
        return Mmap.mmap(io, A, dim)
    else
        return read!(io, A(undef, dim))
    end
end

niopen(file::AbstractString, mode) = niopen(open(file, mode))
function niopen(io)
    if isgz(io)
        return GzipDecompressorStream(io)
    else
        return io
    end
end

