# Avoid method ambiguity
write(io::Base64.Base64EncodePipe, vol::NIVolume{UInt8,1}) = invoke(write, (IO, NIVolume{UInt8,1}), io, vol)

# Write a NIfTI file
function write(io::IO, vol::NIVolume)
    write(io, niupdate(vol).header)
    if isempty(vol.extensions)
        write(io, Int32(0))
    else
        for ex in vol.extensions
            sz = esize(ex)
            write(io, Int32(sz))
            write(io, Int32(ex.ecode))
            write(io, ex.edata)
            write(io, zeros(UInt8, sz - length(ex.edata)))
        end
    end
    write(io, vol.raw)
end

# Convenience function to write a NIfTI file given a path
function niwrite(path::AbstractString, vol::NIVolume)
    if split(path,".")[end] == "gz"
        iogz = gzopen(path, "w9")
        write(iogz, vol)
        close(iogz)
    else
        io = open(path, "w")
        write(io, vol)
        close(io)
    end
end

# Read header from a NIfTI file
function read_header(io::IO)
    header, swapped = read(io, NIfTI1Header)
    if header.sizeof_hdr != SIZEOF_HDR
        error("This is not a NIfTI-1 file")
    end
    header, swapped
end

# Read extension fields following NIfTI header
function read_extensions(io::IO, header::NIfTI1Header)
    if eof(io)
        return NIfTI1Extension[]
    end

    extension = read!(io, Array{UInt8}(undef, 4))
    if extension[1] != 1
        return NIfTI1Extension[]
    end

    extensions = NIfTI1Extension[]
    should_bswap = header.dim[1] > 7
    while header.magic == NP1_MAGIC ? position(io) < header.vox_offset : !eof(io)
        esize = read(io, Int32)
        ecode = read(io, Int32)

        if should_bswap
            esize = bswap(esize)
            ecode = bswap(ecode)
        end

        push!(extensions, NIfTI1Extension(ecode, read!(io, Array{UInt8}(undef, esize-8))))
    end
    extensions
end

# Look for a gzip header in an IOStream
function isgz(io::IO)
    ret = read(io, UInt8) == 0x1F && read(io, UInt8) == 0x8B
    seek(io, 0)
    ret
end

# Read a NIfTI file. The optional mmap argument determines whether the contents are read in full
# (if false) or mmaped from the disk (if true).
"""

"""
function niread(file::AbstractString; mmap::Bool=false)
    file_io = open(file, "r")
    header_gzipped = isgz(file_io)
    header_io = header_gzipped ? gzdopen(file_io) : file_io
    header, swapped = read_header(header_io)
    extensions = read_extensions(header_io, header)
    dims = convert(Tuple{Vararg{Int}}, header.dim[2:header.dim[1]+1])

    if !haskey(NIfTI_DT_BITSTYPES, header.datatype)
        error("data type $(header.datatype) not yet supported")
    end
    dtype = NIfTI_DT_BITSTYPES[header.datatype]

    local volume
    if header.magic == NP1_MAGIC
        if mmap
            if header_gzipped
                close(header_io)
                close(file_io)
                error("cannot mmap a gzipped NIfTI file")
            else
                volume = Mmap.mmap(header_io, Array{dtype,length(dims)}, dims, Int(header.vox_offset))
            end
        else
            seek(header_io, Int(header.vox_offset))
            volume = read!(header_io, Array{dtype}(undef, dims))
            if !eof(header_io)
                warn("file size does not match length of data; some data may be ignored")
            end
            close(header_io)
            !header_gzipped || close(file_io)
        end
    elseif header.magic == NI1_MAGIC
        close(header_io)
        !header_gzipped || close(file_io)

        volume_name = replace(file, r"\.\w+(\.gz)?$" => "")*".img"
        if !isfile(volume_name)
            if isfile(volume_name*".gz")
                volume_name *= ".gz"
            else
                error("NIfTI file is dual file storage, but $volume_name does not exist")
            end
        end

        volume_io = open(volume_name, "r")
        volume_gzipped = isgz(volume_io)
        if mmap
            if volume_gzipped
                close(volume_io)
                error("cannot mmap a gzipped NIfTI file")
            else
                volume = Mmap.mmap(volume_io, Array{dtype,length(dims)}, dims)
            end
        else
            if volume_gzipped
                volume_gz_io = gzdopen(volume_io)
                volume = read!(volume_gz_io, Array{dtype}(undef, dims))
                close(volume_gz_io)
            else
                volume = read!(volume_io, Array{dtype}(undef, dims))
            end
            close(volume_io)
        end
    end
    if swapped && sizeof(eltype(volume)) > 1
        volume = mappedarray((ntoh, hton), volume)
    end

    return NIVolume(header, extensions, volume)
end

# Allow file to be indexed like an array, but with indices yielding scaled data
@inline getindex(f::NIVolume{T}, idx::Vararg{Int}) where {T} =
    getindex(f.raw, idx...,) * f.header.scl_slope + f.header.scl_inter

add1(x::Union{AbstractArray{T},T}) where {T<:Integer} = x + 1
add1(::Colon) = Colon()
@inline vox(f::NIVolume, args...,) = getindex(f, map(add1, args)...,)

size(f::NIVolume) = size(f.raw)
size(f::NIVolume, d) = size(f.raw, d)
ndims(f::NIVolume) = ndims(f.raw)
length(f::NIVolume) = length(f.raw)
lastindex(f::NIVolume) = lastindex(f.raw)
