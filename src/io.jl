# Look for a gzip header in an IOStream
function isgz(io::IO)
    ret = read(io, UInt8) == 0x1F && read(io, UInt8) == 0x8B
    seek(io, 0)
    ret
end

function getimg(f::AbstractString)
    img_file = replace(f, ".hdr" => ".img")
    if isfile(out)
        img_file
    else
        error("NIfTI file is dual file storage, but $img_file does not exist")
    end
end

function gethdr(f::AbstractString)
    hdr_file = replace(f, ".img" => ".hdr")
    if isfile(out)
        img_file
    else
        error("NIfTI file is dual file storage, but $hdr_file does not exist")
    end
end

# load
function read_volume(io::IO, hdr::NiftiHeader, needswap::Bool=false;
                     mmap::Bool=false)
    seek(io, Int(hdr.vox_offset))
    dims = convert(Tuple{Vararg{Int}}, hdr.dim[2:hdr.dim[1]+1])

    if !haskey(NIFTI_DT_BITSTYPES, hdr.datatype)
        error("data type $(hdr.datatype) not yet supported")
    end

    dtype = NIFTI_DT_BITSTYPES[hdr.datatype]
    if mmap
        vol = Mmap.mmap(io, Array{dtype,length(dims)}, dims)
    else
        vol = read!(io, Array{dtype}(undef, dims))
    end
    if needswap && sizeof(eltype(vol)) > 1
        vol = mappedarray((ntoh, hton), vol)
    end

    # scale if slope isn't zero
    if hdr.scl_slope != AbstractFloat(0) # && dtype != RGB24  # TODO finalize RGB types
        vol = vol.*hdr.scl_slope.+hdr.scl_inter
    end
    vol
end

# TODO
#function write_volume(io::IO, vol::AbstractArray)
#end

function load(f::File{format"NII"}; mode="r", mmap::Bool=false,
              returntype::Type=ImageMeta)
    open(f, mode) do s
        load(s; mmap=mmap, returntype=returntype)
    end
end

"""
* `returntype::Type=ImageMeta`: allows specification of what format the NIfTI
file is read into Julia as:
    * ImageMeta: includes header, extension, and volume data
    * NiftiHeader: only loads subtypes of NiftiHeader type
"""
function load(s::Stream{format"NII"}; mmap::Bool=false, returntype::Type=ImageMeta)
    version, needswap = checkfile(stream(s))
    if returntype == ImageMeta
        header = ifelse(version == 1,
                        read(stream(s), Nifti1Header, needswap),
                        read(stream(s), Nifti2Header, needswap))
        extension = read_extension(stream(s), header, needswap)
        volume = read_volume(stream(s), header, needswap; mmap=mmap)
        nii2img(volume, extension, header)
    elseif returntype == NiftiHeader
        ifelse(version == 1,
               read(s.io, Nifti1Header, needswap),
               read(s.io, Nifti2Header, needswap))
    else
        error("$returntype is not currently a supported NIfTI return type.")
    end
end

function load(hdrf::File{format"HDR"}; mode="r", mmap::Bool=false,
              returntype::String="image")
    imgf = getimg(hdrf)
    open(hdrf, mode) do hdrs
        open(imgf) do imgs
            isgzipped = isgz(hdrs)
            if isgzipped
                if mmap
                    error("Cannot mmap a gzipped NIfTI file.")
                end
                gzdopen(hdrs) do gzhdrs
                    gzdopen(imgs) do gzimgs
                        out = load(gzhdrs, gzimgs; mmap=mmap,
                                   returntype=returntype)
                    end
                end
            else
                open(hdrs, imgs; mmap=mmap, returntype=returntype)
            end
        end
    end
end

function load(imgf::File{format"IMG"}; mode="r", mmap::Bool=false,
              returntype::Type=ImageMeta)
    hdrf = gethdr(imgf)
    open(hdrf, mode) do hdrs
        open(imgf) do imgs
            isgzipped = isgz(hdrs)
            if isgzipped
                if mmap
                    error("Cannot mmap a gzipped NIfTI file.")
                end
                gzdopen(hdrs) do gzhdrs
                    gzdopen(imgs) do gzimgs
                        out = load(gzhdrs, gzimgs; mmap=mmap,
                                   returntype=returntype)
                    end
                end
            else
                open(hdrs, imgs; mmap=mmap, returntype=returntype)
            end
        end
    end
end

function load(hdrs::Stream{format"HDR"}, imgs::Stream{format"IMG"};
                    mmap::Bool=false, returntype::Type=ImageMeta)
    version, needswap = checkfile(stream(hdrs))
    if returntype == ImageMeta
        header = ifelse(version == 1,
                        read(stream(hdrs.io), Nifti1Header, needswap),
                        read(stream(hdrs), Nifti2Header, needswap))
        extension = read_extension(stream(hdrs), header, needswap)
        volume = read_volume(stream(imgs), header, needswap; mmap=mmap)
        nii2img(volume, extension, header)
    elseif returntype == NiftiHeader
        ifelse(version == 1,
               read(stream(hdrs), Nifti1Header, needswap),
               read(stream(hdrs), Nifti2Header, needswap))
    else
        error("$returntype is not currently a supported NIfTI return type.")
    end
end

# save TODO
function save(f::File{format"NII"}, img::ImageMeta; version::Int)
    open(f, "w") do s
        hdr, vol, ext = img2nii(img)

        # ensure header, extension, and volume match up
        # Assume all is well and write
        write(stream(io), hdr)
        write_extension(stream(s), ext)
        write(stream(io), vol)
    end
end

# returns version, byteswap::Bool
function checkfile(s::IO)
    ret = read(s, Int32)
    seek(s,0)
    if ret == Int32(348)
       return 1, false
    elseif ret == Int32(540)
       return 2, false
    elseif ret == ntoh(Int32(348))
       return 1, true
    elseif ret == ntoh(Int32(540))
       return 2, true
    else
       return false, false
    end
end


