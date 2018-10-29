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
function read_volume(io::IO, vox_offset::F, datatype::Int16,
                     scl_slope::F, scl_inter::F,
                     needswap::Bool=false; mmap::Bool=false) where {F<:AbstractFloat}
    seek(io, vox_offset)
    dtype = NIfTI_DT_BITSTYPES[datatype]
    if mmap
        vol = Mmap.mmap(io, Array{dtype,length(dims)}, dims)
    else
        vol = read!(io, Array{dtype}(undef, dims))
    end
    if needswap && sizeof(eltype(vol)) > 1
        vol = mappedarray((ntoh, hton), vol)
    end
    # scale if slope isn't zero
    if scl_slope != F(0) && dtype != RGB24  # TODO finalize RGB types
        vol = vol.*f.header.scl_slope.+f.header.scl_inter
    end
    vol
end

# TODO
function write_volume(io::IO, vol::AbstractArray)
end

function load(f::File{format"NII"}; mode="r", mmap::Bool=false,
              returntype::String="image")
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
    version, needswap = checkfile(s)
    if returntype == ImageMeta
        header = ifelse(version == 1,
                        read(s.io, Nifti1Header, needswap),
                        read(s.io, Nifti2Header, needswap))
        intent = makeintent(header.intent_code, header.intent_p1,
                            header.intent_p2, header.intent_p3,
                            header.intent_name, header.dim)
        extension = read_extension(s, header, needswap)
        volume = read_volume(s, header.vox_offset, header.datatype,
                             needswap; mmap=mmap)
        nii2img(volume, extension, header, intent)
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
              returntype::String="image")
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
    version, needswap = checkfile(hdrs)
    if returntype == ImageMeta
        header = ifelse(version == 1,
                read(hdrs.io, Nifti1Header, needswap),
                read(hdrs.io, Nifti2Header, needswap))
        intent = getintent(header)
        extension = read_extension(hdrs, header, needswap)
        volume = read_volume(imgs, header.vox_offset, header.datatype,
                             needswap; mmap=mmap)
        nii2img(volume, extension, header, intent)
    elseif returntype == NiftiHeader
        ifelse(version == 1,
               read(hdrs.io, Nifti1Header, needswap),
               read(hdrs.io, Nifti2Header, needswap))
    else
        error("$returntype is not currently a supported NIfTI return type.")
    end
end

# save TODO
function save(f::File{format"NII"}, img::ImageMeta)
    open(f, "w") do io
        hdr = getheader(img)
        ext = getextension(img, hdr)
        vol = getimage(img, hdr, ext)

        # ensure header, extension, and volume match up
        if imgchecks(volume, extension, hdr)
        end

        # Assume all is well and write
        write_header(io, hdr)
        write_extension(io, extension)
        write_volume(io, volume)
    end
end

# returns version, byteswap::Bool
function checkfile(::Stream)
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

# FileIO stuff

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

function detecthdr(io::IO)
    detectnii(io)
end

function detectimg(io::IO)
    img_file = filename(io)
    hdr_file = gethdr(img_file)
    open(hdr_file, "r") do io
        detectnii(io)
    end
end

add_format(format"NII", detectnii, [".nii"], [:NIfTI])
# HDR is already FileIO registered
#add_format(format"HDR", detecthdr, [".hdr"], [:NIfTI])
add_format(format"IMG", detectimg, [".img"], [:NIfTI])
