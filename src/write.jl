# TODO:
# niupdate!

function voxoffset(A::AbstractArray, v::Val{1})
    if isempty(extension(A))
        SIZEOF_HDR1+4
    else
        mapreduce(esize, +, extension(A))+SIZEOF_HDR1+4
    end
end

function voxoffset(A::AbstractArray, v::Val{2})
    if isempty(extension(A))
        SIZEOF_HDR2+4
    else
        mapreduce(esize, +, extension(A))+SIZEOF_HDR2+4
    end
end

# Gets dim to be used in header
function nidim(x::AbstractArray)
    dim = ones(Int16, 8)
    dim[1] = ndims(x)
    dim[2:dim[1]+1] = [size(x)...]
    (dim...,)
end

pixdim(A::AbstractArray{T,N}) where {T,N} =
    map(i->ustrip(step(i.val)), AxisArrays.axes(A))::NTuple{N,<:Any}
xyztunits(A::AbstractArray) =
    (NiftiUnitsReverse[spatunits(A)[1]] & 0x07) | (NiftiUnitsReverse[timeunits(A)] & 0x38)

function toffset(A::AbstractArray{T,N}) where {T,N}
    ax = timeaxis(A)
    if isa(ax, Nothing)
        return Float64(1)
    else
        Float64(ustrip(ax[1]))
    end
end

function writehdr1(s::ImageStream{T}) where T
    write(s, Int32(348))                                  # sizeof_hdr::Int32
    # data_type::NTuple{10,UInt8}, db_name::NTuple{18,UInt8},
    # extents::Int32, session_error::Int16, regular::Int8
    write(s, fill(Int8(0), 35))
    write(s, Int8(diminfo(s)))                            # dim_info::Int8
    write(s, reinterpret(Int16, nidim(s)))                # dim::NTuple{8,Int16}
    write(s, reinterpret(Float32, intentparams(s)))       # intent_p1/2/3::Float32
    write(s, Int16(NiftiIntentsReverse[intent(s)]))       # intent_code::Int16
    write(s, Int16(NiftiDatatypesReverse[T]))             # datatype::Int16
    write(s, Int16(sizeof(T)*8))                          # bitpix::Int16
    write(s, Int16(slicestart(s)))                        # slice_start::Int16
    write(s, reinterpret(Float32, pixdim(s)))             # pixdim::NTuple{8,Float32}
    write(s, Float32(voxoffset(s, Val(1))))               # vox_offset::Float32
    write(s, Float32(scaleslope(s)))                      # scl_slope::Float32
    write(s, Float32(scaleintercept(s)))                  # scl_inter::Float32
    write(s, Int16(sliceend(s)))                          # slice_end::Int16
    write(s, Int8(NiftiSliceCodesReverse[slicecode(s)]))  # slice_code::Int8
    write(s, Int8(xyztunits(s)))                          # xyzt_units::Int8
    write(s, Float32(calmax(s)))                          # cal_max::Float32
    write(s, Float32(calmax(s)))                          # cal_min::Float32
    write(s, Float32(sliceduration(s)))                   # slice_duration::Float32
    write(s, Float32(toffset(s)))                         # toffset::Float32
    write(s, Int32[0, 0])                                 # glmax::Int32,glmin::Int32

    descrip = description(s)                              # descrip::NTuple{80,UInt8}
    if length(descrip) > 80
        write(s, codeunits(descrip[1:80]))
    else
        write(s, descrip*String(fill(UInt8(0),80-length(descrip))))
    end

    aux = auxfile(s)  
    if length(aux) > 24                                   # aux_file::NTuple{24,UInt8}
        write(s, codeunits(aux[1:80]))
    else
        write(s, aux*String(fill(UInt8(0),24-length(aux))))
    end

    write(s, Int16(NiftiXFormReverse[qformcode(s)]))      # qform_code::Int16
    write(s, Int16(NiftiXFormReverse[sformcode(s)]))      # sform_code::Int16

    # TODO: should probably actually convert back
    write(s, zero(Float32))                               # quatern_b::Float32
    write(s, zero(Float32))                               # quatern_c::Float32
    write(s, zero(Float32))                               # quatern_d::Float32
    write(s, zero(Float32))                               # qoffset_x::Float32
    write(s, zero(Float32))                               # qoffset_y::Float32
    write(s, zero(Float32))                               # qoffset_z::Float32
    write(s, convert(Array{Float32}, sform(s)')[:])       # srows... (switch to row major and unwrap)
    write(s, codeunits(intentname(s)))                    # intent_name::NTuple{16,UInt8}
    write(s, [NP1_MAGIC...])                              # magic::NTuple{8,UInt8}
end

function writehdr2(s::ImageStream{T}) where T
    write(s, Int32(540))                             # sizeof_hdr::Int32
    write(s, NP2_MAGIC)                              # magic::NTuple{8,UInt8}
    write(s, Int16(NiftiDatatypesReverse[T]))        # datatype::Int16
    write(s, Int16(sizeof(T)*8))                     # bitpix::Int16
    write(s, reinterpret(Int64, nidim(s)))           # dim::NTuple{8,Int64}
    write(s, reinterpret(Float64, intentparams(s)))
    write(s, reinterpret(Float64, pixdim(s)))        # pixdim::NTuple{8,Float64}
    write(s, Int64(voxoffset(s, Val(2))))            # vox_offset::Int64
    write(s, Float64(scaleslope(s)))                 # scl_slope::Float64
    write(s, Float64(scaleintercept(s)))             # scl_inter::Float64
    write(s, Float64(calmax(s)))                     # cal_max::Float64
    write(s, Float64(calmin(s)))                     # cal_min::Float64
    write(s, Float64(sliceduration(s)))              # slice_duration::Float64
    write(s, Float64(toffset(s)))                    # toffset::Float64
    write(s, Int64(slicestart(s)))                   # slice_start::Int64
    write(s, Int64(sliceend(s)))                     # slice_end::Int64
                            
    descrip = description(s)                         # descrip::NTuple{80,UInt8}
    if length(descrip) > 80
        write(s, codeunits(descrip[1:80]))
    else
        write(s, descrip*String(fill(UInt8(0),80-length(descrip))))
    end

    aux = auxfile(s)  
    if length(aux) > 24                                     # aux_file::NTuple{24,UInt8}
        write(s, codeunits(aux[1:80]))
    else
        write(s, aux*String(fill(UInt8(0),24-length(aux))))
    end

    write(s, Int32(NiftiXFormReverse[qformcode(s)]))         # qform_code::Int32
    write(s, Int32(NiftiXFormReverse[sformcode(s)]))         # sform_code::Int32
    # TODO: should probably actually convert back
    write(s, zero(Float64))                                  # quatern_b::Float64
    write(s, zero(Float64))                                  # quatern_c::Float64
    write(s, zero(Float64))                                  # quatern_d::Float64
    write(s, zero(Float64))                                  # qoffset_x::Float64
    write(s, zero(Float64))                                  # qoffset_y::Float64
    write(s, zero(Float64))                                  # qoffset_z::Float64
    write(s, convert(Array{Float64}, sform(s)')[:])          # srows... (switch to row major and unwrap)
    write(s, Int32(NiftiSliceCodesReverse[slicecode(s)]))    # slice_code::Int32
    write(s, Int32(xyztunits(s)))                            # xyzt_units::UInt32
    write(s, Int32(NiftiIntentsReverse[intent(s)]))          # intent_code::Int32
    write(s, codeunits(intentname(s)))                       # intent_name::NTuple{16,UInt8}
    write(s, Int8(diminfo(s)))                               # dim_info::Int8
    write(s, fill(UInt8(0), 15))                             # unused_str::NTuple{15,UInt8}
end
