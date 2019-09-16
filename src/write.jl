function save(f::Formatted{DataFormat{:NII}},
              A::AbstractArray,
              version::Int=1;
              copyprops::Bool=false
             )

    anew = nieltransform(A)
    s = savestreaming(f, anew, version, copyprops=copyprops)
    write(s, extension(s))
    write(s, anew)
    flush(stream(s))
    return nothing
end

function savestreaming(s::File{DataFormat{:NII}}, A::AbstractArray, version::Int=1; copyprops::Bool=false)
    return savestreaming(open(s, "w"), A, version, copyprops=copyprops)
end

function savestreaming(s::Stream{DataFormat{:NII}}, A::AbstractArray, version::Int=1; copyprops::Bool=false)
    return savestreaming(s, ArrayInfo(A, copyprops), version)
end

function savestreaming(s::Stream{DataFormat{:NII},IOType}, info::ArrayInfo, version::Int=1) where IOType
    if file_extension(s) == ".gz"
        ret = ArrayStream(gzdopen(stream(s)), info)
    else
        ret = ArrayStream(stream(s), info)
    end
    if version == 1
        writehdr1(ret)
    else
        writehdr2(ret)
    end
    return ret
end


# TODO: this will need to support all the crazy CIFTI formats in the future
nieltransform(A::AbstractArray{<:BitTypes}) = A
function nieltransform(A::AbstractArray{T}) where T
    error("NIfTI.jl current doesn't support writing arrays of eltype $T.")
end

function writehdr1(s::ArrayStream{T}) where T
                                                              # offset - C name :: Type
                                                              # -----------------------

    sf = sform(s)
    a, qb, qc, qd, qx, qy, qz, xd, yd, zd, qfac = mat2quat(s)

    write(s, Int32(348))                                      # 0 - sizeof_hdr::Int32
    # data_type::NTuple{10,UInt8}, db_name::NTuple{18,UInt8},
    # extents::Int32, session_error::Int16, regular::Int8
    write(s, fill(Int8(0), 35))
    # FIXME change to todiminfo
    write(s, Int8(diminfo(s)))                                # 39 - dim_info::Int8
    write(s, convert(Vector{Int16}, nidim(s)))                # 40 - dim::NTuple{8,Int16}
    write(s, convert(Vector{Float32}, [intentparams(s)...]))  # 56 - intent_p1/2/3::Float32
    write(s, Int16(NiftiIntentsReverse[intent(s)]))           # 68 - intent_code::Int16
    write(s, Int16(NiftiDatatypesReverse[T]))                 # 70 - datatype::Int16
    write(s, Int16(sizeof(T)*8))                              # 72 - bitpix::Int16
    write(s, Int16(slicestart(s) - 1))                        # 74 - slice_start::Int16 - NOTE: go to zero based indexing
    write(s, convert(Vector{Float32}, pixdim(s, qfac)))       # 76 - pixdim::NTuple{8,Float32}
    dataoffset!(s, voxoffset(s, Val(1)))
    write(s, Float32(dataoffset(s)))                         # 108 - vox_offset::Float32
    write(s, Float32(scaleslope(s)))                          # 112 - scl_slope::Float32
    write(s, Float32(scaleintercept(s)))                      # 116 - scl_inter::Float32
    write(s, Int16(sliceend(s) - 1))                          # 120 - slice_end::Int16 - NOTE: go to zero based indexing
    write(s, Int8(NiftiSliceCodesReverse[slicecode(s)]))      # 122 - slice_code::Int8
    write(s, Int8(xyztunits(s)))                              # 123 - xyzt_units::Int8
    write(s, Float32(calmax(s)))                              # 124 - cal_max::Float32
    write(s, Float32(calmin(s)))                              # 128 - cal_min::Float32
    write(s, Float32(sliceduration(s)))                       # 132 - slice_duration::Float32
    write(s, Float32(toffset(s)))                             # 136 - toffset::Float32
    write(s, Int32[0, 0])                                     # 144 - glmax::Int32,glmin::Int32

    descrip = description(s)                                  # 148 - descrip::NTuple{80,UInt8}
    if length(descrip) > 80
        write(s, codeunits(descrip[1:80]))
    else
        write(s, descrip*String(fill(UInt8(0),80-length(descrip))))
    end

    aux = auxfiles(s)[1]
    if length(aux) > 24                                   # 228 - aux_file::NTuple{24,UInt8}
        write(s, codeunits(aux[1:80]))
    else
        write(s, aux*String(fill(UInt8(0),24-length(aux))))
    end

    write(s, Int16(NiftiXFormReverse[qformcode(s)]))      # 252 - qform_code::Int16
    write(s, Int16(NiftiXFormReverse[sformcode(s)]))      # 254 - sform_code::Int16

    # TODO: Test this
    write(s, Float32(qb))                               # 256 - quatern_b::Float32
    write(s, Float32(qc))                               # quatern_c::Float32
    write(s, Float32(qd))                               # quatern_d::Float32
    write(s, Float32(qx))                               # qoffset_x::Float32
    write(s, Float32(qy))                               # qoffset_y::Float32
    write(s, Float32(qz))                               # qoffset_z::Float32
    write(s, convert(Array{Float32}, sf[1,1:4]))        # 280 - srows
    write(s, convert(Array{Float32}, sf[2,1:4]))        # 296 - srows
    write(s, convert(Array{Float32}, sf[3,1:4]))        # 312 - srows
    write(s, codeunits(intentname(s)))                    # 328 - intent_name::NTuple{16,UInt8}
    write(s, [NP1_MAGIC...])                              # 344 - magic::NTuple{4,UInt8}

    return nothing
end

function writehdr2(s::ArrayStream{T}) where T
    sf = sform(s)
    a, qb, qc, qd, qx, qy, qz, xd, yd, zd, qfac = getquatern(s)

    write(s, Int32(540))                                 # 0 - sizeof_hdr::Int32
    write(s, [NP2_MAGIC...])                                  # 4 - magic::NTuple{8,UInt8}
    write(s, Int16(NiftiDatatypesReverse[T]))            # 12 - datatype::Int16
    write(s, Int16(sizeof(T)*8))                         # 14 - bitpix::Int16
    write(s, convert(Vector{Int64}, nidim(s)))           # 16 - dim::NTuple{8,Int64}
    write(s, convert(Vector{Float64}, [intentparams(s)...]))  # 80 - intent_[1:3]
    write(s, convert(Vector{Float64}, pixdim(s, qfac)))        # 104 - pixdim::NTuple{8,Float64}
    data_offset!(s, voxoffset(s, Val(2)))
    write(s, Int64(data_offset(s)))                      # 168 - vox_offset::Int64
    write(s, Float64(scaleslope(s)))                     # 176 - scl_slope::Float64
    write(s, Float64(scaleintercept(s)))                 # 184 - scl_inter::Float64
    write(s, Float64(calmax(s)))                         # 192 - cal_max::Float64
    write(s, Float64(calmin(s)))                         # 200 - cal_min::Float64
    write(s, Float64(sliceduration(s)))                  # 208 - slice_duration::Float64
    write(s, Float64(toffset(s)))                        # 216 - toffset::Float64
    write(s, Int64(slicestart(s) - 1))                   # 224 - slice_start::Int64 - NOTE: go to zero based indexing
    write(s, Int64(sliceend(s) - 1))                     # 232 - slice_end::Int64
    descrip = description(s)                             # 240 - descrip::NTuple{80,UInt8}
    if length(descrip) > 80
        write(s, codeunits(descrip[1:80]))
    else
        write(s, descrip*String(fill(UInt8(0),80-length(descrip))))
    end

    aux = auxfiles(s)[1]
    if length(aux) > 24                                  # 320 - aux_file::NTuple{24,UInt8}
        write(s, codeunits(aux[1:80]))
    else
        write(s, aux*String(fill(UInt8(0),24-length(aux))))
    end

    write(s, Int32(NiftiXFormReverse[qformcode(s)]))     # 344 - qform_code::Int32
    write(s, Int32(NiftiXFormReverse[sformcode(s)]))     # 348 - sform_code::Int32
    # TODO: should probably actually convert back
    write(s, Float64(qb))                              # 352 - quatern_b::Float64
    write(s, Float64(qc))                              # 360 - quatern_c::Float64
    write(s, Float64(qd))                              # 368 - quatern_d::Float64
    write(s, Float64(qx))                              # 376 - qoffset_x::Float64
    write(s, Float64(qy))                              # 384 - qoffset_y::Float64
    write(s, Float64(qz))                              # 392 - qoffset_z::Float64

    write(s, convert(Array{Float64}, sf[1,1:4]))      # 400 - srow_x
    write(s, convert(Array{Float64}, sf[2,1:4]))      # 432 - srow_y
    write(s, convert(Array{Float64}, sf[3,1:4]))      # 464 - srow_z

    write(s, Int32(NiftiSliceCodesReverse[slicecode(s)]))  # 496 - slice_code::Int32
    write(s, UInt32(xyztunits(s)))                         # 500 - xyzt_units::UInt32
    write(s, Int32(NiftiIntentsReverse[intent(s)]))        # 504 - intent_code::Int32
    write(s, codeunits(intentname(s)))                     # 508 - intent_name::NTuple{16,UInt8}
    write(s, Int8(diminfo(s)))                             # 524 - dim_info::Int8
    write(s, fill(UInt8(0), 15))                           # 525 - unused_str::NTuple{15,UInt8}
    return nothing
end
