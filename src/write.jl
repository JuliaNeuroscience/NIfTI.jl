# TODO:
# Should probably update the qform stuff too
# (although a lot of software seems to just skip this)


function niwrite(f::String, A::AbstractArray; version::Int=1, copyprops::Bool=false)
    open(f, "w") do io
        atransformed = nieltransform(A)
        s = ImageStream(io, atransformed; copyprops=copyprops)

        if version == 1
            writehdr1(s)
        else
            writehdr2(s)
        end
        write(s, extension(s))
        write(s, atransformed)  # ends on 1180064
    end
    return nothing
end

function niprops(A::AbstractArray)
    p = ImageProperties{:NII}()
    p["spacedirections"] = spacedirections(A)
end

nieltransform(A::AbstractArray{<:BitTypes}) = A
nieltransform(A::AbstractArray{T}) where T =
    error("NIfTI.jl current doesn't support writing arrays of eltype $T.")

function writehdr1(s::ImageStream{T}) where T
    write(s, Int32(348))                                      # 0 - sizeof_hdr::Int32
    # data_type::NTuple{10,UInt8}, db_name::NTuple{18,UInt8},
    # extents::Int32, session_error::Int16, regular::Int8
    write(s, fill(Int8(0), 35))
    write(s, Int8(diminfo(s)))                                # 39 - dim_info::Int8
    write(s, convert(Vector{Int16}, nidim(s)))                # 40 - dim::NTuple{8,Int16}
    write(s, convert(Vector{Float32}, [intentparams(s)...]))  # 56 - intent_p1/2/3::Float32
    write(s, Int16(NiftiIntentsReverse[intent(s)]))           # 68 - intent_code::Int16
    write(s, Int16(NiftiDatatypesReverse[T]))                 # 70 - datatype::Int16
    write(s, Int16(sizeof(T)*8))                              # 72 - bitpix::Int16
    write(s, Int16(slicestart(s)))                            # 74 - slice_start::Int16
    write(s, convert(Vector{Float32}, pixdim(s)))             # 76 - pixdim::NTuple{8,Float32}
    write(s, Float32(voxoffset(s, Val(1))))                   # 108 - vox_offset::Float32
    write(s, Float32(scaleslope(s)))                          # 112 - scl_slope::Float32
    write(s, Float32(scaleintercept(s)))                      # 116 - scl_inter::Float32
    write(s, Int16(sliceend(s)))                              # 120 - slice_end::Int16
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

    aux = auxfile(s)
    if length(aux) > 24                                   # 228 - aux_file::NTuple{24,UInt8}
        write(s, codeunits(aux[1:80]))
    else
        write(s, aux*String(fill(UInt8(0),24-length(aux))))
    end

    write(s, Int16(NiftiXFormReverse[qformcode(s)]))      # 252 - qform_code::Int16
    write(s, Int16(NiftiXFormReverse[sformcode(s)]))      # 254 - sform_code::Int16

    # TODO: should probably actually convert back
    write(s, zero(Float32))                               # 256 - quatern_b::Float32
    write(s, zero(Float32))                               # quatern_c::Float32
    write(s, zero(Float32))                               # quatern_d::Float32
    write(s, zero(Float32))                               # qoffset_x::Float32
    write(s, zero(Float32))                               # qoffset_y::Float32
    write(s, zero(Float32))                               # qoffset_z::Float32
    sf = sform(s)
    write(s, convert(Array{Float32}, [sf[1,:]...]))         # 280 - srows
    write(s, convert(Array{Float32}, [sf[2,:]...]))         # 296 - srows
    write(s, convert(Array{Float32}, [sf[3,:]...]))         # 312 - srows
    write(s, codeunits(intentname(s)))                    # 328 - intent_name::NTuple{16,UInt8}
    write(s, [NP1_MAGIC...])                              # 344 - magic::NTuple{4,UInt8}
    return nothing
end

function writehdr2(s::ImageStream{T}) where T
    write(s, Int32(540))                             # sizeof_hdr::Int32
    write(s, NP2_MAGIC)                              # magic::NTuple{8,UInt8}
    write(s, Int16(NiftiDatatypesReverse[T]))        # datatype::Int16
    write(s, Int16(sizeof(T)*8))                     # bitpix::Int16
    write(s, convert(Vector{Int64}, nidim(s)))           # dim::NTuple{8,Int64}
    write(s, convert(Vector{Float64}, intentparams(s)))
    write(s, convert(Vector{Float64}, pixdim(s)))        # pixdim::NTuple{8,Float64}
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
    if length(aux) > 24                                      # aux_file::NTuple{24,UInt8}
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

    sf = sform(s)
    write(s, convert(Array{Float64}, [sf[1,:]...]))            # srows
    write(s, convert(Array{Float64}, [sf[2,:]...]))            # srows
    write(s, convert(Array{Float64}, [sf[3,:]...]))            # srows

    write(s, Int32(NiftiSliceCodesReverse[slicecode(s)]))    # slice_code::Int32
    write(s, UInt32(xyztunits(s)))                            # xyzt_units::UInt32
    write(s, Int32(NiftiIntentsReverse[intent(s)]))          # intent_code::Int32
    write(s, codeunits(intentname(s)))                       # intent_name::NTuple{16,UInt8}
    write(s, Int8(diminfo(s)))                               # dim_info::Int8
    write(s, fill(UInt8(0), 15))                             # unused_str::NTuple{15,UInt8}
    return nothing
end
