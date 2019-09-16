function readhdr1(s::IO, p::AbstractDict)
    # Uncecessary fields
    skip(s, 35)

#    p["header"] = ImageProperties{:header}()
    d = read(s, Int8)
    freqdim!(p, Int(d & Int8(3) + 1))
    phasedim!(p, Int((d >> 2) & Int8(3) + 1))
    slicedim!(p, Int((d >> 4) + 1))
    N = Int(read(s, Int16))
    sz = ([Int(read(s, Int16)) for i in 1:N]...,)
    skip(s, (7-N)*2)  # skip filler dims

    # intent parameters
    intentparams!(p, Tuple(float.(read!(s, Vector{Int32}(undef, 3)))))
    intent!(p, num2intent(read(s, Int16)))
    T = get(NiftiDatatypes, read(s, Int16), UInt8)

    # skip bitpix
    skip(s, 2)
    slicestart!(p, Int(read(s, Int16)) + 1)  # to 1 based indexing

    qfac = Float64(read(s, Float32))
    pixdim = (Float64.(read!(s, Vector{Float32}(undef, N)))...,)

    skip(s, (7-N)*4)  # skip filler dims

    dataoffset!(p, Int(read(s, Float32)))
    scaleslope!(p, Float64(read(s, Float32)))
    scaleintercept!(p, Float64(read(s, Float32)))
    sliceend!(p, Int(read(s, Int16) + 1))  # to 1 based indexing
    slicecode!(p, numeric2slicecode(read(s, Int8)))

    xyzt_units = Int32(read(s, Int8))
    sp_units = get(NiftiUnits, xyzt_units & 0x07, 1)


    calmax!(p, read(s, Float32))
    calmin!(p, read(s, Float32))
    sliceduration!(p, read(s, Float32))
    toffset = read(s, Float32)

    skip(s, 8)
    description!(p, String(read(s, 80)))
    auxfiles!(p, [String(read(s, 24))])

    qformcode!(p, xform(read(s, Int16)))
    sformcode!(p, xform(read(s, Int16)))
    if qformcode(p) == UnkownSpace
        skip(s, 12)  # skip quaternion b/c space is unkown
        qx = read(s, Float64)
        qy = read(s, Float64)
        qz = read(s, Float64)
        qform!(p, quat2mat(zero(Float64), zero(Float64), zero(Float64),
                           zero(Float64), zero(Float64), zero(Float64),
                           dx, dy, dz, zero(Float64)))
    else
        quaternb!(p, Float64(read(s, Float32)))
        quaternc!(p, Float64(read(s, Float32)))
        quaternd!(p, Float64(read(s, Float32)))
        qx = Float64(read(s, Float32))
        qy = Float64(read(s, Float32))
        qz = Float64(read(s, Float32))

        qform!(p, quat2mat(quaternb(p), quaternc(p), quaternd(p),
                           qx, qy, qz, dx, dy, dz, qfac))

        if sformcode(p) == UnkownSpace
            skip(s, 48)
            p["spacedirections"] = (Tuple(qform(p)[1,1:3]), Tuple(qform(p)[2,1:3]), Tuple(qform(p)[3,1:3]))
            dimnames = orientation(qform(p))
        else
            sform!(p, MMatrix{4,4,Float64}(vcat(Float64.(read!(s, Matrix{Float32}(undef, (1,4)))),
                                                Float64.(read!(s, Matrix{Float32}(undef, (1,4)))),
                                                Float64.(read!(s, Matrix{Float32}(undef, (1,4)))),
                                                Float64[0, 0, 0, 1]')))
            p["spacedirections"] = (Tuple(sform(p)[1,1:3]), Tuple(sform(p)[2,1:3]), Tuple(sform(p)[3,1:3]))
        end
    end
    dimnames = affinematrix(p)
    axs = (range(qx, step=pixdim[1], length=sz[1])*sp_units,)
    if N > 1
        axs = (axs..., range(qy, step=pixdim[2], length=sz[2])*sp_units)
    end
    if N > 2
        axs = (axs..., range(qz, step=pixdim[3], length=sz[3])*sp_units)
    end
    if N > 3
        dimnames = (dimnames..., :time)
        tu = get(NiftiUnits, xyzt_units & 0x38, 1)
        axs = (axs..., range(toffset, step=pixdim[4], length=sz[4])*tu)
    end
    if N > 4
        dimnames = (dimnames..., intentaxis(intent)...)
        axs = (axs..., map(i->range(one(Float64), step=pixdim[i], lenght=sz[i]), 5:N)...)
    end
    intentname!(p, String(read(s, 16)))
    magicbytes!(p, read(s, 4))

    return ArrayInfo{T}(NamedTuple{dimnames}(axs), p)
end
