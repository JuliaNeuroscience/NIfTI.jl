function readhdr2(s::IO, p::AbstractDict)
#    p["header"] = ImageProperties{:header}()
    magicbytes!(p, read(s, 8))
    T = get(NiftiDatatypes, read(s, Int16), UInt8)
    skip(s, 2)  # skip bitpix

    N = read(s, Int64)
    sz = Tuple(read!(s, Vector{Int64}(undef, N)))
    skip(s, (7-N)*8)  # skip filler dims

    intentparams!(p, (read!(s, Vector{Float64}(undef, 3))...,))

    qfac = read(s, Float64)
    pixdim = (read!(s, Vector{Float64}(undef, N))...,)
    skip(s, (7-N)*8)  # skip filler dims

    dataoffset!(p, read(s, Int64))
    scaleslope!(p, read(s, Float64))
    scaleintercept!(p, read(s, Float64))
    calmax!(p, read(s, Float64))
    calmin!(p, read(s, Float64))

    sliceduration!(p, read(s, Float64))
    toffset = read(s, Float64)

    slicestart!(p, read(s, Int64) + 1)  # to 1 based indexing
    sliceend!(p, read(s, Int64) + 1)
    description!(p, String(read(s, 80)))
    auxfiles!(p, [String(read(s, 24))])
    qformcode!(p, xform(read(s, Int32)))
    sformcode!(p, xform(read(s, Int32)))
    if qformcode(p) == UnkownSpace
        skip(s, 12)  # skip quaternion b/c space is unkown
        qx = read(s, Float64)
        qy = read(s, Float64)
        qz = read(s, Float64)
        qform!(p, quat2mat(zero(Float64), zero(Float64), zero(Float64),
                           zero(Float64), zero(Float64), zero(Float64),
                           dx, dy, dz, zero(Float64)))
    else
        quaternb!(p, read(s, Float64))
        quaternc!(p, read(s, Float64))
        quaternd!(p, read(s, Float64))
        qx = read(s, Float64)
        qy = read(s, Float64)
        qz = read(s, Float64)

        qform!(p, quat2mat(quaternb(p), quaternc(p), quaternd(p),
                           qx, qy, qz, dx, dy, dz, qfac))

        if sformcode(p) == UnkownSpace
            skip(s, 48)
            p["spacedirections"] = (Tuple(qform(p)[1,1:3]), Tuple(qform(p)[2,1:3]), Tuple(qform(p)[3,1:3]))
        else
            sform!(p, MMatrix{4,4,Float64}(vcat(read!(s, Matrix{Float64}(undef, (1,4))),
                                                read!(s, Matrix{Float64}(undef, (1,4))),
                                                read!(s, Matrix{Float64}(undef, (1,4))),
                                                Float64[0, 0, 0, 1]')))
            p["spacedirections"] = (Tuple(sform(p)[1,1:3]), Tuple(sform(p)[2,1:3]), Tuple(sform(p)[3,1:3]))
        end
    end
    dimnames = orientation(affinematrix(p))
    slicecode!(p, numeric2slicecode(read(s, Int32)))

    xyzt_units = read(s, Int32)
    sp_units = get(NiftiUnits, xyzt_units & 0x07, 1)
    intent!(p, num2intent(read(s, Int32)))

    intentname!(p, String(read(s, 16)))
    d = read(s, Int8)
    freqdim!(p, Int(d & Int8(3) + 1))
    phasedim!(p, Int((d >> 2) & Int8(3) + 1))
    slicedim!(p, Int((d >> 4) + 1))
 
    skip(s, 15)

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
        dimnames = (dimnames..., intentdimnames(intent, N)...)
        axs = (axs..., map(i->range(one(Float64), step=pixdim[i], lenght=sz[i]), 5:N)...)
    end

    return ArrayInfo{T}(NamedTuple{dimnames}(axs), p)
end
