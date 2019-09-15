function readhdr2(s::IO, p::AbstractDict)
#    p["header"] = ImageProperties{:header}()
    magicbytes!(p, (read(s, 8)...,))
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
        qform!(p, MMatrix{4,4,Float64,16}([pixdim[2] 0         0          0 0
                                                     0 pixdim[3]          0 0
                                                     0         0  pixdim[4] 0
                                                     0         0          0 1]))
        qx = Float64(read(s, Float64))
        qy = Float64(read(s, Float64))
        qz = Float64(read(s, Float64))
    else
        b = Float64(read(s, Float64))
        c = Float64(read(s, Float64))
        d = Float64(read(s, Float64))
        qx = Float64(read(s, Float64))
        qy = Float64(read(s, Float64))
        qz = Float64(read(s, Float64))
        a = 1 - (b*b + c*c + d*d)
        if a < 1.e-7                   # special case
            a = 1 / sqrt(b*b+c*c+d*d)
            b *= a
            c *= a
            d *= a                   # normalize (b,c,d) vector
            a = zero(Float64)        # a = 0 ==> 180 degree rotation
        else
            a = sqrt(a)              # angle = 2*arccos(a)
        end
        # make sure are positive
        xd = pixdim[1] > 0 ? pixdim[1] : one(Float64)
        yd = pixdim[2] > 0 ? pixdim[2] : one(Float64)
        zd = pixdim[3] > 0 ? pixdim[3] : one(Float64)
        zd = qfac < 0 ? -zd : zd
        qform!(p, MMatrix{4,4,Float64}([[((a*a+b*b-c*c-d*d)*xd),       (2*(b*c-a*d)*yd),       (2*(b*d+a*c)*zd),   qx]'
                                        [     (2*(b*c+a*d )*xd), ((a*a+c*c-b*b-d*d)*yd),       (2*(c*d-a*b)*zd),   qy]'
                                        [      (2*(b*d-a*c)*xd),       (2*(c*d+a*b)*yd), ((a*a+d*d-c*c-b*b)*zd),   qz]'
                                        [                     0,                      0,                      0, qfac]']))
        if sformcode(p) == UnkownSpace
            skip(s, 48)
            p["spacedirections"] = (Tuple(qform(p)[1,1:3]), Tuple(qform(p)[1,1:3]), Tuple(qform(p)[3,1:3]))
            dimnames = orientation(qform(p))
        else
            sform!(p, MMatrix{4,4,Float64}(vcat(read!(s, Matrix{Float64}(undef, (1,4))),
                                                read!(s, Matrix{Float64}(undef, (1,4))),
                                                read!(s, Matrix{Float64}(undef, (1,4))),
                                                Float64[0, 0, 0, 1]')))
            p["spacedirections"] = (Tuple(sform(p)[1,1:3]), Tuple(sform(p)[1,1:3]), Tuple(sform(p)[3,1:3]))
            dimnames = orientation(sform(p))
        end
    end
 


    slicecode!(p, numeric2slicecode(read(s, Int32)))

    xyzt_units = read(s, Int32)
    sp_units = get(NiftiUnits, xyzt_units & 0x07, 1)
    p["header"]["intent"] = get(NiftiIntents, read(s, Int32), NoIntent)

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

    extension!(p, read(s, p, NiftiExtension))

    return ArrayInfo{T}(NamedTuple{dimnames}(axs), p)
end
