
get_sform(x::NIVolume) = get_sform(x.header)
function get_sform(hdr::NIfTI1Header)
    if hdr.sform_code > 0
        return @inbounds Float32[
            hdr.srow_x[1]  hdr.srow_x[2]  hdr.srow_x[3]  hdr.srow_x[4]
            hdr.srow_y[1]  hdr.srow_y[2]  hdr.srow_y[3]  hdr.srow_y[4]
            hdr.srow_z[1]  hdr.srow_z[2]  hdr.srow_z[3]  hdr.srow_z[4]
            0              0              0              1
        ]
    else
        # TODO there's no sform in this case, should we just throw an error?
        return nothing
    end
end

get_qform(x::NIVolume) = get_qform(x.header)
function get_qform(hdr::NIfTI1Header)
    if hdr.qform_code <= 0
        return @inbounds Float32[
            hdr.pixdim[2]   0              0              0
            0               hdr.pixdim[3]  0              0
            0               0              hdr.pixdim[4]  0
            0               0              0              1
        ]
    else
        dx = hdr.pixdim[2]
        dy = hdr.pixdim[3]
        # aka qfac left handedness
        if hdr.pixdim[1] < 0  
            dz = -hdr.pixdim[4]
        else
            dz = hdr.pixdim[4]
        end
        b = hdr.quatern_b
        c = hdr.quatern_c
        d = hdr.quatern_d
        b2 = b*b
        c2 = c*c
        d2 = d*d
        a = 1 - b2 - c2 - d2
        if a < 1.e-7
            a = 1 / sqrt(b2 + c2 + d2)
            b *= a
            c *= a
            d *= a       # normalize (b,c,d) vector
            a = zero(a)  # a = 0 ==> 180 degree rotation
        else
            a = sqrt(a)   # angle = 2*arccos(a)
        end
        return Float32[
            (a*a+b*b-c*c-d*d)*dx   (2*b*c-2*a*d)*dy       (2*b*d+2*a*c)*dz      hdr.qoffset_x
            (2*b*c+2*a*d)*dx       (a*a+c*c-b*b-d*d)*dy   (2*c*d-2*a*b)*dz      hdr.qoffset_y
            (2*b*d-2*a*c)*dx       (2*c*d+2*a*b)*dy       (a*a+d*d-c*c-b*b)*dz  hdr.qoffset_z
            0                      0                      0                     1
        ]
    end
end

getaffine(x::NIVolume) = getaffine(x.header)

# Convert a NIfTI header to a 4x4 affine transformation matrix
function getaffine(hdr::NIfTI1Header)
    if hdr.sform_code > 0
        return get_sform(hdr)
    else
        return get_qform(hdr)
    end
end

# Set affine matrix of NIfTI header
function setaffine(h::NIfTI1Header, affine::Array{T,2}) where {T}
    size(affine, 1) == size(affine, 2) == 4 ||
        error("affine matrix must be 4x4")
    affine[4, 1] == affine[4, 2] == affine[4, 3] == 0 && affine[4, 4] == 1 ||
        error("last row of affine matrix must be [0 0 0 1]")
    h.qform_code = one(Int16)
    h.sform_code = one(Int16)
    h.pixdim = (zero(Float32), h.pixdim[2:end]...,)
    h.quatern_b = zero(Float32)
    h.quatern_c = zero(Float32)
    h.quatern_d = zero(Float32)
    h.qoffset_x = zero(Float32)
    h.qoffset_y = zero(Float32)
    h.qoffset_z = zero(Float32)
    h.srow_x = convert(Tuple{Vararg{Float32}}, (affine[1, :]...,))
    h.srow_y = convert(Tuple{Vararg{Float32}}, (affine[2, :]...,))
    h.srow_z = convert(Tuple{Vararg{Float32}}, (affine[3, :]...,))
    h
end

