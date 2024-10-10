
@inline function _det(r11, r12, r13, r21, r22, r23, r31, r32, r33)
    r11 * r22 * r33 -
    r11 * r32 * r23 -
    r21 * r12 * r33 +
    r21 * r32 * r13 +
    r31 * r12 * r23 -
    r31 * r22 * r13
end

@inline function _mul_trace(
    x11, x12, x13, x21, x22, x23, x31, x32, x33,
    y11, y12, y13, y21, y22, y23, y31, y32, y33
)

    return (x11 * y11 + x12 * y21 + x13 * y31) +  # z11
           (x21 * y12 + x22 * y22 + x23 * y32) +  # z22
           (x31 * y13 + x32 * y23 + x33 * y33)    # z33
end

get_sform(x::NIVolume) = get_sform(x.header)
function get_sform(hdr::NIfTIHeader)
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
function get_qform(hdr::NIfTIHeader)
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
"""
    getaffine(x::NIVolume)

Gets a 4x4 affine transformation volume's header's sform
"""
getaffine(x::NIVolume) = getaffine(x.header)

# Convert a NIfTI header to a 4x4 affine transformation matrix
"""
    getaffine(hdr::NIfTIHeader)

Gets a 4x4 affine transformation matrix from a header's sform
"""
function getaffine(hdr::NIfTIHeader)
    if hdr.sform_code > 0
        return get_sform(hdr)
    else
        return get_qform(hdr)
    end
end

"""
    function setaffine(h::NIfTIHeader, affine::Array{T,2}) where {T}

Set the affine of a `NIfTIHeader` to 4x4 affine matrix `affine`
"""
function setaffine(h::NIfTIHeader, affine::Array{T,2}) where {T}
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

"""
    orientation(img)::Tuple{Symbol,Symbol,Symbol}

Returns a tuple providing the orientation of a NIfTI image.
"""
orientation(x) = orientation(x.header)
function orientation(hdr::NIfTIHeader)
    if hdr.sform_code > 0
        return @inbounds _dir2ori(
            hdr.srow_x[1], hdr.srow_x[2], hdr.srow_x[3],
            hdr.srow_y[1], hdr.srow_y[2], hdr.srow_y[3],
            hdr.srow_z[1], hdr.srow_z[2], hdr.srow_z[3]
        )
    elseif hdr.qform_code <= 0
        return @inbounds _dir2ori(
            hdr.pixdim[2], 0, 0,
            0, hdr.pixdim[3], 0,
            0, 0, hdr.pixdim[4]
        )
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
        return _dir2ori(
            (a*a+b*b-c*c-d*d)*dx,  (2*b*c-2*a*d)*dy,      (2*b*d+2*a*c)*dz,
            (2*b*c+2*a*d)*dx,      (a*a+c*c-b*b-d*d)*dy,  (2*c*d-2*a*b)*dz,
            (2*b*d-2*a*c)*dx,      (2*c*d+2*a*b)*dy,      (a*a+d*d-c*c-b*b)*dz,
        )
    end
end

_encoding_name(x) = _encoding_name(Int(x))
@inline function _encoding_name(x::Int)
    if x === 1
        return :left
    elseif x === -1
        return :right
    elseif x === 2
        return :posterior
    elseif x === -2
        return :anterior
    elseif x === 3
        return :inferior
    elseif x === -3
        return :superior
    else
        error("$x does not map to a dimension name.")
    end
end

function _dir2ori(xi, xj, xk, yi, yj, yk, zi, zj, zk)
    # Normalize column vectors to get unit vectors along each ijk-axis
    # normalize i axis
    val = sqrt(xi*xi + yi*yi + zi*zi)
    if val == 0
        error("Invalid rotation directions.")
    end
    xi /= val
    yi /= val
    zi /= val

    # normalize j axis
    val = sqrt(xj*xj + yj*yj + zj*zj)
    if val == 0
        error("Invalid rotation directions.")
    end
    xj /= val
    yj /= val
    zj /= val

    # orthogonalize j axis to i axis, if needed
    val = xi*xj + yi*yj + zi* zj  # dot product between i and j
    if abs(val) > .0001
        xj -= val*xi
        yj -= val*yi
        zj -= val*zi

        val = sqrt(xj*xj + yj*yj + zj*zj)  # must renormalize
        if val == 0
            error("The first and second dimensions cannot be parallel.")
        end
        xj /= val
        yj /= val
        zj /= val
    end

    # normalize k axis; if it is zero, make it the cross product i x j
    val = sqrt(xk*xk + yk*yk + zk*zk)
    if val == 0
        xk = yi*zj-zi*yj
        yk = zi*xj-zj*xi
        zk = xi*yj-yi*xj
    else
        xk = xk/val
        yk = yk/val
        zk = zk/val
    end

    # orthogonalize k to i
    val = xi*xk + yi*yk + zi*zk  # dot product between i and k
    if abs(val) > 0.0001
        xk -= val*xi
        yk -= val*yi
        zk -= val*zi

        # must renormalize
        val = sqrt(xk*xk + yk*yk + zk*zk)
        if val == 0
            return 0  # I think this is suppose to be an error output
        end
        xk /= val
        yk /= val
        zk /= val
    end

    # orthogonalize k to j */
    val = xj*xk + yj*yk + zj*zk  # dot product between j and k
    if abs(val) > 0.0001
        xk -= val*xj
        yk -= val*yj
        zk -= val*zj

        val = sqrt(xk*xk + yk*yk + zk*zk)
        if val == 0
            return 0  # bad
        end
        xk /= val
        yk /= val
        zk /= val
    end

    # at this point Q is the rotation matrix from the (i,j,k) to (x,y,z) axes
    detQ = _det(xi, xj, xk, yi, yj, yk, zi, zj, zk)
    # if( detQ == 0.0 ) return ; /* shouldn't happen unless user is a DUFIS */

    # Build and test all possible +1/-1 coordinate permutation matrices P;
    # then find the P such that the rotation matrix M=PQ is closest to the
    # identity, in the sense of M having the smallest total rotation angle.

    # Despite the formidable looking 6 nested loops, there are
    # only 3*3*3*2*2*2 = 216 passes, which will run very quickly.
    vbest = -666
    ibest = pbest=qbest=rbest= 1
    jbest = 2
    kbest = 3
    for (i, j, k) in ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))
        for p in (-1, 1)           # p,q,r are -1 or +1
            for q in (-1, 1)       # and go into rows 1,2,3
                for r in (-1, 1)
                    p11, p12, p13 = _nval_other_zero(i, p)
                    p21, p22, p23 = _nval_other_zero(j, q)
                    p31, p32, p33 = _nval_other_zero(k, r)
                    #=
                    P[1,i] = p
                    P[2,j] = q
                    P[3,k] = r
                    detP = det(P)  # sign of permutation
                    =#
                    detP = _det(p11, p12, p13, p21, p22, p23, p31, p32, p33)
                    # doesn't match sign of Q
                    if detP * detQ >= 0.0
                        # angle of M rotation = 2.0 * acos(0.5 * sqrt(1.0 + trace(M)))
                        # we want largest trace(M) == smallest angle == M nearest to I
                        val = _mul_trace(
                            p11, p12, p13, p21, p22, p23, p31, p32, p33,
                            xi, xj, xk, yi, yj, yk, zi, zj, zk
                        )
                        if val > vbest
                            vbest = val
                            ibest = i
                            jbest = j
                            kbest = k
                            pbest = p
                            qbest = q
                            rbest = r
                        end
                    end
                end
            end
        end
    end
    # At this point ibest is 1 or 2 or 3; pbest is -1 or +1; etc.

    # The matrix P that corresponds is the best permutation approximation
    # to Q-inverse; that is, P (approximately) takes (x,y,z) coordinates
    # to the (i,j,k) axes.

    # For example, the first row of P (which contains pbest in column ibest)
    # determines the way the i axis points relative to the anatomical
    # (x,y,z) axes.  If ibest is 2, then the i axis is along the y axis,
    # which is direction P2A (if pbest > 0) or A2P (if pbest < 0).

    # So, using ibest and pbest, we can assign the output code for
    # the i axis.  Mutatis mutandis for the j and k axes, of course.

    return (_encoding_name(ibest*pbest), _encoding_name(jbest*qbest), _encoding_name(kbest*rbest))
end

@inline function _nval_other_zero(n, val)
    if n === 1
        return val, 0, 0
    elseif n === 2
        return 0, val, 0
    else
        return 0, 0, val
    end
end


