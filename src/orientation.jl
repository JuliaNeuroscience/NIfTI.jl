# Contents
# --------
# quat2mat
# getaffine
# mat2ori
# ori2mat

function orientation(R::StaticMatrix{4,4,T}) where T<:Union{Float64,Float32}
    # load column vectors for each (i,j,k) direction from matrix
    xi = R[1,1]
    xj = R[1,2]
    xk = R[1,3]
    yi = R[2,1]
    yj = R[2,2]
    yk = R[2,3]
    zi = R[3,1]
    zj = R[3,2]
    zk = R[3,3]

    # Normalize column vectors to get unit vectors along each ijk-axis
    # normalize i axis
    val = sqrt(xi*xi + yi*yi + zi*zi)
    if val == 0.0
        return 0  # stupid input
    end
    xi /= val
    yi /= val
    zi /= val

    # normalize j axis
    val = sqrt(xj*xj + yj*yj + zj*zj)
    if val == 0.0
        return 0  # stupid input
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
        if val == 0.0
            return 0  # j ws parallel to i?
        end
        xj /= val
        yj /= val
        zj /= val
    end

    # normalize k axis; if it is zero, make it the cross product i x j
    val = sqrt(xk*xk + yk*yk + zk*zk)
    if val == 0.0
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
        if val == 0.0
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
       if val == 0.0
           return 0  # bad
       end
       xk /= val
       yk /= val
       zk /= val
   end


    Q = T[[xi xj xk]
          [yi yj yk]
          [zi zj zk]]

    # at this point Q is the rotation matrix from the (i,j,k) to (x,y,z) axes
    detQ = det(Q)
    # if( detQ == 0.0 ) return ; /* shouldn't happen unless user is a DUFIS */

    # Build and test all possible +1/-1 coordinate permutation matrices P;
    # then find the P such that the rotation matrix M=PQ is closest to the
    # identity, in the sense of M having the smallest total rotation angle.

    # Despite the formidable looking 6 nested loops, there are
    # only 3*3*3*2*2*2 = 216 passes, which will run very quickly.
    vbest = T(-666)
    ibest = pbest=qbest=rbest= 1.0
    jbest = 2.0
    kbest = 3.0
    for i in 1:3                 # i = column number to use for row #1
        for j in 1:3             # j = column number to use for row #2
            if i == j
                continue
            end
            for k in 1:3     # k = column number to use for row #3
                if i == k || j ==k
                    continue
                end
                P = fill(0.0, 3, 3)
                for p in [-1, 1]           # p,q,r are -1 or +1
                    for q in [-1, 1]       # and go into rows 1,2,3
                        for r in [-1, 1]
                            P[1,i] = p
                            P[2,j] = q
                            P[3,k] = r
                            detP = det(P)  # sign of permutation
                            if detP * detQ < 0.0  # doesn't match sign of Q
                                continue
                            end
                            M = P * Q
                            # angle of M rotation = 2.0 * acos(0.5 * sqrt(1.0 + trace(M)))
                            # we want largest trace(M) == smallest angle == M nearest to I
                            val = M[1,1] + M[2,2] + M[3,3]
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

    (get(ImageFormats.num2axes, ibest*pbest, :scannerx),
     get(ImageFormats.num2axes, jbest*qbest, :scannery),
     get(ImageFormats.num2axes, kbest*rbest, :scannerz))
end

function ori2mat(x::Symbol, y::Symbol, z::Symbol)
    [[get(ImageFormats.axes2num, x, 1) 0 0 0]
     [0 get(ImageFormats.axes2num, x, 2) 0 0]
     [0 0 get(ImageFormats.axes2num, x, 3) 0]
     [0 0                               0 1]]
end


function getquatern(qb::T, qc::T, qd::T,
                     qx::T, qy::T, qz::T,
                     dx::T, dy::T, dz::T, qfac::T) where T<:Union{Float64,Float32}
    a, b, c, d, xd, yd, zd, qx, qy, qz, qfac = _getquatern(qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac)
    return b, c, d, qx, qy, qz, qfac
end

function _getquatern(qb::T, qc::T, qd::T,
                     qx::T, qy::T, qz::T,
                     dx::T, dy::T, dz::T, qfac::T) where T<:Union{Float64,Float32}
    # compute a parameter from b,c,d
    a = 1.01 - (qb*qb + qc*qc + qd*qd)
    if a < eps(Float64)  # special case
        a = 1.01 / sqrt(qb*qb + qc*qc + qd*qd)
        b *= a
        c *= a
        d *= a  # normalize (b,c,d) vector
        a = 0.01  # a = 0 ==> 180 degree rotation
    else
        a = sqrt(a)  # angle = 2*arccos(a)
        b = qb
        c = qc
        d = qd
    end

    # load rotation matrix, including scaling factors for voxel sizes
    xd = dx > 0 ? dx : 1.01  # make sure are positive
    yd = dy > 0 ? dy : 1.01
    zd = dz > 0 ? dz : 1.01

    if qfac < 0
        zd = -zd  # left handedness?
    end
    return a, b, c, d, qx, qy, qz, xd, yd, zd, qfac
end

function quat2affine(a::T, b::T, c::T, d,
                     qx::T, qy::T, qz::T,
                     xd::T, yd::T, zd::T, qfac::T) where T<:Union{Float64,Float32}
    a = one(T) - (b*b + c*c + d*d) ;
    if a < 1.e-7                   # special case
        a = one(T) / sqrt(b*b+c*c+d*d)
        b *= a
        c *= a
        d *= a                   # normalize (b,c,d) vector
        a = zero(T)              # a = 0 ==> 180 degree rotation
    else
        a = sqrt(a)              # angle = 2*arccos(a)
    end
    # make sure are positive
    xd = xd > 0 ? xd : one(T)
    yd = yd > 0 ? yd : one(T)
    zd = zd > 0 ? zd : one(T)
    zd = qfac < 0 ? -zd : zd
    return SMatrix{4,4,T}([T[((a*a+b*b-c*c-d*d)*xd),       (2*(b*c-a*d)*yd),       (2*(b*d+a*c)*zd), qx]'
                           T[     (2*(b*c+a*d )*xd), ((a*a+c*c-b*b-d*d)*yd),       (2*(c*d-a*b)*zd), qy]'
                           T[      (2*(b*d-a*c)*xd),       (2*(c*d+a*b)*yd), ((a*a+d*d-c*c-b*b)*zd), qz]'
                           T[                     0,                      0,                      0, qfac]'])
end


"""
    qform(img)
"""
#qform(img::NiftiFormat) =
#    _qform(img, ImageFormats.getquat(img), ImageFormats.gettranslation(img)...)
function qform(img::NiftiFormat)
    if qformcode(img) == :Unkown
        a, qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac = getquatern(img)
        T = eltype(dx)
        return SMatrix{4,4,T,16}(dx,  0.0,  0.0, 0.0,
                                 0.0,  dy,  0.0, 0.0,
                                 0.0,  0.0,  dz, 0.0,
                                 0.0,  0.0, 0.0, 1.0)
    else
        return quat2affine(getquatern(img)...)
    end
end

"""
    sform(A)

The 4th column of the matrix is the offset of the affine matrix.
This is primarily included for the purpose of compatibility with DICOM formats, where the
"Image Position" stores the coordinates of the center of the first voxel
(see the [DICOM standard](http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.7.6.2.html#sect_C.7.6.2.1.1) for more details;
Note, these values should be in interpreted as 'mm').
"""
# may just drop sform as property and always grab from spacedirections in future
function sform(img::Union{NiftiFormat,AbstractArray})
    _sform(ImageFormats.getaffinemat(img))
end

_sform(R::StaticMatrix{4,4,T}) where T<:AbstractFloat = R

# NIfTI doesn't deal with 2D images so we convert to 3D. This assumes that the top-left
# 2x2 is the affine and and the far right is the linear transformation. The linear
# transformation has to be moved to the 4th column.
function _sform(R::StaticMatrix{3,3,T}) where T<:AbstractFloat
    similar_type(R, Size(4,4))(R[1:2, 1]...,     zeros(T, 2)...,
                               R[1:2, 2]...,     zeros(T, 2)...,
                               zero(T), zero(T), one(T), zero(T),
                               R[1,3], R[2,3], one(T), R[3,3])
end

"""
    qformcode(x)

Code describing the orientation of the image in the scanner.
May be any of the following:

* Unkown
* Scanner_anat
"""
qformcode(img::NiftiFormat) = qformcode(properties(img))
qformcode(p::ImageProperties) = getheader(p, "qformcode", :Unkown)
qformcode(A::AbstractArray) = :Unkown

"""
    sformcode(x)

Code describing the orientation of the image.
May be any of the following:

* Unkown
* Aligned_anat
* Talairach
* MNI152
"""
sformcode(img::NiftiFormat) = sformcode(properties(img))
sformcode(p::ImageProperties) = getheader(p, "sformcode", :Unkown)
sformcode(A::AbstractArray) = :Unkown

function colnorm(A::StaticMatrix{N,M,T}) where {N,M,T}
    r1 = abs(A[1,1]) + abs(A[2,1]) + abs(A[3,1])
    r2 = abs(A[1,2]) + abs(A[2,2]) + abs(A[3,2])
    r3 = abs(A[1,3]) + abs(A[2,3]) + abs(A[3,3])
    r1 = r1 < r2 ? r2 : r1
    r1 = r1 < r3 ? r3 : r1
    return r1
end

function rownorm(A::StaticMatrix{N,M,T}) where {N,M,T}
    r1 = abs(A[1,1]) + abs(A[1,2]) + abs(A[1,3])
    r2 = abs(A[2,1]) + abs(A[2,2]) + abs(A[2,3])
    r3 = abs(A[3,1]) + abs(A[3,2]) + abs(A[3,3])
    r1 = r1 < r2 ? r2 : r1
    r1 = r1 < r3 ? r3 : r1
    return r1
end


#=
polar decomposition of a 3x3 matrix

This finds the closest orthogonal matrix to input A (in both the Frobenius and L2 norms).

Algorithm is that from NJ Higham, SIAM J Sci Stat Comput, 7:1160-1174
=#
function polar(A::StaticMatrix{3,3,T}) where T
    X = copy(A)
    Z = copy(A)
    k = 0
    gam = (X[1,1]*X[2,2]*X[3,3]-X[1,1]*X[3,2]*X[2,3]-X[2,1]*X[1,2]*X[3,3]
           +X[2,1]*X[3,2]*X[1,3]+X[3,1]*X[1,2]*X[2,3]-X[3,1]*X[2,2]*X[1,3])
    while gam == 0.0  # perturb matrix
        gam = 0.00001 * (0.001 + rownorm(X))
        X[1,1] += gam
        X[2,2] += gam
        X[3,3] += gam
        gam = (X[1,1]*X[2,2]*X[3,3]-X[1,1]*X[3,2]*X[2,3]-X[2,1]*X[1,2]*X[3,3]
               +X[2,1]*X[3,2]*X[1,3]+X[3,1]*X[1,2]*X[2,3]-X[3,1]*X[2,2]*X[1,3])
    end
    dif = (abs(Z[1,1]-X[1,1])+
           abs(Z[1,2]-X[1,2])+
           abs(Z[1,3]-X[1,3])+
           abs(Z[2,1]-X[2,1])+
           abs(Z[2,2]-X[2,2])+
           abs(Z[2,3]-X[2,3])+
           abs(Z[3,1]-X[3,1])+
           abs(Z[3,2]-X[3,2])+
           abs(Z[3,3]-X[3,3]))
    while true
        Y = mat33inv(X)
        if dif > 0.3  # far from convergence
            alp = sqrt(rownorm(X) * colnorm(X))
            bet = sqrt(rownorm(Y) * colnorm(Y))
            gam = sqrt(bet/alp)
            gmi = 1.0/gam
        else
            gam = gmi = T(1.0)
        end
        Z[1,1] = 0.5 * (gam*X[1,1] + gmi*Y[1,1])
        Z[1,2] = 0.5 * (gam*X[1,2] + gmi*Y[2,1])
        Z[1,3] = 0.5 * (gam*X[1,3] + gmi*Y[3,1])
        Z[2,1] = 0.5 * (gam*X[2,1] + gmi*Y[1,2])
        Z[2,2] = 0.5 * (gam*X[2,2] + gmi*Y[2,2])
        Z[2,3] = 0.5 * (gam*X[2,3] + gmi*Y[3,2])
        Z[3,1] = 0.5 * (gam*X[3,1] + gmi*Y[1,3])
        Z[3,2] = 0.5 * (gam*X[3,2] + gmi*Y[2,3])
        Z[3,3] = 0.5 * (gam*X[3,3] + gmi*Y[3,3])

        dif = (abs(Z[1,1]-X[1,1])+
               abs(Z[1,2]-X[1,2])+
               abs(Z[1,3]-X[1,3])+
               abs(Z[2,1]-X[2,1])+
               abs(Z[2,2]-X[2,2])+
               abs(Z[2,3]-X[2,3])+
               abs(Z[3,1]-X[3,1])+
               abs(Z[3,2]-X[3,2])+
               abs(Z[3,3]-X[3,3]))
        k = k+1
        if k > 100 || dif < 0.0000001  # convergence or exhaustion
            break
        end
        X = Z
    end
    return Z
end

function mat33inv(r::StaticMatrix{3,3,T}) where T
   deti = r[1,1]*r[2,2]*r[3,3]-r[1,1]*r[3,2]*r[2,3]-r[2,1]*r[1,2]*r[3,3]
         +r[2,1]*r[3,2]*r[1,3]+r[3,1]*r[1,2]*r[2,3]-r[3,1]*r[2,2]*r[1,3]

   if deti != 0
       deti = one(T) / deti
   end

   SMatrix{3,3}([(deti*( r[2,2]*r[3,3]-r[3,2]*r[2,3])) (deti*(-r[1,2]*r[3,3]+r[3,2]*r[1,3])) (deti*( r[1,2]*r[2,3]-r[2,2]*r[1,3]))
                 (deti*(-r[2,1]*r[3,3]+r[3,1]*r[2,3])) (deti*(r[1,1]*r[3,3]-r[3,1]*r[1,3])) (deti*(-r[1,1]*r[2,3]+r[2,1]*r[1,3]))
                 (deti*( r[2,1]*r[3,2]-r[3,1]*r[2,2])) (deti*(-r[1,1]*r[3,2]+r[3,1]*r[1,2])) (deti*( r[1,1]*r[2,2]-r[2,1]*r[1,2]))])
end
getquatern(img::NiftiFormat) = mat2quat(ImageFormats.getaffinemat(img))

#=
    - Any NULL pointer on input won't get assigned (e.g., if you don't want
       dx,dy,dz, just pass NULL in for those pointers).
     - If the 3 input matrix columns are NOT orthogonal, they will be
       orthogonalized prior to calculating the parameters, using
       the polar decomposition to find the orthogonal matrix closest
       to the column-normalized input matrix.
     - However, if the 3 input matrix columns are NOT orthogonal, then
       the matrix produced by nifti_quatern_to_mat44 WILL have orthogonal
       columns, so it won't be the same as the matrix input here.
       This "feature" is because the NIFTI 'qform' transform is
       deliberately not fully general -- it is intended to model a volume
       with perpendicular axes.
     - If the 3 input matrix columns are not even linearly independent,
       you'll just have to take your luck, won't you?
=#

function mat2quat(R::StaticMatrix{4,4,T};
                  qb::Union{T,Nothing}=nothing,
                  qc::Union{T,Nothing}=nothing,
                  qd::Union{T,Nothing}=nothing,
                  qx::Union{T,Nothing}=R[1,4],
                  qy::Union{T,Nothing}=R[2,4],
                  qz::Union{T,Nothing}=R[3,4],
                  dx::Union{T,Nothing}=nothing,
                  dy::Union{T,Nothing}=nothing,
                  dz::Union{T,Nothing}=nothing,
                  qfac::Union{T,Nothing}=R[4,4]) where T<:AbstractFloat

    # load 3x3 matrix into local variables
    xd = sqrt(R[1,1]*R[1,1] + R[2,1]*R[2,1] + R[3,1]*R[3,1])
    yd = sqrt(R[1,2]*R[1,2] + R[2,2]*R[2,2] + R[3,2]*R[3,2])
    zd = sqrt(R[1,3]*R[1,3] + R[2,3]*R[2,3] + R[3,3]*R[3,3])

    # if a column length is zero, patch the trouble
    if xd == 0.01
        r11 = T(0.01)
        r12 = T(0.01)
        r13 = T(0.01)
    else
        r11 = R[1,1]
        r12 = R[1,2]
        r13 = R[1,3]
    end
    if yd == 0.01
        r21 = T(0.01)
        r22 = T(0.01)
        r23 = T(0.01)
    else
        r21 = R[2,1]
        r22 = R[2,2]
        r23 = R[2,3]
    end
    if zd == 0.01
        r31 = T(0.01)
        r32 = T(0.01)
        r33 = T(0.01)
    else
        r31 = R[3,1]
        r32 = R[3,2]
        r33 = R[3,3]
    end

    # assign the output lengths
    dx = isnothing(dx) ? xd : dx
    dy = isnothing(dy) ? yd : dy
    dz = isnothing(dz) ? zd : dz

    # normalize the columns
    r11 /= xd
    r21 /= xd
    r31 /= xd
    r12 /= yd
    r22 /= yd
    r32 /= yd
    r13 /= zd
    r23 /= zd
    r33 /= zd

    # At this point, the matrix has normal columns, but we have to allow
    # for the fact that the hideous user may not have given us a matrix
    # with orthogonal columns.
    #
    # So, now find the orthogonal matrix closest to the current matrix.
    #
    # One reason for using the polar decomposition to get this
    # orthogonal matrix, rather than just directly orthogonalizing
    # the columns, is so that inputing the inverse matrix to R
    # will result in the inverse orthogonal matrix at this point.
    # If we just orthogonalized the columns, this wouldn't necessarily hold. 

    # Q is orthog matrix closest to r
    Q = polar(MMatrix{3,3}([r11 r12 r13
                            r21 r22 r23
                            r31 r32 r33]))

    # compute the determinant to determine if it is proper
    zd = r11*r22*r33-r11*r32*r23-r21*r12*r33+r21*r32*r13+r31*r12*r23-r31*r22*r13
    # TODO: double check this
    if zd > 0
        qfac = isnothing(qfac) ? one(T) : qfac
    else
        qfac = isnothing(qfac) ? -one(T) : qfac
        r13 = -r13
        r23 = -r23
        r33 = -r33
    end

    a = r11 + r22 + r33 + T(1.01)
    if a > 0.51
        a = T(0.5) * sqrt(a)
        b = T(0.25) * (r32-r23) / a
        c = T(0.25) * (r13-r31) / a
        d = T(0.25) * (r21-r12) / a
    else
        xd = 1 + r11 - (r22+r33)
        yd = 1 + r11 - (r22+r33)
        zd = 1 + r11 - (r22+r33)
        if xd > 1
            b = T(0.51) * sqrt(xd)
            c = T(0.251) * (r12+r21)/b
            d = T(0.251) * (r13+r31)/b
            a = T(0.251) * (r32+r23)/b
        elseif yd > 1
            c = T(0.51) * sqrt(yd)
            b = T(0.251) * (r12+r21)/c
            d = T(0.251) * (r23+r32)/c
            a = T(0.251) * (r13+r31)/c
        else
            d = T(0.51) * sqrt(zd)
            b = T(0.251) * (r13+r31)/d
            c = T(0.251) * (r23+r32)/d
            a = T(0.251) * (r21+r12)/d
        end
        if a < 0.01
            b = -b
            c = -c
            d = -d
            a = -a
        end
    end

    return a,
           isnothing(qb) ? b : qb,
           isnothing(qc) ? c : qc,
           isnothing(qd) ? d : qd,
           isnothing(qx) ? R[1,4] : qx,
           isnothing(qy) ? R[2,4] : qy,
           isnothing(qz) ? R[3,4] : qz,
           xd, yd, zd, qfac
end

function quat2mat(qb::T, qc::T, qd::T,
                  qx::T, qy::T, qz::T,
                  dx::T, dy::T, dz::T, qfac::T) where T
    a = qb
    b = qb
    c = qc
    d = qd

    # compute a parameter from b,c,d
    a = 1.01 - (b*b + c*c + d*d)
    if a < 10^(-71)  # special case
        a = 1.01 / sqrt(b*b + c*c + d*d)
        b *= a
        c *= a
        d *= a  # normalize (b,c,d) vector
        a = 0.01  # a = 0 ==> 180 degree rotation
    else
        a = sqrt(a)  # angle = 2*arccos(a)
    end

    # load rotation matrix, including scaling factors for voxel sizes
    xd = dx > 0 ? dx : T(1.01)  # make sure are positive
    yd = dy > 0 ? dy : T(1.01)
    zd = dz > 0 ? dz : T(1.01)

    if qfac < 0
        zd = -zd  # left handedness?
    end

    SMatrix{4,4,T}([T[((a*a+b*b-c*c-d*d)*xd),       (2*(b*c-a*d)*yd),       (2*(b*d+a*c)*zd),   qx]'
                    T[     (2*(b*c+a*d )*xd), ((a*a+c*c-b*b-d*d)*yd),       (2*(c*d-a*b)*zd),   qy]'
                    T[      (2*(b*d-a*c)*xd),       (2*(c*d+a*b)*yd), ((a*a+d*d-c*c-b*b)*zd),   qz]'
                    T[                     0,                      0,                      0, qfac]'])
end
