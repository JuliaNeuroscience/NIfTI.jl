function colnorm(A::StaticMatrix{3,3,Float64})
    r1 = abs(A[1,1]) + abs(A[2,1]) + abs(A[3,1])
    r2 = abs(A[1,2]) + abs(A[2,2]) + abs(A[3,2])
    r3 = abs(A[1,3]) + abs(A[2,3]) + abs(A[3,3])
    r1 = r1 < r2 ? r2 : r1
    r1 = r1 < r3 ? r3 : r1
    return r1
end

function rownorm(A::StaticMatrix{3,3,Float64})
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
function polar(A::MMatrix{3,3,Float64})
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

function mat33inv(r::MMatrix{3,3,T}) where T<:Union{Float64,Float32}
   deti = r[1,1]*r[2,2]*r[3,3]-r[1,1]*r[3,2]*r[2,3]-r[2,1]*r[1,2]*r[3,3]
         +r[2,1]*r[3,2]*r[1,3]+r[3,1]*r[1,2]*r[2,3]-r[3,1]*r[2,2]*r[1,3]

   if deti != 0
       deti = one(T) / deti
   end

   MMatrix{3,3}([(deti*( r[2,2]*r[3,3]-r[3,2]*r[2,3])) (deti*(-r[1,2]*r[3,3]+r[3,2]*r[1,3])) (deti*( r[1,2]*r[2,3]-r[2,2]*r[1,3]))
                 (deti*(-r[2,1]*r[3,3]+r[3,1]*r[2,3])) (deti*(r[1,1]*r[3,3]-r[3,1]*r[1,3])) (deti*(-r[1,1]*r[2,3]+r[2,1]*r[1,3]))
                 (deti*( r[2,1]*r[3,2]-r[3,1]*r[2,2])) (deti*(-r[1,1]*r[3,2]+r[3,1]*r[1,2])) (deti*( r[1,1]*r[2,2]-r[2,1]*r[1,2]))])
end

function mat2quat(R::MMatrix{4,4,Float64};
                  qb::Union{Float64,Nothing}=quaternb(x),
                  qc::Union{Float64,Nothing}=quaternc(x),
                  qd::Union{Float64,Nothing}=quaternd(x),
                  qx::Union{Float64,Nothing}=R[1,4],
                  qy::Union{Float64,Nothing}=R[2,4],
                  qz::Union{Float64,Nothing}=R[3,4],
                  dx::Union{Float64,Nothing}=nothing,
                  dy::Union{Float64,Nothing}=nothing,
                  dz::Union{Float64,Nothing}=nothing,
                  qfac::Union{Float64,Nothing}=R[4,4])

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

