###
### Orientation utils
###
@inline function colnorm(A::StaticMatrix{3,3,Float64})
    r1 = abs(A[1,1]) + abs(A[2,1]) + abs(A[3,1])
    r2 = abs(A[1,2]) + abs(A[2,2]) + abs(A[3,2])
    r3 = abs(A[1,3]) + abs(A[2,3]) + abs(A[3,3])
    r1 = r1 < r2 ? r2 : r1
    r1 = r1 < r3 ? r3 : r1
    return r1
end

@inline function rownorm(A::StaticMatrix{3,3,Float64})
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
function polar(A::StaticMatrix{3,3,Float64})
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

function mat33inv(r::StaticMatrix{3,3,Float64})
   deti = r[1,1]*r[2,2]*r[3,3]-r[1,1]*r[3,2]*r[2,3]-r[2,1]*r[1,2]*r[3,3]
         +r[2,1]*r[3,2]*r[1,3]+r[3,1]*r[1,2]*r[2,3]-r[3,1]*r[2,2]*r[1,3]

   if deti != 0
       deti = one(T) / deti
   end

   SMatrix{3,3}([(deti*( r[2,2]*r[3,3]-r[3,2]*r[2,3])) (deti*(-r[1,2]*r[3,3]+r[3,2]*r[1,3])) (deti*( r[1,2]*r[2,3]-r[2,2]*r[1,3]))
                 (deti*(-r[2,1]*r[3,3]+r[3,1]*r[2,3])) (deti*(r[1,1]*r[3,3]-r[3,1]*r[1,3])) (deti*(-r[1,1]*r[2,3]+r[2,1]*r[1,3]))
                 (deti*( r[2,1]*r[3,2]-r[3,1]*r[2,2])) (deti*(-r[1,1]*r[3,2]+r[3,1]*r[1,2])) (deti*( r[1,1]*r[2,2]-r[2,1]*r[1,2]))])
end
#getquatern(img::NiftiFormat) = mat2quat(getaffinemat(img))
