function eltransform(to::ImageMeta, from::AbstractArray)
end

# something like StructOfArrays.jl may be better than these
expand2indices(A::AbstractArray{T,NA}, idx::Tuple{Vararg{<:Any,NI}}) where {T,NA,NI} =
    to_indices(A,(idx..., ntuple(i->Colon(), NA-NI)...))
expand2indices(A::AbstractArray{T,N}, idx::Tuple{Vararg{<:Any,N}}) where {T,N} =
    to_indices(A,idx)

AxType{Ax} = Union{AxisArray{T,N,D,Ax},ImageMeta{T,N,AxisArray{T,N,D,Ax}}} where {T,N,D}

# SMatrix
function transform(::Type{SMatrix{Tt}}, A::AbstractArray{Ta}) where {Tt,Ta}
    if Ta<:Tt
        _matview(A)
    elseif isabstracttype(Tt)
        _matview(A)
    else
        _matview(convert(Tt, A))
    end
end
_matview(A::AxType) = matview(PermutedDimsArray(A, permutation((:row, :col, setdiff(axisnames(A), (:row,:col))...,), axisnames(A))))
_matview(A::AxType{Tuple{Axis{:row}, Axis{:col}, Vararg{<:Axis}}}) = matview(A)
_matview(A::AbstractArray) = matview(A)
function matview(A::AbstractArray{T,N}) where {T,N}
    n = size(A,1)
    m = size(A,2)
    L = n*m
    reinterpret(SArray{Tuple{n,m},T,2,L}, reshape(A, (L, map(i->size(A,i),3:N)...,)))
end
extractmat(A::AbstractArray{SMatrix{R,C,T,L},N}) where {R,C,T,N,L} = reshape(reinterpret(T, A), (R,C, Base.tail(size(A))...,))

# SVector
function transform(::Type{SVector{Tt}}, A::AbstractArray{Ta}) where {Tt,Ta}
    if Ta<:Tt
        _vecview(A)
    elseif isabstracttype(Tt)
        _vecview(A)
    else
        _vecview(convert(Tt, A))
    end
end

_vecview(A::AxType) = vecview(PermutedDimsArray(A, permutation((:vecdim, setdiff(axisnames(A), :vecdim)...,), axisnames(A))))
_vecview(A::AxType{Tuple{Axis{:vecdim,<:Any}, Vararg{<:Axis}}}) = vecview(A)
_vecview(A::AbstractArray) = vecview(A)
vecview(A::AbstractArray{T,N}) where {T,N} = reinterpret(SVector{size(A,1),T},A)

# Point
function transform(::Type{Point{Tt}}, A::AbstractArray{Ta}) where {Tt,Ta}
    if Ta<:Tt
        _pointview(A)
    elseif isabstracttype(Tt)
        _pointview(A)
    else
        _pointview(convert(Tt, A))
    end
end
_pointview(A::AxType) = pointview(PermutedDimsArray(A, permutation((:pointdim, setdiff(axisnames(A), :pointdim)...,), axisnames(A))))
_pointview(A::AxType{Tuple{Axis{:pointdim,<:Any}, Vararg{<:Axis}}}) = pointview(A)
_pointview(A::AbstractArray) = pointview(A)
function pointview(A::AbstractArray)
    @assert size(A,1) == 3 "First dimension must be of size 3."
    mappedarray(Point, extractpoint,
                A[expand2indices(A,(1,))...],
                A[expand2indices(A,(2,))...],
                A[expand2indices(A,(3,))...])
end

# Quat
function transform(::Type{Quat{Tt}}, A::AbstractArray{Ta}) where {Tt,Ta}
    if Ta<:Tt
        _quatview(A)
    elseif isabstracttype(Tt)
        _quatview(A)
    else
        _quatview(convert(Tt, A))
    end
end
_quatview(A::AxType) = quatview(PermutedDimsArray(A, permutation((:quatdim, setdiff(axisnames(A), :quatdim)...,), axisnames(A))))
_quatview(A::AxType{Tuple{Axis{:quatdim,<:Any}, Vararg{<:Axis}}}) = quatview(A)
_quatview(A::AbstractArray) = quatview(A)
extractquat(a::Quat) = (a.w, a.x, a.y, a.z)
extractquat(A::AbstractArray{Quat}) = extractquat.(A)
function quatview(A::AbstractArray)
    if size(A, 1) != 4
        @error "Selected dimension, $i, must be of size 4"
    else
        return mappedarray(Quat, extractquat, A[expand2indices(A,(1,))...],
                           A[expand2indices(A,(2,))...], A[expand2indices(A,(3,))...],
                           A[expand2indices(A,(4,))...])
    end
end

# Distributions
for D in (TDist,Chi,Chisq,Poisson,FDist,Beta,Binomial,Gamma,Normal,NoncentralT,
          NoncentralChisq,Logistic,Uniform,NoncentralF,GeneralizedExtremeValue)
    fn = fieldnames(D)
    L = length(fn)

    @eval begin
        extractstat(d::$D) = map(i->getfield(d,i), $fn)
        extractstat(A::AbstractArray{$D}) = extractstat.(A)
        statview(::Type{$D}, A::AbstractArray{T,N}) where {T,N} =
            mappedarray($D, extractstat, ntuple(i->A[expand2indices(A,(i,))...], $L)...)
        @inline function Base.read(io::IO, ::Type{$(D){T}}) where {T}
            elements = Ref{NTuple{$L,T}}()
            read!(io, elements)
            $D(elements[])
        end
        @inline Base.write(io::IO, d::$D) = write(io, Ref(extractstat(d)))
    end
end
function transform(D::Type{<:Distribution}, A::AbstractArray{T}) where {T}
    Tt = eltype(D)
    if T<:Tt
        _statview(D, A)
    elseif isabstracttype(Tt)
        _statview(D, A)
    else
        _statview(D, convert(Tt, A))
    end
end
_statview(D::Type{<:Distribution}, A::AxType) = statview(D, PermutedDimsArray(A, permutation((:statdim, setdiff(axisnames(A), :statdim)...,), axisnames(A))))
_statview(D::Type{<:Distribution}, A::AxType{Tuple{Axis{:statdim,<:Any}, Vararg{<:Axis}}}) = statview(D, A)
_stattview(D::Type{<:Distribution}, A::AbstractArray) = statview(D, A)
