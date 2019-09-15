"""
# NIfTI Intent

## Distribution Images

Correlation
TDist
FDist
ZScore
Chisq
Beta
Binomial
Gamma
Poisson
Normal
NoncentralF
NoncentralChisq
Logistic
Laplace
Uniform
NoncentralT
Weibull
Chi
InverseGaussian
GeneralizedExtremeValue
PValue
LogPValue
Log10PValue
ParametricEstimate

## Label Images

* CustomLabels: element is an index into some set of labels. If present, the
  filename with the labels may be found in with `auxfile(img)`
* NeurNames: This is not yet fully implemented in the NIfTI.jl package or the
  NIfTI standard.

## Vector Images

* General Matrix: To store an M x N StaticMatrix at each voxel.
* Symmetric Matrix: To store an NxN symmetric matrix at each voxel:
    - dataset must have a 5th dimension
    - size(hdr, 5) must be N*(N+1)/2
    - intent_p1 must be N (in float format)
    - the matrix values A[i][[j] are stored in row-order:
        - A[1][1]
        - A[2][1] A[2][2]
        - A[3][1] A[3][2] A[3][3]
        - etc.: row-by-row
* Displacement Vector: Each element is a StaticVector representing displacement.
  The position of each element in the StaticVector is the dimension of
  displacement (e.g., 3 for spatial displacement, 2 for in-plane).
* General Vector: These images contain a StaticVector at each element of the array.
* Point: To signify that the vector value at each voxel is really a spatial
  coordinate (e.g., the vertices or nodes of a surface mesh):
    - intent_name may describe the object these points come from (e.g., "pial",
      "gray/white" , "EEG", "MEG").
* Triangle: To signify that the vector value at each voxel is really a triple of
  indexes (e.g., forming a triangle) from a pointset dataset.
* Quaternion: To signify that the vector value at each voxel is a quaternion:
   - datatype should be a floating point type
* Dimensionless: The name of the parameter may be stored in `intentname`

 Node Index: 
* Shape: To signify that the value at each location is a shape value, such as
  the curvature. Specified for the GIfTI format. loaded as a StaticVector element.
* Dimensionless Intent: 

## GIfTI Images

* GiftiTimeSeries: To signify that the value at each location is from a time series.
* GiftiNodeIndex: To signify that the value at each location is a node index, from a
  complete surface dataset. Specified for the GIfTI format. Loaded as a Point element.
* GiftiRGB
* GiftiRGBA
* GiftiShape


"""
extracttriangle(a::Triangle) = (a.data...,)
extractpoint(a::Point) = (a.data...,)


# TODO: check extractquat for proper output shape
extractquat(a::Quat) = (a.w, a.x, a.y, a.z)
extractquat(A::AbstractArray{Quat}) = extractquat.(A)

extractmat(A::AbstractArray{SMatrix{R,C,T,L},N}) where {R,C,T,N,L} =
    reshape(reinterpret(T, A), (R,C, Base.tail(size(A))...,))
function matview(A::AbstractArray{T,N}) where {T,N}
    n = size(A,1)
    m = size(A,2)
    L = n*m
    reinterpret(SArray{Tuple{n,m},T,2,L}, reshape(A, (L, map(i->size(A,i),3:N)...,)))
end

# Array of Vectors
vecview(A::AbstractArray{T,N}) where {T,N} = reinterpret(SVector{size(A,1),T},A)
extractvec(A::AbstractArray{SVector{L,T}}) where {L,T} = reshape(reinterpret(T, A), (L, Base.tail(size(A))...,))

function eltransform()
end

#############
# transform #
#############
function transform(A::AbstractArray, intent::Type{T}) where {T<:BitTypes}
    if isbitstype(T)
        if eltype(A) != T
            return reinterpret(T, A)
        else
            return A
        end
    end
end
# Distribution #
# TODO handle extractstat
transform(A::AbstractArray{T,5}, D::Type{<:StatP1}) where T = mappedarray(D, A)
transform(A::AbstractArray{T,5}, D::Type{<:StatP2}) where T = mappedarray(D, A[:,:,:,:,1], A[:,:,:,:,2])
transform(A::AbstractArray{T,5}, D::Type{<:StatP3}) where T = mappedarray(D, A[:,:,:,:,1], A[:,:,:,:,2], A[:,:,:,:,3])
# if no 5th dim then intent_params compose single distribution underlying image
# this is accomplished by 'parse_type_intent'
transform(A::AbstractArray{T,N}, D::Type{<:Distribution})  where {T,N} = A

# Labels #
# TODO: Probably should implement IndirectArrays
transform(A::AbstractArray, ::Type{CustomLabels}) = transform(A, eltype(A))
transform(A::AbstractArray, ::Type{NeuroNames}) = transform(A, eltype(A))

# Vector #
transform(A::AbstractArray, ::Type{<:SMatrix}) = matview(PermutedDimsArray(A, (1,2,3,4,6,5)))
transform(A::AbstractArray, ::Type{<:SVector}) = vecview(PermutedDimsArray(A, (5,1,2,3,4)))
transform(A::AbstractArray, ::Type{Quat}) = mappedarray(Quat, a->(a.w, a.x, a.y, a.z), A[:,:,:,:,1],A[:,:,:,:,2], A[:,:,:,:,3], A[:,:,:,:,4])
transform(A::AbstractArray, ::Type{Point}) = mappedarray(Point, A[:,:,:,:,1],A[:,:,:,:,2], A[:,:,:,:,3])
transform(A::AbstractArray, ::Type{Triangle}) = mappedarray(Triangle, A[:,:,:,:,1],A[:,:,:,:,2], A[:,:,:,:,3])
transform(A::AbstractArray, ::Type{Dimensionless}) = transform(A, eltype(A))
transform(A::AbstractArray, ::Type{DisplacementVector}) = transform(A, SVector)
# TODO: transform(A::AbstractArray, ::Type{SymmetricMatrix})
transform(A::AbstractArray, ::Type{GiftiNodeIndex}) = transform(A, Point)
transform(A::AbstractArray, ::Type{GiftiTimeSeries}) = transform(A, eltype(A))

# Gifti #
transform(A::AbstractArray, ::Type{GiftiShape}) =
    vecview(PermutedDimsArray(A, (5,1,2,3,4)))
transform(A::AbstractArray, ::Type{GiftiRGB}) =
    colorview(RGB, PermutedDimsArray(A, (5,1,2,3,4)))
transform(A::AbstractArray, ::Type{GiftiRGBA}) =
    colorview(RGBA, PermutedDimsArray(A, (5,1,2,3,4)))
transform(A::AbstractArray, ::Type{RGB}) = colorview(RGB, A)
transform(A::AbstractArray, ::Type{RGBA}) = colorview(RGBA, A)

##############
# rtransform #
##############
dim1to5(A::AbstractArray{T,7}) where T = reshape(PermutedDimsArray(A, (2,3,4,5,1,6,7)), )
dim1to5(A::AbstractArray{T,6}) where T = PermutedDimsArray(A, (2,3,4,5,1,6))
function dim1to5(A::AbstractArray{T,5}) where T
    sz = size(A)
    reshape(PermutedDimsArray(A, (2,3,4,5,1)), sz[2],sz[3],sz[4],sz[5],sz[1])
end
function dim1to5(A::AbstractArray{T,4}) where T
    sz = size(A)
    reshape(PermutedDimsArray(A, (2,3,4,1)), sz[2],sz[3],sz[4],1,sz[1])
end
function dim1to5(A::AbstractArray{T,3}) where T
    sz = size(A)
    reshape(PermutedDimsArray(A, (2,3,1)), sz[2],sz[3],1,1,sz[1])
end
function dim1to5(A::AbstractArray{T,2}) where T
    sz = size(A)
    reshape(PermutedDimsArray(A, (2,3,1)), sz[2],1,1,1,sz[1])
end

rtransform(A::AbstractArray{<:BitTypes}) = A
# Distribution #
# TODO handle extractstat
rtransform(A::AbstractArray{T}) where {T<:Distribution} = (T, dim1to5(extractstat(A)))

# Labels #
# TODO: Probably should implement IndirectArrays
#transform(A::AbstractArray, ::Type{CustomLabels}) = transform(A, eltype(A))
#transform(A::AbstractArray, ::Type{NeuroNames}) = transform(A, eltype(A))

# Vector #
rtransform(A::AbstractArray{T}) where {T<:SMatrix} = (T,dim1to5(extractmat(A)))
rtransform(A::AbstractArray{T}) where {T<:SVector} = (T,dim1to5(extractvec(A)))
rtransform(A::AbstractArray{T}) where {T<:Quat} = (T,dim1to5(extractquat(A)))
rtransform(A::AbstractArray{T}) where {T<:Point} = (T,dim1to5(extractpoint(A)))
rtransform(A::AbstractArray{T}) where {T<:Triangle} = (T,dim1to5(extracttriangle(A)))
# TODO: transform(A::AbstractArray, ::Type{DisplacementVector}) = transform(A, SVector)
# TODO: transform(A::AbstractArray, ::Type{Dimensionless}) = transform(A, eltype(A))
# TODO: transform(A::AbstractArray, ::Type{SymmetricMatrix})

# Gifti #
# TODO: transform(A::AbstractArray, ::Type{GiftiShape}) = vecview(PermutedDimsArray(A, (5,1,2,3,4))
# TODO: transform(A::AbstractArray, ::Type{GiftiNodeIndex}) = transform(A, Point)
# TODO: transform(A::AbstractArray, ::Type{GiftiTimeSeries}) = transform(A, eltype(A))
rtransform(A::AbstractArray{RGB,1}) = (GiftiRGB, dim1to5(colorview(RGB, A)))
rtransform(A::AbstractArray{RGBA,1}) = (GiftiRGBA, dim1to5(colorview(RGBA, A)))
rtransform(A::AbstractArray{RGB,N}) where N = (RGBA, reshape(channelview(A), size(A,1)*3, Base.tail(size(A))...))
rtransform(A::AbstractArray{RGBA,N}) where N = (RGBA, reshape(channelview(A), size(A,1)*4, Base.tail(size(A))...))
