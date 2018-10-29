abstract type AbstractVectorImage end


struct DisplacementVectorImage <: AbstractVectorImage end
struct GeneralVectorImage <: AbstractVectorImage end

"""
These images contain a dimension for storing additional per voxel/pixel data
"""
struct VectorImage <: AbstractVectorImage end

"""
To store an M x N matrix at each voxel:
- dataset must have a 5th dimension (dim[0]=5 and dim[5]>1)
- dim[5] must be M*N
- intent_p1 must be M (in float format)
- intent_p2 must be N (ditto)
- the matrix values A[i][[j] are stored in row-order:
    - A[0][0] A[0][1] ... A[0][N-1]
    - A[1][0] A[1][1] ... A[1][N-1]
    - etc., until
    - A[M-1][0] A[M-1][1] ... A[M-1][N-1]
"""
struct GeneralMatrixImage <: AbstractVectorImage
    M::Integer
    N::Integer
end

"""
To store an NxN symmetric matrix at each voxel:
- dataset must have a 5th dimension
- dim[5] must be N*(N+1)/2
- intent_p1 must be N (in float format)
- the matrix values A[i][[j] are stored in row-order:
    - A[0][0]
    - A[1][0] A[1][1]
    - A[2][0] A[2][1] A[2][2]
    - etc.: row-by-row
"""
struct SymmetricMatrixImage <: AbstractVectorImage
    N::AbstractFloat
end

"""
To signify that the vector value at each voxel is really a
spatial coordinate (e.g., the vertices or nodes of a surface mesh):
- dataset must have a 5th dimension
- dim[1] = 5
- dim[2] = number of points
- dim[3] = dim[4] = dim[5] = 1
- dim[6] must be the dimensionality of space (e.g., 3 => 3D space).
- intent_name may describe the object these points come from
 (e.g., "pial", "gray/white" , "EEG", "MEG").
"""
struct PointSetImage <: AbstractVectorImage
    npoints::Integer
    sdim::Integer
end

"""
To signify that the vector value at each voxel is really a triple
of indexes (e.g., forming a triangle) from a pointset dataset:
- dataset must have a 5th dimension
- dim[1] = 5
- dim[2] = number of triangles
- dim[3] = dim[4] = dim[5] = 1
- dim[6] = 3
- datatype should be an integer type (preferably DT_INT32)
- the data values are indexes (0,1,...) into a pointset dataset.
"""
struct TriangleImage <: AbstractVectorImage
    ntriangles::Integer
end

function intent_vector(intent_code::AbstractString, img::ImageMeta, intent_p1::T, intent_p2::T, intent_p3::T) where {T<:AbstractFloat}
    if intent_code == "GeneralMatrix"
        if intent_p1 * intent_p2 == ndims(vol)[5]
            img.properties["header"]["intent"] = GeneralMatrixImage(intent_p1, intent_p2)
        else
            error("GeneralMatrixImage must have:\n
                  dim[5] == M*N (where M = intent_p1, N = intent_p2)")
        end
    elseif intent_code == "SymmetricMatrix"
        if intent_p1 * (intent_p1+1)/2 == ndims(vol)[5]
            img.properties["header"]["intent"] = SymmetricMatrixImage(intent_p1)
        else
            error("SymmetricMatrixImage must have dim[5] == N*(N+1)/2
                   (where N=intent_p1)")
        end
    elseif intent_code == "PontSet"
        if dim[1] == 5 && dim[3] == 1 && dim[4] == 1 && dim[5] == 1
            img.properties["header"]["intent"] = PointSetImage(dim[2], dim[6])
        else
            error("PointSetImage must have:\n
                   dim[1] == 5\n
                   dim[3] == 1\n
                   dim[4] == 1\n
                   dim[5] == 1\n")
        end
    elseif intent_code == "GeneralVector"
        img.properties["header"]["intent"] = GeneralVectorImage()
    elseif intent_code == "TriangleSet"
        if dim[1] == 5 && dim[3] == 1 && dim[4] == 1 && dim[5] == 1 && dim[6] == 3
            img.properties["header"]["intent"] = TriangleImage(dim[2])
        else
             error("TriangleImage must have:\n
                   dim[1] == 5\n
                   dim[3] == 1\n
                   dim[4] == 1\n
                   dim[5] == 1\n
                   dim[6] == 3")
        end
    end
end

