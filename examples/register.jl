# register.jl
# Register NIfTI volumes based on header information

using NIfTI, Grid

# Affine matrices that add or subtract one
# Needed because Julia and Grid.jl assume zero-based indexing
const ADD1 = [1 0 0 1
              0 1 0 1
              0 0 1 1
              0 0 0 1]
const SUB1 = [1 0 0 -1
              0 1 0 -1
              0 0 1 -1
              0 0 0 1]

# Round to an integer if the second argument is an integer type
mayberound{T<:Integer}(x, ::Type{T}) = round(T, x)
mayberound(x, ::DataType) = x

function apply_rotation!(out, A, data)
    A = ADD1*A*SUB1

    # Get pointers to each pixel dim of data
    s = stride(data, 4)
    ptrs = [pointer_to_array(pointer(data)+(i-1)*s, (size(data, 1), size(data, 2), size(data, 3)))
            for i = 1:size(data, 4)]

    # Allocate B, which represents voxels in targ, and C, which
    # represents voxels in mov
    B = Array(Float32, 4)
    B[4] = 1
    C = Array(Float32, 4)
    Csub = pointer_to_array(pointer(C), 3)

    ic = InterpGridCoefs(ptrs[1], InterpLinear)
    for i = 1:size(out, 3)
        B[3] = i
        for j = 1:size(out, 2)
            B[2] = j
            for k = 1:size(out, 1)
                B[1] = k
                # Convert voxels in targ coords to voxels in mov coords
                A_mul_B!(C, A, B)

                # Interpolate voxels in targ coords from voxels in mov
                # coords
                set_position(ic, BCnan, false, Csub)
                v = interp(ic, ptrs[1])
                if !isnan(v)
                    out[k, j, i, 1] = mayberound(v, eltype(out))
                    for l = 2:size(data, 4)
                        out[k, j, i, l] = mayberound(interp(ic, ptrs[l]), eltype(out))
                    end
                end
            end
        end
    end
end

function register{T<:FloatingPoint}(targ::NIVolume, mov::NIVolume{T}, outtype::DataType=T)
    # Construct the affine matrix that maps voxels of targ to voxels of
    # mov
    targaffine = getaffine(targ.header)
    movaffine = getaffine(mov.header)
    A = inv(movaffine)*targaffine

    # Allocate output
    out = zeros(outtype, size(targ, 1), size(targ, 2), size(targ, 3), size(mov, 4))

    # Perform rotation
    apply_rotation!(out, A, mov.raw)

    # Set header of mov
    newheader = deepcopy(mov.header)
    setaffine(newheader, targaffine)
    newheader.pixdim = (newheader.pixdim[1], targ.header.pixdim[2:4]..., newheader.pixdim[5:end]...)
    NIVolume(newheader, mov.extensions, out)
end
register(targ::NIVolume, mov::NIVolume) =
    register(targ, NIVolume(mov.header, mov.extensions, map(Float32, mov.raw)), eltype(mov))
register(targ::AbstractString, mov::AbstractString, out::AbstractString) =
    niwrite(out, register(niread(targ, mmap=true), niread(mov, mmap=true)))
