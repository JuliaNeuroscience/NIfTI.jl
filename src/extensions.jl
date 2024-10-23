
@inline function to_ecode(x::Int32)
    if x === Int32(2)
        return :DICOM
    elseif x === Int32(4)
        return :AFNI
    elseif x === Int32(6)
        return :Comment
    elseif x === Int32(8)
        return :XCEDE
    elseif x === Int32(10)
        return :JimDimInfo
    elseif x === Int32(12)
        return :WorkflowFWDS
    elseif x === Int32(14)
        return :Freesurfer
    elseif x === Int32(16)
        return :PyPickly
    elseif x === Int32(18)
        return :MindIdent
    elseif x === Int32(20)
        return :BValue
    elseif x === Int32(22)
        return :SphericalDirection
    elseif x === Int32(24)
        return :DTComponent
    elseif x === Int32(26)
        return :SHCDegreeOrder
    elseif x === Int32(28)
        return :Voxbo
    elseif x === Int32(30)
        return :Caret
    elseif x === Int32(32)
        return :CIfTI
    elseif x === Int32(34)
        return :VariableFrameTiming
    elseif x === Int32(38)
        return :Eval
    elseif x === Int32(40)
        return :Matlab
    else
        return :Ignore
    end
end

@inline function ecode_to_int(x::Symbol)
    if x === :DICOM
        return Int32(2)
    elseif x === :AFNI
        return Int32(4)
    elseif x === :Comment
        return Int32(6)
    elseif x === :XCEDE
        return Int32(8)
    elseif x === :JimDimInfo
        return Int32(10)
    elseif x === :WorkflowFWDS
        return Int32(12)
    elseif x === :Freesurfer
        return Int32(14)
    elseif x === :PyPickly
        return Int32(16)
    elseif x === :MindIdent
        return Int32(18)
    elseif x === :BValue
        return Int32(20)
    elseif x === :SphericalDirection
        return Int32(22)
    elseif x === :DTComponent
        return Int32(24)
    elseif x === :SHCDegreeOrder
        return Int32(26)
    elseif x === :Voxbo
        return Int32(28)
    elseif x === :Caret
        return Int32(30)
    elseif x === :CIfTI
        return Int32(32)
    elseif x === :VariableFrameTiming
        return Int32(34)
    elseif x === :Eval
        return Int32(38)
    elseif x === :Matlab
        return Int32(40)
    else
        return Int32(0)
    end
end

"""
    NIfTIExtension

* `ecode::Int32`: extension code that can be interpreted used `ecode(::NIfTIExtension)`
* `edata::Vector{UInt8}`: raw data with no byte-swapping
"""
struct NIfTIExtension
    ecode::Int32
    edata::Vector{UInt8}
end

"""
    NIfTI Extension Codes

* Ignore
* DICOM: intended for raw DICOM attributes
* AFNI: Robert W Cox: rwcox@nih.gov https://afni.nimh.nih.gov/afni
* Comment: plain ASCII text only
* XCEDE: http://www.nbirn.net/Resources/Users/Applications/xcede/index.html
* JIMDimInfo: Dimensionalinformation for the JIM software (XML format)
* WorkflowFWDS:
* Freesurfer:
* PyPickle: embedded Python objects
* MiNDIdent: LONI MiND codes: http://www.loni.ucla.edu/twiki/bin/view/Main/MiND
* BValue
* SphericalDirection
* DTComponent
* SHCDegreeOrder
* Voxbo: www.voxbo.org
* Caret: http://brainvis.wustl.edu/wiki/index.php/Caret:Documentation
* CIfTI
* VariableFrameTiming
* AgilentProcpar
* Eval
* Matlab
"""
ecode(x::NIfTIExtension) = to_ecode(x.ecode)
ecode(x::Vector{NIfTIExtension}) = map(ecode, x)


# Calculates the size of a NIfTI extension
"""
    esize(ex::NIfTIExtension)

NIfTIExtensions should be of a byte size that is a mulitple of 16. This includes
raw encoding of the the `ecode` (as an Int32) and the esize itself (also as an
Int32). Therefore, `8 + sizeof(ex.edata)` should be divisible by 16.
"""
function esize(ex::NIfTIExtension)
    ret = 8 + length(ex.edata)
    @assert ret%16 != 0 "NIfTIExtension has innapropriate size. See docstrings for more details."
    return ret
end

function read_extensions(io, n; swapbyte::Bool=false, skip_extension_flag::Bool=false)
    ret = NIfTIExtension[]
    if eof(io)
        return ret
    end

    if !skip_extension_flag
        b1 = read(io, UInt8)
        # GZIP doesn't skip so we read and throw away
        read(io, UInt8)
        read(io, UInt8)
        read(io, UInt8)
        if b1 === zero(UInt8)
            return ret
        end
    end

    swap_int32 = swapbyte ? bswap : x -> x
    counter = skip_extension_flag ? 0 : 4
    while counter < (n - 1)
        esize = swap_int32(read(io, Int32))
        ec = swap_int32(read(io, Int32))
        push!(ret, NIfTIExtension(ec, read!(io, Array{UInt8}(undef, esize - 8))))
        counter += esize
    end
    return ret

end

write(path::AbstractString, x::Vector{NIfTIExtension}; swapbyte::Bool=false, skip_extension_flag::Bool=false) =
open(path, "w") do io
    write(io, x, swapbyte=swapbyte, skip_extension_flag=skip_extension_flag)
end

function write(io::IO, x::Vector{NIfTIExtension}; swapbyte::Bool=false, skip_extension_flag::Bool=false)
    if isempty(x)
        if !skip_extension_flag
            write(io, fill(UInt8(0), 4))
        end
        return
    end
    if !skip_extension_flag
        write(io, UInt8[1, 0, 0, 0])
    end

    swap_int32 = swapbyte ? bswap : x -> x

    for ex in x
        write(io, swap_int32(Int32(esize(ex))))
        write(io, swap_int32(ex.ecode))
        write(io, ex.edata)
        #       write(io, zeros(UInt8, esize(ex) - length(ex.edata)))
    end

end

