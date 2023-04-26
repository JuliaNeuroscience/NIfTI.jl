
"""
    NIfTI.Extension(code::Int32, data::Vector{UInt8})

* `ecode::Int32`: extension code that can be interpreted used `ecode(::NIfTI.Extension)`
* `edata::Vector{UInt8}`: raw data with no byte-swapping
"""
struct Extension
    ecode::Int32
    edata::Vector{UInt8}

    function Extension(code::Int32, data::Vector{UInt8})
        @assert (8 + length(data))%16 != 0 "NIfTI.Extension has innapropriate size. See docstrings for more details."
        new(code, data)
    end
    function Extension(x, data::Vector{UInt8})
        @assert (8 + length(data))%16 != 0 "NIfTI.Extension has innapropriate size. See docstrings for more details."
        if x === :DICOM
            code = Int32(2)
        elseif x === :AFNI
            code = Int32(4)
        elseif x === :Comment
            code = Int32(6)
        elseif x === :XCEDE
            code = Int32(8)
        elseif x === :JimDimInfo
            code = Int32(10)
        elseif x === :WorkflowFWDS
            code = Int32(12)
        elseif x === :Freesurfer
            code = Int32(14)
        elseif x === :PyPickly
            code = Int32(16)
        elseif x === :MindIdent
            code = Int32(18)
        elseif x === :BValue
            code = Int32(20)
        elseif x === :SphericalDirection
            code = Int32(22)
        elseif x === :DTComponent
            code = Int32(24)
        elseif x === :SHCDegreeOrder
            code = Int32(26)
        elseif x === :Voxbo
            code = Int32(28)
        elseif x === :Caret
            code = Int32(30)
        elseif x === :CIfTI
            code = Int32(32)
        elseif x === :VariableFrameTiming
            code = Int32(34)
        elseif x === :Eval
            code = Int32(38)
        elseif x === :Matlab
            code = Int32(40)
        else
            code = Int32(0)
        end
        new(code, data)
    end
end

"""
    NIfTI.ecode(e::NIfTI.Extension) -> Symbol
    NIfTI.ecode(e::Vector{NIfTI.Extension}) -> Vector{Symbol}

Converts stored `ecode` fields from `NIfTI.Extension` to there associated semantic labels.

NIfTI Extension Code Labels
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
@inline function ecode(e::Extension)
    x = e.ecode
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
ecode(x::Vector{Extension}) = map(ecode, x)

# Calculates the size of a NIfTI extension
"""
    NIfTI.esize(ex::NIfTI.Extension)

NIfTI.Extension should be of a byte size that is a mulitple of 16. This includes
raw encoding of the the `ecode` (as an Int32) and the esize itself (also as an
Int32). Therefore, `8 + sizeof(ex.edata)` should be divisible by 16.
"""
esize(ex::Extension) = length(ex.edata) + 8

function write(io::IO, x::Vector{Extension})
    if isempty(x)
        write(io, fill(UInt8(0), 4))
    else
        write(io, UInt8[1, 0, 0, 0])
        for ex in x
            write(io, Int32(esize(ex)))
            write(io, ex.ecode)
            write(io, ex.edata)
#            write(io, zeros(UInt8, esize(ex) - length(ex.edata)))
        end
    end
end

"""
    NIfTI.extensions(img::AbstractArray) -> Vector{NIfTI.Extension}

Returns a vector of `NIfTI.Extension`s associated with `img`. If no extensions are found,
an empty vector of `NIfTI.Extension`s is returned.
"""
function extensions(x::MetadataArray)
    exts = get(metadata(x), :extensions, nothing)
    return exts isa Vector{NIfTI.Extension} ? exts : NIfTI.Extension[]
end
extensions(x::AbstractArray) = NIfTI.Extension[]

