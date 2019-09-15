"""
    NiftiExtension

"`edata::Vector{UInt8}`: raw data with no byte-swapping"
"""
mutable struct NiftiExtension
    ecode::Symbol
    edata::Vector{UInt8}
end

# hi juliaString ctermfg orange
const NiftiEcode = Dict{Int,Symbol}(
    0  => :Ignore,
    2  => :DICOM,
    4  => :AFNI,
    6  => :Comment,
    8  => :XCEDE,
    10 => :JimDimInfo,
    12 => :WorkflowFWDS,
    14 => :Freesurfer,
    16 => :PyPickly,
    18 => :MindIdent,
    20 => :BValue,
    22 => :SphericalDirection,
    24 => :DTComponent,
    26 => :SHCDegreeOrder,
    28 => :Voxbo,
    30 => :Caret,
    32 => :CIfTI,
    34 => :VariableFrameTiming,
    38 => :Eval,
    40 => :Matlab)

const NiftiEcodeReverse = Dict{Symbol,Int}()
for (k, v) in NiftiEcode
    NiftiEcodeReverse[v] = k
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
ecode(x::Int32) = get(NiftiEcode, x, :Ignore)
ecode(x::NiftiExtension) = x.ecode
ecode(x::Vector{NiftiExtension}) = map(ecode, x)
ecode(img::ImageMeta) = ecode(extension(img))

# Calculates the size of a NIfTI extension
"""
    esize(ex::NiftiExtension)

NiftiExtensions should be of a byte size that is a mulitple of 16. This includes
raw encoding of the the `ecode` (as an Int32) and the esize itself (also as an
Int32). Therefore, `8 + sizeof(ex.edata)` should be divisible by 16.
"""
function esize(ex::NiftiExtension)
    ret = 8 + length(ex.edata)
    if ret%16 != 0
        error("NiftiExtension has innapropriate size. See docstrings for more details.")
    else
        return ret
    end
end


extension(x::Any) = getter(x, "extension", Vector{NiftiExtension}, NiftiExtension[])
extension!(x::Any, val::Vector{NiftiExtension}) = setter!(x, "extension", val, Vector{NiftiExtension})


function read(io::IO, p::AbstractDict, ::Type{NiftiExtension})
    ret = NiftiExtension[]
    if eof(io)
        return ret
    else
        ext = read!(io, Array{UInt8}(undef, 4))
        if ext[1] == 0
            return ret
        else
            counter = position(io)
            while counter < (data_offset(p)-1)
                esize = read(io, Int32)
                ec = read(io, Int32)
                push!(ret, NiftiExtension(ecode(ec), read!(io, Array{UInt8}(undef, esize-8))))
                counter += esize
            end
            return ret
        end
    end
end

function write(io::IO, x::Vector{NiftiExtension})
    if isempty(x)
        write(io, fill(UInt8(0), 4))
    else
        write(io, UInt8[1, 0, 0, 0])
        for ex in x
            write(io, Int32(esize(ex)))
            write(io, Int32(get(NiftiEcodeReverse, ecode(ex), 0)))
            write(io, ex.edata)
#            write(io, zeros(UInt8, esize(ex) - length(ex.edata)))
        end
    end
end
