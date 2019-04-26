"""
    NiftiExtension

"`edata::Vector{UInt8}`: raw data with no byte-swapping"
"""
mutable struct NiftiExtension
    ecode::Symbol
    edata::Vector{UInt8}
end

const NiftiEcode = Dict{Int16, Symbol}([
    (Int16(0), :Ignore),
    (Int16(2), :DICOM),
    (Int16(4), :AFNI),
    (Int16(6), :Comment),
    (Int16(8), :XCEDE),
    (Int16(10), :JimDimInfo),
    (Int16(12), :WorkflowFWDS),
    (Int16(14), :Freesurfer),
    (Int16(16), :PyPickly),
    (Int16(18), :MindIdent),
    (Int16(20), :BValue),
    (Int16(22), :SphericalDirection),
    (Int16(24), :DTComponent),
    (Int16(26), :SHCDegreeOrder),
    (Int16(28), :Voxbo),
    (Int16(30), :Caret),
    (Int16(32), :CIfTI),
    (Int16(34), :VariableFrameTiming),
    (Int16(38), :Eval),
    (Int16(40), :Matlab)
])

const NiftiEcodeReverse = Dict{Symbol,Int16}()
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
function ecode(img::ImageMeta)
    if haskey(img.properties, "nifti")
        return ecode(get(img["nifti"], "extension", extension(data(img))))
    else
        return ecode(extension(data(img)))
    end
end


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


# makes empty extension
function extension(img::ImageMeta)
    if haskey(img.properties, "nifti")
        return get(img["nifti"], "extension", extension(data(img)))
    else
        return extension(data(img))
    end
end
extension(A::AbstractArray) = NiftiExtension[]

function read(io::IOMeta{F}, ::Type{NiftiExtension}) where F
    ret = NiftiExtension[]
    if eof(io)
        return ret
    else
        extension = read!(io, Array{UInt8}(undef, 4))
        if extension[1] == 0
            return ret
        else
            counter = position(io)
            while counter < (data_offset(io)-1)
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
            write(io, Int32(esize(x)))
            write(io, Int32(get(NiftiEcodeReverse, ecode(x), 0)))
            write(io, x.edata)
            write(io, zeros(UInt8, esize(ex) - length(x.edata)))
        end
    end
end
