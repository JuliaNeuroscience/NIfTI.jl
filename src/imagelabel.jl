"""
Label Images

If intent_code is "Label" the value at each voxel is an index into some set of
labels. The filename with the labels may be stored in the `aux_file` field of
the NIfTI header. If `aux_file` is found an attempt is made to read it in using
FileIO. Arguments for reading these are passed via `args...` and `kwargs` from
the original call to `load` for a NIfTI file.

If intent_code is "NeuroName" value at each voxel is an index into the
[NeuroNames](http://braininfo.rprc.washington.edu/Nnont.aspx) label set. This
feature is not officially supported yet.


"""
abstract type AbstractLabelImage end

struct LabelImage <: AbstractLabelImage
    labels::Dict{Integer,String}
end

function write_labelimage()
end

function LabelImage(img::ImageMeta, hdr::NiftiHeader, args...; kwargs...)
    if isfile(hdr.aux_file)
        # TODO what format should it import as?
        load(aux_file, args...; kwargs...)
    else
        error("$aux_file not found")
    end
end

struct NeuroNameImage <: AbstractLabelImage
    labels::Dict{Integer,String}
end

function NeuroNameImage(ontology::String)
    # TODO this is much too large of a file to load into memory for each
    # image that may use it. Maybe provide function in `properties` which
    # appropriately indexes it, or extrat subset based on whats found in
    # the volume + user args. Also should be able to specify which species
    # atlas to use
    labels = load("~/.julia/packages/*/NeuroNames.xsl")
    NueroNameImage(ontology, labels)
end

function intent_labels(img::ImageMeta, hdr::NiftiHeader, args...; kwargs...)
    intent = get(NIFTI_INTENT, hdr.intent_code, false)
    img.properties["headaer"]["labels"] = intent(img, hdr)
    img
end

function setintent_labels!(hdr::NiftiHeader, img::ImageMeta;
                           aux_file::Union{String,Integer}=0)
    hdr.intent_code = NIFTI_INTENT_REVERSE[typeof(img.properties["header"]["intent"])]
    hdr.aux_file = aux_file
    if aux_file != 0
        # TODO need to figure out how it should be imported
        write(aux_file, img.properties["header"]["intent"].labels)
    end
end
