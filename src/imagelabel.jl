
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
function intent_labels(aux_file::AbstractString, img::ImageMeta, args...; kwargs...)
    if aux_file == "NeuroName"
        # TODO this is much too large of a file to load into memory for each
        # image that may use it. Maybe provide function in `properties` which
        # appropriately indexes it, or extrat subset based on whats found in
        # the volume + user args
        img.properties["header"]["labels"] = load("~/.julia/packages/*/NeuroNames.xsl", args...; kwargs...)
    elseif isfile(aux_file)
        img.properties["header"]["labels"] = load(aux_file, args...; kwargs...)
    else
        error("$aux_file not found")
    end
    img
end

