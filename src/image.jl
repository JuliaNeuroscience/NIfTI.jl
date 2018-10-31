# functions for organizing NIfTI -> Julia and Julia -> NIfTI

"""
"""
function nii2img(vol::AbstractArray,
                 ext::Vector{NiftiExtension},
                 hdr::NiftiHeader)
    header = Dict{String,Any}()

    sto_xyz, qto_xyz = getaffine(hdr)
    # TODO keep both?
    to_xyz = hdr.sform_code > 0 ? sto_xyz : qto_xyz
    xform_code = hdr.sform_code > 0 ? hdr.sform_code : hdr.qform_code
    # TODO
    # - ensure orientation works correctly
    # - check sto_xyz vs qto_xyz implementation
    #   - both have value (sometimes)
    ori = mat2ori(sto_xyz)

    intent_code = get(NIFTI_INTENT, hdr.intent_code, 0)
    # create AxisArray
    # TODO consolidate this step
    if intent_code == DimensionlessImage
        axsimg = vol
    else
        su = NIFTI_UNITS[hdr.xyzt_units & 0x07]
        axs = (Axis{ori[1]}(range(1, step=hdr.pixdim[2], length=hdr.dim[2])*su))
        if hdr.dim[1] > 1
            axs = (axs..., Axis{ori[2]}(range(1, step=hdr.pixdim[3], length=hdr.dim[3])*su))
        end
        if hdr.dim[1] > 2
            axs = (axs..., Axis{ori[3]}(range(1, step=hdr.pixdim[4], length=hdr.dim[4])*su))
        end
        if hdr.dim[1] > 3
            # time unit
            axs = (axs..., Axis{:time}(timedim(hdr)))
        end
        #axsimg = AxisArray(vol,
        #                   ([or[i] for i in 1:nd]...,),  # axisnames
        #                   ([]...,),
        #                   ([]...,))
        axsimg = AxisArray(vol, axs...,)
    end
    img = ImageMeta(axsimg,
                    properties = Dict{String,Any}(
                         "filetype" => "NIfTI",
                         "spacedirections" => to_xyz,
                         "descrip" => getdescription(hdr),
                         "header" => Dict{String,Any}(
                                     # can extract from data{T}
                                     # header["bitpix"] => hdr.bitpix,
                                     # Stuff that doesn't fit anywhere else
                                     header["slice_code"] => hdr.slice_code,
                                     header["slice_start"] => hdrslice_start,
                                     header["slice_end"] => hdr.slice_end,
                                     header["slice_duration"] => hdr.slice_duration,
                                     header["cal_max"] => hdr.cal_max,
                                     header["cal_min"] => hdr.cal_min,
                                     header["dim_info"] => diminfo(hdr.dim_info),
                                     header["extension"] => ext,
                                     header["xform"] => get(NIFTI_XFORM, xform_code, nothing))))
    intent_nifti(img, hdr)
end

"""
    img2nii(img::ImageMeta, version::Int) -> NiftiHeader, NiftiExtension, vol::AbstractArray
"""
function img2nii(img::ImageMeta, version::Int; description::String="")
    if version == 1
        hdr = Nifti1Header()
    elseif version == 2
        hdr = Nifti2Header()
    end

    # set units/dims/stepsize
    setdim!(hdr, img)
    setunits!(hdr, img)
    setpixdim!(hdr, img)

    # get base element type of img
    hdr.datatype = getdatatype(T)
    hdr.bitpix = nibitpix(T)

    # set intent values
    hdr.intent_name = img.properties["header"]["intent_name"]
    hdr.intent_code = get(NIFTI_INTENT_REVERSE,
                          img.properties["header"]["intent_name"],
                          nothing)
    hdr.intent_p1 = img.properties["header"]["intent_p1"][1]
    hdr.intent_p2 = img.properties["header"]["intent_p2"][1]
    hdr.intent_p3 = img.properties["header"]["intent_p3"][1]

    # set slice values
    setslicing!(hdr, img)

    # set orientation values
    setaffine!(hdr, img)
    setoffset!(hdr, img)

    setdescription!(hdr, get(img.properties["header"], "descrip", description))

    hdr.scl_slope = 0
    hdr.scl_inter = 0
    hdr.cal_max = img.properties["header"]["cal_max"]
    hdr.cal_min = img.properties["header"]["cal_min"]
    hdr.aux_file = img.properties["header"]["aux_file"]
    hdr.descrip = img.properties["descrip"]
    return hdr, getvolume(img, version), getextension(img)
end
