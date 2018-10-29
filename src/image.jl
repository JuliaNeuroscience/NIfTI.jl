# functions for organizing NIfTI -> Julia and Julia -> NIfTI

"""
"""
function nii2img(vol::Array,
                 ext::NiftiExtension,
                 hdr::NiftiHeader)
    header = Dict{String,Any}()

    qfac = (hdr.pixdim[1] < 0.0) ? -1.0 : 1.0
    sto_xyz, qto_xyz = getaffine(hdr.qform_code, hdr.sform_code,
                                 hdr.srow_x, hdr.srow_y, hdr.srow_x,
                                 hdr.quatern_b, hdr.quatern_c,
                                 hdr.quatern_d, hdr.quatern_x,
                                 hdr.quatern_y, hdr.quatern_z,
                                 hdr.pixdim[2], hdr.pixdim[3],
                                 hdr.pixdim[4], qfac)
    # TODO
    # - ensure orientation works correctly
    # - check sto_xyz vs qto_xyz implementation
    #   - both have value (sometimes)
    ori = mat2ori(sto_xyz)

    intent_code = NIFTI_INTENT[hdr.intent_code]
    # create AxisArray
    # TODO consolidate this step
    if intent_code == "Dimensionless"
        img = vol
    else
        su = NIFTI_UNITS[hdr.xyzt_units & 0x07]
        axs = (Axis{ori[1]}(hdr.pixdim[1]*su))
        if dim[1] > 1
            axs = (axs..., Axis{ori[2]}(hdr.pixdim[2]*su))
        end
        if dim[1] > 2
            axs = (axs..., Axis{ori[3]}(hdr.pixdim[3]*su))
        end
        if dim[1] > 3
            # time unit
            axs = (axs..., Axis{:time}(timing(volinfo)))
        end
        img = ImageMeta(AxisArray(vol, axs...),
                        properties = Dict{String,Any}(
                                         "filetype" => ftype,
                                         "header" => Dict{}(
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
                                                     header["extension"] => ext)))
    end

    if hdr.intent_code > 1 && hdr.intent_code < 25
        intent_stats(intent_code, img)
    elseif intent_code == "Label"
        intent_label(hdr.aux_file, img, args...; kwargs...)
    elseif intent_code == "NeuroName"
        intent_label("NeuroName", img, args...; kwargs...)
    elseif intent_code > 1003 && intent_code < 1010
        # TODO VectorImages
        intent_vector(intent_code, img)
    elseif intent_code == "Quaternion"
        # TODO
        error("$intent_code is not yet supported")
    elseif intent_code > 2000 && intent_code < 2006
        # TODO
        intent_gifti(intent_code, img)
    end
end

"""
    img2nii(img::ImageMeta, version::Int) -> NiftiHeader, NiftiExtension, vol::AbstractArray
"""
function img2nii(img::ImageMeta, version::Int)
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
                          img.properties["header"]["intent_name"], nothing)
    hdr.intent_p1 = img.properties["header"]["intent_p1"][1]
    hdr.intent_p2 = img.properties["header"]["intent_p2"][1]
    hdr.intent_p3 = img.properties["header"]["intent_p3"][1]

    # set slice values
    setslicing!(hdr, img)

    # set orientation values
    setaffine!(hdr, img)
    setoffset!(hdr, img)

    hdr.scl_slope = 0
    hdr.scl_inter = 0
    hdr.cal_max = img.properties["header"]["cal_max"]
    hdr.cal_min = img.properties["header"]["cal_min"]
    hdr.aux_file = img.properties["header"]["aux_file"]
    hdr.descrip = img.properties["descrip"]
    hdr
end
