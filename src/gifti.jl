abstract type GiftiImage end

struct GiftiTimeSeries <: GiftiImage end
struct GiftiNodeIndex <: GiftiImage end
struct GiftiRGBVector <: GiftiImage end
struct GiftiRGBAVector <: GiftiImage end
struct GiftiShape <: GiftiImage end

function intent_gifti(intent_code::AbstractString, img::ImageMeta)
    if intent_code == "TimeSeries"
        """
        To signify that the value at each location is from a time series.
        """
        img.properties["header"]["intent"] = GiftiTimeSeries()
    elseif intent_code == "NodeIndex"
        """
        To signify that the value at each location is a node index, from a complete
        surface dataset.
        """
        img.properties["header"]["intent"] = GiftiNodeIndex()
    elseif intent_code == "RGBVector"
        """
        To signify that the vector value at each location is a 4 valued RGBA
        vector, of whatever type.
        - dataset must have a 5th dimension
        - dim[1] = 5
        - dim[2] = number of nodes
        - dim[3]= dim[4] = dim[5] = 1
        - dim[6] = 4
        """
        if ndims(img) == 5 && dim[3] == 1 && dim[4] == 1 &&
            dim[5] == 1 && dim[6] == 3
            # colorview(RGB, vol)
            img.properties["header"]["intent"] = GiftiRGBVector()
        else
            error("RGBVector intent must have:\n
                   dim[1] == 5\n
                   dim[3] == 1\n
                   dim[4] == 1\n
                   dim[5] == 1\n
                   dim[6] == 3")
        end
    elseif intent_code == "RGBAVector"
        """
        To signify that the vector value at each location is a 4 valued RGBA
        vector, of whatever type.
        - dataset must have a 5th dimension
        - dim[1] = 5
        - dim[2] = number of nodes
        - dim[3] = dim[4] = dim[5] = 1
        - dim[6] = 4
        """
        if dim[1] == 5 && dim[3] == 1 && dim[4] == 1 && dim[5] == 1 && dim[6] == 4
            #colorview(RGBA, vol)
            img.properties["header"]["intent"] = GiftiRGBAVector()
        else
            error("RGBVector intent must have:\n
                   dim[1] == 5\n
                   dim[3] == 1\n
                   dim[4] == 1\n
                   dim[5] == 1\n
                   dim[6] == 4")
        end
    elseif intent_code == "Shape"
        """
        To signify that the value at each location is a shape value,
        such as the curvature.
        """
        img.properties["header"]["intent"] = GiftiShape()
    end
end

function setintent_gifti!(hdr::NiftiHeader, img::ImageMeta, intenttype::Type)
    if intenttype == GiftiTimeSeries
        hdr.intent_code = NIFTI_INTENT_REVERSE[GiftiTimeSeries]
        hdr.intent_name = ""
        hdr.intent_p1 = 0
        hdr.intent_p2 = 0
        hdr.intent_p3 = 0
    elseif intenttype == GiftiNodeIndex
        hdr.intent_code = NIFTI_INTENT_REVERSE[GiftiNodeIndex]
        hdr.intent_name = ""
        hdr.intent_p1 = 0
        hdr.intent_p2 = 0
        hdr.intent_p3 = 0
    elseif intenttype == GiftiRGBVector
        hdr.intent_code = NIFTI_INTENT_REVERSE[GiftiRGBVector]
        hdr.intent_name = ""
        hdr.intent_p1 = 0
        hdr.intent_p2 = 0
        hdr.intent_p3 = 0
    elseif intenttype == GiftiRGBAVector
        hdr.intent_code = NIFTI_INTENT_REVERSE[GiftiRGBAVector]
        hdr.intent_name = ""
        hdr.intent_p1 = 0
        hdr.intent_p2 = 0
        hdr.intent_p3 = 0
    elseif intenttype == GiftiShape
        hdr.intent_code = NIFTI_INTENT_REVERSE[intenttype]
        hdr.intent_name = ""
        hdr.intent_p1 = 0
        hdr.intent_p2 = 0
        hdr.intent_p3 = 0
    end
end
