function intent_gifti(intent_code::AbstractString, img::ImageMeta)

    if intent_code == "TimeSeries"
        """
        To signify that the value at each location is from a time series.
        """
    elseif intent_code == "NodeIndex"
        """
        To signify that the value at each location is a node index, from a complete
        surface dataset.
        """

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
            RGBVector()
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
            RGBAVector()
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
    end
end
