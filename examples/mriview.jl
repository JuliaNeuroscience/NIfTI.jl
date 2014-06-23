# mriview.jl
# Simple MRI viewer

# This script displays a NIfTI volume using ImageView.

using NIfTI, Tk, Images, ImageView

function mriview(ni::NIVolume, xy=["x", "z"], zoom::Int=5)
    img = Image(ni.raw, ["spatialorder" => ["x", "y", "z"], "timedim" => 4, "colorspace" => "Gray",
                         "pixelspacing"=>voxel_size(ni.header)])
    img = scale(scaleminmax(img), img)
    imgc, imgslice = display(img, xy=xy)
    tksize = get_size(toplevel(imgc))
    controlheight = tksize[2] - size(img, xy[2])
    set_size(toplevel(imgc), max(tksize[1], size(img, xy[1])*zoom),
                             size(img, xy[2])*zoom+controlheight)
    (imgc, imgslice)
end
mriview(path::String) = mriview(niread(path, mmap=true))
