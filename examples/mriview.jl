# mriview.jl
# Simple MRI viewer

# This script displays a NIfTI volume using ImageView.

# Copyright (C) 2013   Simon Kornblith

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
