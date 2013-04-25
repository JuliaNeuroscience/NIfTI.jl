# dicom2nifti.jl
# Simple DICOM to NIfTI converter

# This script walks through the current working directory looking for DICOM
# files and sorts them into volumes. It then creates directories corresponding
# to each protocol, and saves nii files named by series number. It has
# (barely) been tested with DICOM files from a Siemens scanner. It may not
# work with other scanners. At present, it doesn't save orientation or slice
# information in the resulting NIfTI volumes, so DO NOT USE IT FOR ANYTHING
# BUT TESTING.

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

using DICOM, NIfTI

function walk(fn::Function, path::String)
    for fname in readdir(path)
    	if fname == "." || fname == ".."
    		continue
    	end
    	fpath = joinpath(path, fname)
    	if isdir(lstat(fpath))
    		walk(fn, fpath)
    	else
    		fn(fpath)
    	end
    end
end

function findtag(d, tag, from, what)
	idx = findfirst((x) -> x.tag == tag, d)
	if idx == 0
		error("$what missing from $from")
	end
	d[idx].data
end

function run()
	# Load all DICOMs into an array, indexed by series number
	dicoms = Dict{Int, Any}()
	walk(".") do fpath
		local d
		f = open(fpath, "r")
		try
			d = DICOM.dcm_parse(f)
		catch
			return
		finally
			close(f)
		end

		series_number_index = findfirst((x) -> x.tag == (0x0020, 0x0011), d)
		if series_number_index == 0
			println("No series number found for $fpath; ignoring")
			return
		end
		series_number = d[series_number_index].data[1]

		dicom_arr = get(dicoms, series_number, nothing)
		if dicom_arr == nothing
			dicoms[series_number] = dicom_arr = {}
		end
		push!(dicom_arr, d)
	end

	# Sort DICOMs
	for (series_number, slices) in dicoms
		d = slices[1]

		protocol_name = findtag(d, (0x0018, 0x1030), series_number, "Protocol Name")[1]
		image_position = findtag(d, (0x0020, 0x0032), series_number,
			"Image Position (Patient)")
		image_orientation = findtag(d, (0x0020, 0x0037), series_number,
			"Image Orientation (Patient)")

		# Sort slices
		# Formulae taken from http://nipy.sourceforge.net/nibabel/dicom/dicom_orientation.html
		z_dir_cos = cross(image_orientation[1:3], image_orientation[4:6])
		z_coord = sortby!(slices,
			(d) -> dot(findtag(d, (0x0020, 0x0032), series_number, "Image Position (Patient)"),
				z_dir_cos))

		raw = cat(3,
			[findtag(d, (0x7FE0, 0x0010), series_number, "Pixel Data")[1] for d in slices]...)

		# Create a name for each protocol
		protocol_dir = replace(protocol_name, "/", "")
		if !isdir(protocol_dir)
			mkdir(protocol_dir)
		end

		# Write NIfTI volumes
		ni = NIfTIVolume(raw)
		niftiwrite(joinpath(protocol_dir, "$(series_number).nii"), ni)
	end
end

run()