# NIFTI.jl
# Methods for reading NIFTI MRI files in Julia

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

module NIFTI

using StrPack
import Base.getindex, Base.size, Base.ndims, Base.length, Base.endof
export niftiread, voxelsize, tr, vox

macro a(xpr...)
	println(xpr)
end

@struct type NIFTI1Header
	sizeof_hdr::Int32
	data_type::ASCIIString(10)
	db_name::ASCIIString(18)
	extents::Int32
	session_error::Int16
	regular::Int8
	dim_info::Int8

	dim::Vector{Int16}(8)
	intent_p1::Float32
	intent_p2::Float32
	intent_p3::Float32
	intent_code::Int16
	datatype::Int16
	bitpix::Int16
	slice_start::Int16
	pixdim::Vector{Float32}(8)
	vox_offset::Float32
	scl_slope::Float32
	scl_inter::Float32
	slice_end::Int16
	slice_code::Int8
	xyzt_units::Int8
	cal_max::Float32
	cal_min::Float32
	slice_duration::Float32
	toffset::Float32
	glmax::Int32
	glmin::Int32

	descrip::ASCIIString(80)
	aux_file::ASCIIString(24)

	qform_code::Int16
	sform_code::Int16
	quatern_b::Float32
	quatern_c::Float32
	quatern_d::Float32
	qoffset_x::Float32
	qoffset_y::Float32
	qoffset_z::Float32

	srow_x::Vector{Float32}(4)
	srow_y::Vector{Float32}(4)
	srow_z::Vector{Float32}(4)

	intent_name::ASCIIString(16)

	magic::ASCIIString(4)
end align_packed

type NIFTI1Extension
	code::Int32
	edata::Vector{Uint8}
end

type NIFTIFile{T}
	header::NIFTI1Header
	extensions::Vector{NIFTI1Extension}
	raw::Array{T}
end

const NIFTI_DT_BITSTYPES = (Int8=>Type)[
	2 => Uint8,
	4 => Int16,
	8 => Int32,
	16 => Float32,
	32 => Complex64,
	64 => Float64,
	256 => Int8,
	512 => Uint16, 
	768 => Uint32,
	1024 => Int64,
	1280 => Uint64,
	1792 => Complex128
]

# Conversion factors to mm/ms
# http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/xyzt_units.html
const SPATIAL_UNIT_MULTIPLIERS = [
	1000,	# 1 => NIFTI_UNITS_METER
	1,		# 2 => NIFTI_UNITS_MM
	0.001	# 3 => NIFTI_UNITS_MICRON
]
const TIME_UNIT_MULTIPLIERS = [
	1000,	# NIFTI_UNITS_SEC
	1,		# NIFTI_UNITS_MSEC
	0.001,	# NIFTI_UNITS_USEC
	1,		# NIFTI_UNITS_HZ
	1,		# NIFTI_UNITS_PPM
	1		# NIFTI_UNITS_RADS
]

# Read header from a NIFTI file
function read_header(io::IO)
	header = unpack(io, NIFTI1Header)
	if header.sizeof_hdr != 348
		error("This is not a NIFTI-1 file")
	end
	if header.dim[1] > 7
		error("Byte swapping not yet supported")
	end
	header
end

# Read extension fields following NIFTI header
function read_extensions(io::IO, header::NIFTI1Header)
	if eof(io)
		return NIFTI1Extension[]
	end

	extension = read(io, Uint8, 4)
	if extension[1] != 1
		return NIFTI1Extension[]
	end

	extensions = NIFTI1Extension[]
	should_bswap = header.dim[1] > 7
	while header.magic == "n+1" ? position(io) < header.vox_offset : !eof(io)
		esize = read(io, Int32)
		code = read(io, Int32)
		
		if should_bswap
			esize = bswap(esize)
			code = bswap(code)
		end

		push!(extensions, NIFTI1Extension(code, read(io, Uint8, esize-8)))
	end
	extensions
end

# Read a NIFTI file. The optional mmap argument determines whether the contents are read in full
# (if false) or mmaped from the disk (if true).
function niftiread(file::String; mmap::Bool=false)
	header_io = open(file, "r")
	header = read_header(header_io)
	extensions = read_extensions(header_io, header)
	dims = tuple(int(header.dim[2:header.dim[1]+1])...)

	if !has(NIFTI_DT_BITSTYPES, header.datatype)
		error("Data type $(header.datatype) not yet supported")
	end
	dtype = NIFTI_DT_BITSTYPES[header.datatype]

	local volume
	if header.magic == "n+1"
		if mmap
			volume = mmap_array(dtype, dims, header_io, int(header.vox_offset))
		else
			seek(header_io, int(header.vox_offset))
			volume = read(header_io, dtype, dims)
			@assert eof(header_io)
			close(header_io)
		end
	elseif header.magic == "nil"
		close(header_io)
		volume_io = open(replace(file, r"\.\w+$", "")*".img", "r")
		if mmap
			volume = mmap_array(dtype, dims, volume_io)
		else
			volume = read(volume_io, dtype, dims)
			close(volume_io)
		end
	end

	return NIFTIFile(header, extensions, volume)
end

# Always in mm
voxelsize(header::NIFTI1Header) =
	header.pixdim[2:4] * SPATIAL_UNIT_MULTIPLIERS[header.xyzt_units & 7]

# Always in ms
tr(header::NIFTI1Header) = header.pixdim[5] * TIME_UNIT_MULTIPLIERS[header.xyzt_units >> 3]

# Allow file to be indexed like an array, but with indices yielding scaled data
getindex(f::NIFTIFile, args...) = getindex(f.raw, args...) * f.header.scl_slope + f.header.scl_inter
vox(f::NIFTIFile, args...) =
	getindex(f, [isa(args[i], Colon) ? (1:size(f.raw, i)) : args[i] + 1 for i = 1:length(args)]...)
size(f::NIFTIFile) = size(f.raw)
size(f::NIFTIFile, d) = size(f.raw, d)
ndims(f::NIFTIFile) = ndims(f.raw)
length(f::NIFTIFile) = length(f.raw)
endof(f::NIFTIFile) = endof(f.raw)

end