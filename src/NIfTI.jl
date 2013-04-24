# NIfTI.jl
# Methods for reading NIfTI MRI files in Julia

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

module NIfTI

using StrPack
import Base.getindex, Base.size, Base.ndims, Base.length, Base.endof, Base.write
export NIfTIVolume, niftiread, niftiwrite, voxel_size, time_step, vox

@struct type NIfTI1Header
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

type NIfTI1Extension
	ecode::Int32
	edata::Vector{Uint8}
end

type NIfTIVolume{T <: Number}
	header::NIfTI1Header
	extensions::Vector{NIfTI1Extension}
	raw::Array{T}
end

const SIZEOF_HDR = int32(348)

const NIfTI_DT_BITSTYPES = (Int16=>Type)[
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
const NIfTI_DT_BITSTYPES_REVERSE = (Type=>Int16)[]
for (k, v) in NIfTI_DT_BITSTYPES
	NIfTI_DT_BITSTYPES_REVERSE[v] = k
end

# Conversion factors to mm/ms
# http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/xyzt_units.html
const SPATIAL_UNIT_MULTIPLIERS = [
	1000,	# 1 => NIfTI_UNITS_METER
	1,		# 2 => NIfTI_UNITS_MM
	0.001	# 3 => NIfTI_UNITS_MICRON
]
const TIME_UNIT_MULTIPLIERS = [
	1000,	# NIfTI_UNITS_SEC
	1,		# NIfTI_UNITS_MSEC
	0.001,	# NIfTI_UNITS_USEC
	1,		# NIfTI_UNITS_HZ
	1,		# NIfTI_UNITS_PPM
	1		# NIfTI_UNITS_RADS
]

# Always in mm
voxel_size(header::NIfTI1Header) =
	header.pixdim[2:min(header.dim[1], 3)+1] * SPATIAL_UNIT_MULTIPLIERS[header.xyzt_units & int8(3)]

# Always in ms
time_step(header::NIfTI1Header) =
	header.pixdim[5] * TIME_UNIT_MULTIPLIERS[header.xyzt_units >> int8(3)]

function to_dim_info(dim_info::(Integer, Integer, Integer))
	if dim_info[1] > 3 || dim_info[1] < 0
		error("Invalid frequency dimension $(dim_info[1])")
	elseif dim_info[2] > 3 || dim_info[2] < 0
		error("Invalid phase dimension $(dim_info[2])")
	elseif dim_info[3] > 3 || dim_info[3] < 0
		error("Invalid slice dimension $(dim_info[3])")
	end

	int8(dim_info[1] | (dim_info[2] << int8(2)) | (dim_info[3] << int8(4)))
end

# Returns or sets dim_info as a tuple whose values are the frequency, phase, and slice dimensions
dim_info(header::NIfTI1Header) = (header.dim_info & int8(3), (header.dim_info >> 2) & int8(3),
	(header.dim_info >> 4) & int8(3))
dim_info{T <: Integer}(header::NIfTI1Header, dim_info::(T, T, T)) =
	header.dim_info = to_dim_info(dim_info)

# Gets dim to be used in header
function niftidim(x::Array)
	dim = zeros(Int16, 8)
	dim[1] = ndims(x)
	dim[2:dim[1]+1] = [size(x)...]
	dim
end

# Gets datatype to be used in header
function niftidatatype(t::Type)
	t = get(NIfTI_DT_BITSTYPES_REVERSE, t, nothing)
	if t == nothing
		error("Unsupported data type $T")
	end
	t
end

# Gets the size of a type in bits
niftibitpix(t::Type) = int16(sizeof(t)*8)

# Validates srow values
function validate_srow(srow_x::Vector{Float32}, srow_y::Vector{Float32}, srow_z::Vector{Float32})
	if length(srow_x) != 4
		error("srow_x must have 4 values; $(length(srow_x)) encountered")
	elseif length(srow_y) != 4
		error("srow_y must have 4 values; $(length(srow_y)) encountered")
	elseif length(srow_z) != 4
		error("srow_z must have 4 values; $(length(srow_z)) encountered")
	end
end

# Constructor
function NIfTIVolume{T <: Union(Number, Nothing)}(
# Optional MRI volume; if not given, an empty volume is used
raw::T=nothing,
extensions::Union(Vector{NIfTI1Extension}, Nothing)=nothing;

# Fields specified as UNUSED in NIfTI1 spec
data_type::String="", db_name::String="", extents::Integer=int32(0),
session_error::Integer=int16(0), regular::Integer=int8(0), glmax::Integer=int32(0),
glmin::Integer=int16(0),

# The frequency encoding, phase encoding and slice dimensions.
dim_info::NTuple{3, Integer}=(0, 0, 0),
# Describes data contained in the file; for valid values, see
# http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/group__NIfTI1__INTENT__CODES.html
intent_p1::Real=0f0, intent_p2::Real=0f0, intent_p3::Real=0f0,
intent_code::Integer=int16(0), intent_name::String="",
# Information about which slices were acquired, in case the volume has been padded
slice_start::Integer=int16(0), slice_end::Integer=int16(0), slice_code::Int8=int8(0),
# The size of each voxel and the time step. These are formulated in mm unless xyzt_units is
# explicitly specified
voxel_size::NTuple{3, Real}=(0f0, 0f0, 0f0), time_step::Real=0f0, xyzt_units::Int8=int8(18),
# Slope and intercept by which volume shoudl be scaled
scl_slope::Real=0f0, scl_inter::Real=0f0,
# These describe how data should be scaled when displayed on the screen. They are probably
# rarely used
cal_max::Real=0f0, cal_min::Real=0f0,
# The amount of time taken to acquire a slice
slice_duration::Real=0f0,
# Indicates a non-zero start point for time axis
toffset::Real=0f0,

# "Any text you like"
descrip::String="",
# Name of auxiliary file
aux_file::String="",

# Transform of NIfTI file. See NIfTI spec
qform_code::Integer=int16(0), sform_code::Integer=int16(0), quatern_b::Real=0f0,
quatern_c::Real=0f0, quatern_d::Real=0f0, qoffset_x::Real=0f0, qoffset_y::Real=0f0,
qoffset_z::Real=0f0, srow_x::Vector{Float32}=Float32[0, 0, 0, 0],
srow_y::Vector{Float32}=Float32[0, 0, 0, 0], srow_z::Vector{Float32}=Float32[0, 0, 0, 0])
	local t
	if raw == nothing
		raw = Int16[]
		t = Int16
	else
		t = T
	end

	if extensions == nothing
		extensions = NIfTI1Extension[]
	end

	validate_srow(srow_x, srow_y, srow_z)

	NIfTIVolume(NIfTI1Header(SIZEOF_HDR, data_type, db_name, int32(extents), int16(session_error),
		int8(regular), to_dim_info(dim_info), niftidim(raw), float32(intent_p1), float32(intent_p2),
		float32(intent_p3), int16(intent_code), niftidatatype(t), niftibitpix(t), 
		int16(slice_start), float32([voxel_size..., time_step, 0, 0, 0, 0]), float32(352),
		float32(scl_slope), float32(scl_inter), int16(slice_end), int8(slice_code),
		int8(xyzt_units), float32(cal_max), float32(cal_min), float32(slice_duration),
		float32(toffset), int32(glmax), int32(glmin), descrip, aux_file, int16(qform_code),
		int16(sform_code), float32(quatern_b), float32(quatern_c), float32(quatern_d),
		float32(qoffset_x), float32(qoffset_y), float32(qoffset_z), srow_x, srow_y, srow_z,
		intent_name, "n+1"), extensions, raw)
end

# Calculates the size of a NIfTI extension
esize(ex::NIfTI1Extension) = 8 + iceil(length(ex.edata)/16)*16

# Validates the header of a volume and updates it to match the volume's contents
function niftiupdate{T}(vol::NIfTIVolume{T})
	validate_srow(vol.header.srow_x, vol.header.srow_y, vol.header.srow_z)
	vol.header.sizeof_hdr = SIZEOF_HDR
	vol.header.dim = niftidim(vol.raw)
	vol.header.datatype = niftidatatype(T)
	vol.header.bitpix = niftibitpix(T)
	vol.header.vox_offset = isempty(vol.extensions) ? int32(352) :
		float32(reduce((v, ex) -> v + esize(ex), SIZEOF_HDR, vol.extensions))
	vol
end

# Write the NIfTI header
writeheader(io::IO, vol::NIfTIVolume) = pack(io, niftiupdate(vol).header)

# Write a NIfTI file
function write(io::IO, vol::NIfTIVolume)
	writeheader(io, vol)
	if isempty(vol.extensions)
		write(io, int32(0))
	else
		for ex in vol.extensions
			sz = esize(ex)
			write(io, int32(sz))
			write(io, int32(ex.ecode))
			write(io, ex.edata)
			write(io, zeros(Uint8, sz - length(ex.edata)))
		end
	end
	write(io, vol.raw)
end

# Convenience function to write a NIfTI file given a patch
function niftiwrite(path::String, vol::NIfTIVolume)
	io = open(path, "w")
	write(io, vol)
	close(io)
end

# Read header from a NIfTI file
function read_header(io::IO)
	header = unpack(io, NIfTI1Header)
	if header.sizeof_hdr != SIZEOF_HDR
		error("This is not a NIfTI-1 file")
	end
	if header.dim[1] > 7
		error("Byte swapping not yet supported")
	end
	header
end

# Read extension fields following NIfTI header
function read_extensions(io::IO, header::NIfTI1Header)
	if eof(io)
		return NIfTI1Extension[]
	end

	extension = read(io, Uint8, 4)
	if extension[1] != 1
		return NIfTI1Extension[]
	end

	extensions = NIfTI1Extension[]
	should_bswap = header.dim[1] > 7
	while header.magic == "n+1" ? position(io) < header.vox_offset : !eof(io)
		esize = read(io, Int32)
		ecode = read(io, Int32)
		
		if should_bswap
			esize = bswap(esize)
			ecode = bswap(ecode)
		end

		push!(extensions, NIfTI1Extension(ecode, read(io, Uint8, esize-8)))
	end
	extensions
end

# Read a NIfTI file. The optional mmap argument determines whether the contents are read in full
# (if false) or mmaped from the disk (if true).
function niftiread(file::String; mmap::Bool=false)
	header_io = open(file, "r")
	header = read_header(header_io)
	extensions = read_extensions(header_io, header)
	dims = tuple(int(header.dim[2:header.dim[1]+1])...)

	if !has(NIfTI_DT_BITSTYPES, header.datatype)
		error("Data type $(header.datatype) not yet supported")
	end
	dtype = NIfTI_DT_BITSTYPES[header.datatype]

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

	return NIfTIVolume(header, extensions, volume)
end

# Allow file to be indexed like an array, but with indices yielding scaled data
getindex(f::NIfTIVolume, args...) = getindex(f.raw, args...) * f.header.scl_slope + f.header.scl_inter
vox(f::NIfTIVolume, args...) =
	getindex(f, [isa(args[i], Colon) ? (1:size(f.raw, i)) : args[i] + 1 for i = 1:length(args)]...)
size(f::NIfTIVolume) = size(f.raw)
size(f::NIfTIVolume, d) = size(f.raw, d)
ndims(f::NIfTIVolume) = ndims(f.raw)
length(f::NIfTIVolume) = length(f.raw)
endof(f::NIfTIVolume) = endof(f.raw)

end