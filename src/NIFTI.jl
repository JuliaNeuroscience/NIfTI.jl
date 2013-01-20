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

import Base.read
export niftiread

function struct_string(fn::Function, bytes::Vector{Uint8})
    last_byte = findfirst(bytes, 0)
    fn(last_byte == 0 ? bytes : bytes[1:last_byte-1])
end

macro struct(typename, contents)
	if contents.head != :block
		error("Invalid struct declaration")
	end

	blk = expr(:call, typename)

	for typedecl in contents.args
		if isa(typedecl, LineNumberNode)
			continue
		elseif typedecl.head != :(::)
			error("Invalid struct declaration")
		end
		fieldname = typedecl.args[1]
		fieldtype = typedecl.args[2]
		if isa(fieldtype, Expr)
			if fieldtype.head != :call
				error("Invalid struct declaration")
			end
			# Type has size parameters

			if contains((:ASCIIString, :UTF8String, :String), fieldtype.args[1])
				typedecl.args[2] = fieldtype.args[1]
				push!(blk.args, :(struct_string($(fieldtype.args[1] == :ASCIIString ? :ascii : :utf8), read(ios, Uint8, ($(fieldtype.args[2:end]...))))))
			else
				typedecl.args[2] = :(Array{$(fieldtype.args[1]), $(length(fieldtype.args)-1)})
				push!(blk.args, :(read(ios, $(fieldtype.args[1]), ($(transpose(fieldtype.args[2:end])...)))))
			end
		else
			push!(blk.args, :(read(ios, $fieldtype)))
		end
	end
	
	quote
		global read
		type $typename
			$contents
		end
		read(ios::IOStream, ::Type{$typename}) = $blk
	end
end

@struct NIFTI1Header begin
	sizeof_hdr::Int32
	data_type::ASCIIString(10)
	db_name::ASCIIString(18)
	extents::Int32
	session_error::Int16
	regular::Int8
	dim_info::Int8

	dim::Int16(8)
	intent_p1::Float32
	intent_p2::Float32
	intent_p3::Float32
	intent_code::Int16
	datatype::Int16
	bitpix::Int16
	slice_start::Int16
	pixdim::Float32(8)
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

	srow_x::Float32(4)
	srow_y::Float32(4)
	srow_z::Float32(4)

	intent_name::ASCIIString(16)

	magic::ASCIIString(4)
end

type MRIVolume{T}
	raw::Array{T}
end

type NIFTI1Extension
	code::Int32
	edata::Vector{Uint8}
end

type NIFTIFile
	header::NIFTI1Header
	extensions::Vector{NIFTI1Extension}
	volume::MRIVolume
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

function niftiread_header(io::IO)
	header = read(io, NIFTI1Header)
	if header.sizeof_hdr != 348
		error("This is not a NIFTI-1 file")
	end
	if header.dim[1] > 7
		error("Byte swapping not yet supported")
	end
	header
end

function nifitread_extensions(io::IO, header::NIFTI1Header)
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

function niftiread(file::String, mmap::Bool)
	header_io = open(file, "r")
	header = niftiread_header(header_io)
	extensions = nifitread_extensions(header_io, header)
	dims = tuple(int(reverse(header.dim[2:header.dim[1]+1]))...)

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

	return NIFTIFile(header, extensions, MRIVolume(volume))
end
niftiread(file::String) = niftiread(file, false)

end