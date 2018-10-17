mutable struct NIfTI1Extension
    ecode::Int32
    edata::Vector{UInt8}
end

# Calculates the size of a NIfTI extension
esize(ex::NIfTI1Extension) = 8 + ceil(Int, length(ex.edata)/16)*16


