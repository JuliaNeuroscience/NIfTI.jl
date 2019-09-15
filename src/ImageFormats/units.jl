import Unitful
using Unitful: @unit, @u_str

# MR units
@unit Gauss "Gauss" Gauss 1e-4u"T" false # Gs occupied by Giga Sec in Unitful
 
# http://painterqubits.github.io/Unitful.jl/stable/extending/#Precompilation
const localunits = Unitful.basefactors
function __init__()
    merge!(Unitful.basefactors, localunits)
    Unitful.register(UnitfulMR)
end