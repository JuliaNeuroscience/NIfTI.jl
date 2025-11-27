module NIfTIMakieExt

import Makie
import NIfTI: NIVolume

Makie.convert_single_argument(v::NIVolume) = Array(v)

end