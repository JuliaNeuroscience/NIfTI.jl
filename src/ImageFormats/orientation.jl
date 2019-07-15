# TODO would be neat to have some sort of function that transforms between these views
# (e.g. "right-anterior-superior" => (:R,:A,:S) => (:R2L, :A2P, :S2I))

const space2axes = Dict(
    "right-anterior-superior" => (:R,:A,:S),
    "ras" => (:R,:A,:S),
    "left-anterior-superior" => (:L,:A,:S),
    "las" => (:L,:A,:S),
    "left-posterior-superior" => (:L,:P,:S),
    "lps" => (:L,:P,:S),
    "right-anterior-superior-time" => (:R,:A,:S,:time),
    "rast" => (:R,:A,:S,:time),
    "left-anterior-superior-time" => (:L,:A,:S,:time),
    "last" => (:L,:A,:S,:time),
    "left-posterior-superior-time" => (:L,:P,:S,:time),
    "lpst" => (:L,:P,:S,:time),
    "scanner-xyz" => (:scannerx,:scannery,:scannerz),
    "scanner-xyz-time" => (:scannerx,:scannery,:scannerz,:time),
    "3d-right-handed" => (:xrcs,:yrcs,:zrcs),
    "3d-left-handed" => (:xlcs,:ylcs,:zlcs),
    "3d-right-handed-time" => (:xrcs,:yrcs,:zrcs,:time),
    "3d-left-handed-time" => (:xlcs,:ylcs,:zlcs,:time))

const axes2space = Dict(
    (:R,:A,:S) => "right-anterior-superior",
    (:L,:A,:S) => "left-anterior-superior",
    (:L,:P,:S) => "left-posterior-superior",
    (:R,:A,:S,:time) => "right-anterior-superior-time",
    (:L,:A,:S,:time) => "left-anterior-superior-time",
    (:L,:P,:S,:time) => "left-posterior-superior-time",
    (:scannerx,:scannery,:scannerz) => "scanner-xyz",
    (:scannerx,:scannery,:scannerz,:time) => "scanner-xyz-time",
    (:xrcs,:yrcs,:zrcs) => "3d-right-handed",
    (:xlcs,:ylcs,:zlcs) => "3d-left-handed",
    (:xrcs,:yrcs,:zrcs,:time) => "3d-right-handed-time",
    (:xlcs,:ylcs,:zlcs,:time) => "3d-left-handed-time")

const num2axes = Dict(
    1  => :L2R,
    -1 => :R2L,
    2  => :P2A,
    -2 => :A2P,
    3  => :I2S,
    -3 => :S2I)

const axes2num = Dict(
    :L2R => 1,
    :R2L => -1,
    :P2A => 2,
    :A2P => -2,
    :I2S => 3,
    :S2I => -3)

function isradview(axnames::NTuple{3,Symbol})
    axnames == (:L, :A, :S) | axnames == (:L2R, :A2P, :S2I)
end

function isneuroview(axnames::NTuple{3,Symbol})
    axnames == (:R, :A, :S) | axnames == (:R2L, :A2P, :S2I)
end

function isradview(axnames::String)
    axnames == "left-anterior-superior" | axnames == "left-anterior-superior-time"
end

function isneuroview(axnames::String)
    axnames == "right-anterior-superior" | axnames == "right-anterior-superior-time"
end

getaffine(img::Union{AbstractArray,ImageStream,ImageInfo}) =
    AffineMap(getlinear(img), gettranslation(img))

# linear component of AffineMap
getlinear(img::Union{AbstractArray,ImageStream,ImageInfo}) =
    _getlinear(spacedirections(img))

function _getlinear(sd::NTuple{N,NTuple{N,T}}) where {N,T<:Integer}
    _getlinear(Tuple([float.(x) for x in sd]))
end

function _getlinear(sd::NTuple{2,NTuple{2,T}}) where T<:AbstractFloat
    SMatrix{2,2,T}(sd[1][1], sd[2][1],
                   sd[1][2], sd[2][2])
end

function _getlinear(sd::NTuple{3,NTuple{3,T}}) where T<:AbstractFloat
    SMatrix{3,3,T}(sd[1][1], sd[2][1], sd[3][1],
                   sd[1][2], sd[2][2], sd[3][2],
                   sd[1][3], sd[2][3], sd[3][3])
end

# translation component of AffineMap
gettranslation(img::Union{AbstractArray,ImageStream,ImageInfo}) =
    SVector(ustrip.(first.(axisvalues(img)[1:sdims(img)])))

# quaternion of spacedirections
getquat(img::Union{AbstractArray,ImageStream,ImageInfo}) = _getquat(getlinear(img))

function _getquat(sd::SMatrix{2,2,T}) where T<:AbstractFloat
    _getquat(SMatrix{3,3,T}(sd[1:2, 1]..., zero(T),
                            sd[1:2, 2]..., zero(T),
                            zero(T),  zero(T),  one(T)))
end

_getquat(sd::StaticMatrix{3,3,T}) where T<:AbstractFloat = Quat(sd)

_getquat(sd::NTuple{N,NTuple{N,T}}) where {N,T} = _getquat(map(i->float.(i), sd))

getaffinemat(img::Union{AbstractArray,ImageStream,ImageInfo}) =
    _getaffinemat(getlinear(img), gettranslation(img))


function _getaffinemat(lin::SMatrix{N,N,Tlin}, trans::SVector{N,Ttrans}) where {N,Tlin,Ttrans}
    T = promote_type(Tlin, Ttrans)
    _getaffinemat(T.(lin), T.(trans))
end

function _getaffinemat(lin::SMatrix{2,2,T}, trans::SVector{2,T}) where T
    SMatrix{3,3,T}(lin[1:2, 1]..., zero(T),
                   lin[1:2, 2]..., zero(T),
                   trans..., sign(det(lin)))
end

function _getaffinemat(lin::SMatrix{3,3,T}, trans::SVector{3,T}) where T
    SMatrix{4,4,T}(lin[1:3, 1]..., zero(T),
                   lin[1:3, 2]..., zero(T),
                   lin[1:3, 3]..., zero(T),
                   trans..., sign(det(lin)))
end
