# UnitfulPositions?
module SemanticPositions

export # types
       SemanticDimension,
       SemanticPosition,
       SemanticOrder,
       SemanticOrientation,
       # constants
       Sagittal,
       Axial,
       Coronal,
       Left,
       Right,
       Superior,
       Inferior,
       Anterior,
       Posterior,
       LeftToRight,
       RightToLeft,
       SuperiorToInferior,
       AnteriorToPosterior,
       PosteriorToAnterior,
       # methods
       isradiologic,
       isneurologic,
       num2axis,
       num2ori,
       ori2numb


"""
    SemanticDimension
"""
struct SemanticDimension{D} end

const Sagittal = SemanticDimension{:sagittal}()
const Axial = SemanticDimension{:axial}()
const Coronal = SemanticDimension{:coronal}()

"""
    SemanticPosition
"""
struct SemanticPosition{P,D} end

dimension(::S) where {S<:SemanticPosition} = dimension(S)
dimension(::Type{SemanticPosition{S,D}}) where {S,D} = D

const Left = SemanticPosition{:left,Sagittal}()
const Right = SemanticPosition{:right,Sagittal}()

const Superior = SemanticPosition{:superior,Axial}()
const Inferior = SemanticPosition{:inferior,Axial}()

const Anterior = SemanticPosition{:anterior,Coronal}()
const Posterior = SemanticPosition{:posterior,Coronal}()

eqdim(::SemanticPosition{P1,D}, ::SemanticPosition{P2,D}) where {P1,P2,D} = true
eqdim(::SemanticPosition{P1,D1}, ::SemanticPosition{P2,D2}) where {P1,P2,D1,D2} = false


"""
    SemanticOrder
"""
struct SemanticOrder{F,L}

    function SemanticOrder{F,L}() where {F,L}
        eqdim(F,L) || error("first and last semantic positions must come from same dimension.")
        new{F,L}()
    end
end

SemanticOrder(s::SemanticOrder) = s
SemanticOrder(f::SemanticPosition, l::SemanticPosition) = SemanticOrder{f,l}()

Base.first(::SemanticOrder{F,L}) where {F,L} = F
Base.last(::SemanticOrder{F,L}) where {F,L} = L

const LeftToRight = SemanticOrder(Left, Right)
const RightToLeft = SemanticOrder(Right, Left)

const SuperiorToInferior = SemanticOrder(Superior, Inferior)
const InferiorToSuperior = SemanticOrder(Inferior, Superior)

const AnteriorToPosterior = SemanticOrder(Anterior, Posterior)
const PosteriorToAnterior = SemanticOrder(Posterior, Anterior)

"""
    SemanticOrientation
"""
struct SemanticOrientation{T} end


const RadiologicalView = SemanticOrientation{(LeftToRight,AnteriorToPosterior,SuperiorToInferior)}()

const NeurologicalView = SemanticOrientation{(RightToLeft,AnteriorToPosterior,SuperiorToInferior)}()

"""
    isradiologic(x) -> Bool

Test to see if `x` is in radiological orientation.
"""
isradiologic(::SemanticOrientation{T}) where {T} = isradiologic(T)

isradiologic(x::Tuple{3,<:SemanticOrder}) = x === (LeftToRight, AnteriorToPosterior, SuperiorToInferior)

isradiologic(x::NTuple{3,Symbol}) = isradiologic(SemanticOrder.(x))
isradiologic(x::NTuple{3,AbstractString}) = isradiologic(SemanticOrder.(x))

isradiologic(x...) = isradiologic(Tuple(x))


"""
    isneurologic(x) -> Bool

Test to see if `x` is in neurological orientation.
"""
isneurologic(::SemanticOrientation{T}) where {T} = isneurologic(T)

isneurologic(x::Tuple{3,<:SemanticOrder}) = x === (LeftToRight, AnteriorToPosterior, SuperiorToInferior)

isneurologic(x::NTuple{3,Symbol}) = isneurologic(SemanticOrder.(x))
isneurologic(x::NTuple{3,AbstractString}) = isneurologic(SemanticOrder.(x))

isneurologic(x...) = isneurologic(Tuple(x))

function num2axis(i::Real)
    if i == 1 || i == -1
        return :sagittal
    elseif i == 2 || i == -2
        return :coronal
    elseif i == 3 || i == -3
        return :axial
    else
        return :unkown
    end
end


function num2ori(x::Real)
    if x == 1
        return :left_to_right  # LeftToRight
    elseif x == -1
        return :right_to_left  # RightToLeft
    elseif x == 2
        return :posterior_to_anterior  # PosteriorToAnterior
    elseif x == -2
        return :anterior_to_posterior # AnteriorToPosterior
    elseif x == 3
        return :inferior_to_superior  # InferiorToSuperior
    elseif x == -3
        return :superior_to_inferior  # SuperiorToInferior
    end
end

ori2num(::SemanticOrder{Left,Right}) = 1
ori2num(::SemanticOrder{Right,Left}) = -1

ori2num(::SemanticOrder{Posterior,Anterior}) = 2
ori2num(::SemanticOrder{Anterior,Posterior}) = -2

ori2num(::SemanticOrder{Inferior,Superior}) = 3
ori2num(::SemanticOrder{Superior,Inferior}) = -3

end
