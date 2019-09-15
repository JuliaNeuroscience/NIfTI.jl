
struct NodeIndex{Dim,T,Element,A} <: AbstractVector{Element}
    nodes::Mesh{Dim,T,Element,A}
    space::CoordinateSpace
end

Base.getindex(ni::NodeIndex, i) = ni.nodes[i]
Base.size(ni::NodeIndex) = size(ni.points)
