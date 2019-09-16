# These convert between axis arrays format and namedaxes.
namedaxes2axisarray(x::NamedTuple{Names}) where {Names} = (Axis{first(Names)}(first(x)), namedaxes2axisarray(tail(x))...)
namedaxes2axisarray(x::NamedTuple{(),Tuple{}}) = ()

ImageCore.namedaxes(a::AxisArray) = NamedTuple{axisnames(a)}(axisvalues(a))

ImageCore.HasProperties(::Type{<:ImageMeta}) = HasProperties{true}()

## Time dimension
timedim(x::Union{ArrayStream,ArrayInfo}) = finddim(x, :time)

timeaxis(x::Union{ArrayStream,ArrayInfo}) = findaxis(x, :time)

colordim(x::Union{ArrayStream,ArrayInfo}) = finddim(x, :color)

