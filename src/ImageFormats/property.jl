struct MissingProperty end
const MProperty = MissingProperty()

"""
    getter(x, key[, property_type, default_value])
"""
function getter(x::X, key::Any, property_type::Any, default_value::Any) where {X}
    filter_getproperty(x, key, property_type, default_value, _getter(HasProperties(X), x, key))
end
_getter(::HasProperties{true}, x::Any, key::Any) = get(properties(x), key, MProperty)
_getter(::HasProperties{false}, x::Any, key::Any) = MProperty
_getter(::HasProperties{false}, x::AbstractDict{K,Any}, key::K) where {K} = get(x, key, MProperty)

function filter_getproperty(
    x::Any,
    key::Any,
    property_type::Type{T},
    default_value::Any,
    value::V
   ) where {T,V<:T}
    return value
end

function filter_getproperty(
    x::Any,
    key::Any,
    property_type::Type{T},
    default_value::Any,
    value::V
   ) where {T,V}
    return convert(T, value)
end

function filter_getproperty(
    x::Any,
    key::Any,
    property_type::Function,
    default_value::Any,
    value::V
   ) where {T,V}
    filter_getproperty(x, key, property_type(x), default_value, value)
end

function filter_getproperty(
    x::Any,
    key::Any,
    property_type::Type,
    default_value::Any,
    value::MissingProperty
   )
    return default_value
end

# if `default_value` is a function then run `x` through it and filter again to
# ensure it is `<: property_type`
function filter_getproperty(
    x::Any,
    key::Any,
    property_type::Type,
    default_value::Function,
    value::MissingProperty
   )
    return filter_getproperty(x, key, property_type, default_value, default_value(x))
end

function filter_getproperty(
    x::Any,
    key::Any,
    property_type::Function,
    default_value::Function,
    value::MissingProperty
   )
    return filter_getproperty(x, key, property_type, default_value, default_value(x))
end



"""
    setter(x, key, value[, property_type])

"""
function setter!(x::X, key::Any, value::Any, property_type::Any) where {X}
    _setter(HasProperties(X), x, key, filter_setproperty(x, key, value, property_type))
end
function _setter(::HasProperties{true}, x::Any, key::Any, value::Any)
    setindex!(properties(x), value, key)
end

function _setter(::HasProperties{false}, x::AbstractDict{K,Any}, key::K, value::Any) where {K}
    setindex!(x, value, key)
end

function _setter(::HasProperties{false}, x::Any, key::Any, value::Any)
    error("$(typeof(x)) does not have `properties` method. Cannot set $key property.")
end

function filter_setproperty(
    x::Any,
    key::Any,
    value::V,
    property_type::Type{T},
   ) where {T,V<:T}

    return value
end

function filter_setproperty(
    x::Any,
    key::Any,
    value::V,
    property_type::Type{T},
   ) where {T,V}

    return convert(T, value)
end

function filter_setproperty(
    x::Any,
    key::Any,
    value::Any,
    property_type::Function,
   ) where {T,V}

    return filter_setproperty(x, key, value, property_type(x))
end


