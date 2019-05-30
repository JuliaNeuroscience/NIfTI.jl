#= quick thought

struct ImageProperties{S} <: AbstractDict{String,Any}
    header::H
    spatial::S
end

=#
"""
    ImageProperties

```jldoctest
julia> ImageProperties()
ImageProperties{:Nothing} with 0 entries

julia> ImageProperties(Dict{String,Any}("important_stuff" => 1))
ImageProperties{:NOTHING} with 1 entry:
  "important_stuff" => 1

julia>
```
"""
struct ImageProperties{F} <: AbstractDict{String,Any}
    properties::AbstractDict{String,Any}

    function ImageProperties(d::AbstractDict{String,Any}; copyprops::Bool=false)
        if copyprops
            new{format"NOTHING"}(deepcopy(d))
        else
            new{format"NOTHING"}(d)
        end
    end

    function ImageProperties{F}(d::AbstractDict{String,Any}; copyprops::Bool=false) where {F}
        if copyprops
            new{F}(deepcopy(d))
        else
            new{F}(d)
        end
    end

    function ImageProperties{Fnew}(d::ImageProperties{F}; copyprops::Bool=false) where {F,Fnew}
        if copyprops
            new{F}(deepcopy(d))
        else
            new{F}(d)
        end
    end

    function ImageProperties(d::ImageProperties{F}; copyprops::Bool=false) where F
        if copyprops
            new{F}(deepcopy(properties(d)))
        else
            new{F}(properties(d))
        end
    end
end

ImageProperties() = ImageProperties{format"Nothing"}()
ImageProperties{F}() where {F} = ImageProperties{F}(Dict{String,Any}())

struct IMNothing end   # to avoid confusion in the case where dict[key] === nothing
macro get(p, k, default)
    quote
        p, k = $(esc(p)), $(esc(k))
        val = get(p.properties, k, IMNothing())
        return isa(val, IMNothing) ? $(esc(default)) : val
    end
end

getheader(A::AbstractArray, k::String, default) = default
getheader(img::ImageMeta{T,N,A,<:ImageProperties}, k::String, default) where {T,N,A} =
    getheader(properties(img), k, default)
getheader(p::ImageProperties, k::String, default) = getheader(header(p), k, default)
getheader(p::ImageProperties{:header}, k::String, default) = @get p k default
getheader(::Nothing, k::String, default) = default

function ImageProperties(img::ImageMeta; copyprops::Bool=false)
    if copyprops
        _add_spacedirections(ImageProperties(deepcopy(properties(img))), img)
    else
        _add_spacedirections(ImageProperties(properties(img)), img)
    end
end

function ImageProperties{F}(img::ImageMeta; copyprops::Bool=false) where F
    if copyprops
        _add_spacedirections(ImageProperties{F}(deepcopy(properties(img))), img)
    else
        _add_spacedirections(ImageProperties{F}(properties(img)), img)
    end
end

ImageProperties(a::AbstractArray, props::AbstractDict{String,Any}=Dict{String,Any}(); copyprops::Bool=false) =
    _add_spacedirections(ImageProperties(props), a)

ImageProperties{F}(a::AbstractArray, props::AbstractDict{String,Any}=Dict{String,Any}(); copyprops::Bool=false) where F =
    _add_spacedirections(ImageProperties{F}(props), a)


function _add_spacedirections(p::ImageProperties, A::AbstractArray)
    if haskey(p, "spacedirections")
        return p
    else
        p["spacedirections"] = spacedirections(A)
        return p
    end
end

ImageMetadata.properties(d::ImageProperties) = d.properties

metatype(d::ImageProperties{F}) where {F} = F

Base.setindex!(d::ImageProperties, val, key::String) = setindex!(properties(d), val, key)
Base.getindex(d::ImageProperties, key::String) = properties(d)[key]

Base.copy(d::ImageProperties{F}) where F = ImageProperties{F}(deepcopy(properties(d)))
Base.delete!(d::ImageProperties, key::String) = (delete!(properties(d), key); d)
# If certain dicts are suppose to have certain elements this may need customization per data type
#Base.empty!(d::ImageProperties) = (empty!(properties(d)); d)
#Base.empty(d::ImageProperties{F}) where F = ImageProperties{F}()
Base.isempty(d::ImageProperties) = isempty(properties(d))

Base.haskey(d::ImageProperties, k::String) = haskey(properties(d), k)
Base.keys(d::ImageProperties) = keys(properties(d))
Base.getkey(d::ImageProperties, key::String, default) = getkey(properties(d), key, default)
Base.get!(f::Function, d::ImageProperties, key::String) = get!(f, properties(d), key)

# iteration interface
Base.iterate(d::ImageProperties) = iterate(properties(d))
Base.iterate(d::ImageProperties, state) = iterate(properties(d), state)

Base.length(d::ImageProperties) = length(properties(d))

Base.get(d::ImageProperties, k::String, v) = get(properties(d), k, v)

Base.filter!(f, d::ImageProperties) = filter!(f, properties(d))

Base.pop!(d::ImageProperties, key::String, default) = pop!(properties(d), key, default)
Base.pop!(d::ImageProperties, key::String) = pop!(properties(d), key, Base.secret_table_token)
Base.push!(d::ImageProperties, kv::Pair) = insert!(d, kv[1], kv[2])
Base.push!(d::ImageProperties, kv) = insert!(properties(d), kv[1], kv[2])
Base.in(item, d::ImageProperties) = in(item, properties(d))
