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
struct ImageProperties{sym} <: AbstractDict{String,Any}
    properties::AbstractDict{String,Any}

    ImageProperties(d::AbstractDict{String,Any}) = new{:NOTHING}(d)
    ImageProperties{sym}(d::AbstractDict{String,Any}) where {sym} = new{sym}(d)

    ImageProperties(d::ImageProperties{sym}) where sym = new{sym}(properties(d))
    ImageProperties{sym}(d::AbstractDict{sym}) where sym = new{sym}(properties(d))
end

ImageProperties() = ImageProperties{:Nothing}()
ImageProperties{sym}() where {sym} = ImageProperties{sym}(Dict{String,Any}())

struct IMNothing end   # to avoid confusion in the case where dict[key] === nothing
macro get(d, k, default)
    quote
        d, k = $(esc(d)), $(esc(k))
        val = get(d.properties, k, IMNothing())
        return isa(val, IMNothing) ? $(esc(default)) : val
    end
end


function ImageProperties(img::ImageMeta; copyprops::Bool=false)
    if copyprops
        ImageProperties(deepcopy(properties(img)))
    else
        ImageProperties(properties(img))
    end
end

function ImageProperties{sym}(img::ImageMeta; copyprops::Bool=false) where sym
    if copyprops
        ImageProperties{sym}(deepcopy(properties(img)))
    else
        ImageProperties{sym}(properties(img))
    end
end

ImageMetadata.properties(d::ImageProperties) = d.properties

metatype(d::ImageProperties{sym}) where {sym} = sym

Base.setindex!(d::ImageProperties, val, key::String) = setindex!(properties(d), val, key)
Base.getindex(d::ImageProperties, key::String) = properties(d)[key]

Base.copy(d::ImageProperties{sym}) where sym = ImageProperties{sym}(deepcopy(properties(d)))
Base.delete!(d::ImageProperties, key::String) = (delete!(properties(d), key); d)
# If certain dicts are suppose to have certain elements this may need customization per data type
#Base.empty!(d::ImageProperties) = (empty!(properties(d)); d)
#Base.empty(d::ImageProperties{sym}) where F = ImageProperties{sym}()
Base.isempty(d::ImageProperties) = isempty(properties(d))

Base.in(item, d::ImageProperties) = in(item, properties(d))

Base.pop!(d::ImageProperties, key, default) = pop!(properties(d), key, default)
Base.pop!(d::ImageProperties, key) = pop!(properties(d), key, Base.secret_table_token)
Base.push!(d::ImageProperties, kv::Pair) = insert!(d, kv[1], kv[2])
Base.push!(d::ImageProperties, kv) = insert!(properties(d), kv[1], kv[2])

Base.haskey(d::ImageProperties, k::String) = haskey(properties(d), k)
Base.keys(d::ImageProperties) = keys(properties(d))
Base.getkey(d::ImageProperties, key, default) = getkey(properties(d), key, default)
Base.get!(f::Function, d::ImageProperties, key) = get!(f, properties(d), key)


# iteration interface
Base.iterate(d::ImageProperties) = iterate(properties(d))
Base.iterate(d::ImageProperties, state) = iterate(properties(d), state)

Base.length(d::ImageProperties) = length(properties(d))

Base.filter!(f, d::ImageProperties) = filter!(f, properties(d))
Base.get(d::ImageProperties, k, v) = get(properties(d), k, v)
