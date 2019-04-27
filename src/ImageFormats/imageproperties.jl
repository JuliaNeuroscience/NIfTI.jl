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
    data::AbstractDict{String,Any}

    ImageProperties(d::AbstractDict{String,Any}) = new{:NOTHING}(d)
    ImageProperties{sym}(d::AbstractDict{String,Any}) where {sym} = new{sym}(d)
end

ImageProperties() = ImageProperties{:Nothing}()
ImageProperties{sym}() where {sym} = ImageProperties{sym}(Dict{String,Any}())

function ImageProperties(ImageMeta)
end

metatype(d::ImageProperties{sym}) where {sym} = sym

Base.setindex!(d::ImageProperties, val, key::String) = setindex!(d.data, val, key)
Base.getindex(d::ImageProperties, key::String) = d.data[key]

Base.copy(d::ImageProperties{sym}) where sym = ImageProperties{sym}(deepcopy(d.data))
Base.delete!(d::ImageProperties, key::String) = (delete!(d.data, key); d)
# If certain dicts are suppose to have certain elements this may need customization per data type
#Base.empty!(d::ImageProperties) = (empty!(d.data); d)
#Base.empty(d::ImageProperties{sym}) where F = ImageProperties{sym}()
Base.isempty(d::ImageProperties) = isempty(d.data)

Base.in(item, d::ImageProperties) = in(item, d.data)

Base.pop!(d::ImageProperties, key, default) = pop!(d.data, key, default)
Base.pop!(d::ImageProperties, key) = pop!(d.data, key, Base.secret_table_token)
Base.push!(d::ImageProperties, kv::Pair) = insert!(d, kv[1], kv[2])
Base.push!(d::ImageProperties, kv) = insert!(d.data, kv[1], kv[2])

Base.haskey(d::ImageProperties, k::String) = haskey(d.data, k)
Base.keys(d::ImageProperties) = keys(d.data)
Base.getkey(d::ImageProperties, key, default) = getkey(d.data, key, default)
Base.get!(f::Function, d::ImageProperties, key) = get!(f, d.data, key)

Base.get(d::ImageProperties, k, v) = get(d.data, k, v)

# iteration interface
Base.iterate(d::ImageProperties) = iterate(d.data)
Base.iterate(d::ImageProperties, state) = iterate(d.data, state)

Base.length(d::ImageProperties) = length(d.data)

Base.filter!(f, d::ImageProperties) = filter!(f, d.data)
