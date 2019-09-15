# ImageFormats

## Image I/O for Developers

The following examples demonstrate the flow of data/types that make up ImageFormats at the
I/O level. Each example starts at a level that is intended for public use and goes through
the internals. In what follows, the type parameter `F` is a FileIO registered `DataFormat`.

### Loading Images

We'll begin with the simplest case and slowly increase complexity that allows for flexible
needs and output. We assume two things in most of what follows

1. Some sort of image header info preceeds the actual image data
2. Image data is of a single `Type` (note that I mean `Type` in the Julia sense.
   An image of colors specified by some sort of float must use the same floats throughout.)

Let's begin with the simplest case and implementation, reading in an image that is
`Images.jl` compatible. We need two functions to accomplish this, a `Stream` loader and
something that reads the metadata/header for each image.

#### Super Simple Implementation

The `Stream` loader is for `FileIO.jl` compatibility and can be as simple as:
```julia
load(s::Stream{F}, args...; kwargs...) where F = read(ImageStream(s; kwargs...), args...; kwargs...)
```

Reading metadata is a bit more complicated, as it requires a unique method for each format.
For the sake of simplicity, we'll leave it at this for now:
```julia
function ImageFormats.unsafe_metadata(s::Stream{F}, args...; kwargs...) where F
    ...
end --> ImageInfo
```
Overloading the `unsafe_metadata` requires two things. 1) You must have a `FileIO.jl`
registered format for this to integrate properly with everything else. 2) It must return an
`ImageInfo` object. For more details on what this accomplishes see later sections on
implementing the `metadata` method for your package/format.

#### Extensible Implementation

Now we assume the previous methods have been implemented but we want to implement some
extra features. Here we implement directly loading from a file string, let the user decide
the mode of the stream produced, and a default sink.

```julia
function load(f::File{F}, sink=Array; mode::String="r", kwargs...) where F
    open(f, mode) do s
        load(s, sink; kwargs...)
    end
end
```
The `mode` argument isn't unique to `ImageFormats`. Any package write the above function.
However, we also add the `sink` argument, which dictates the type of array produced. Now
the default array produced by the format `F` will be an `Array`. If you don't want to change
the default (which is `ImageMeta{T,N,AxisArray{T,N,Ax}}`) then you don't have to do
anything. The user will still have the option of changing the sink, as long as it is the
second positional argument in the `read(::ImageStream, sink)` method.

### Loading Metadata

If you'd like to use the `FileIO` method  `metadata` for just loading information about the
image you can implement the following functions.

```julia
metadata(f::File{F}, args...; kwargs...) where F = ImageInfo(f, args...; kwargs...)

metadata(s::Stream{F}, args...; kwargs...) where F = ImageInfo(s, args...; kwargs...)

function ImageFormats.unsafe_metadata(s::Stream{F}, args...; kwargs...) where F
    ...
end --> ImageInfo
```
Note that the first two methods are simply allowing `FileIO` registered functions to call
`ImageInfo`, while the last method is exactly the same as the one used previously for
loading images. It only needs to be implemented once, but is shown here for clarity.


## Internal Methods

Internally, both `ImageStream` and `ImageInfo` rely on `substream`. The only purpose of
substream is to ensure that the `File` or `Stream` provided is registerred and to
decompress a stream, if necessary. Because there currently is no native support for passing
compressed files through a decompression and then a specified loader, `substream` takes
whatever `Stream` it's fed and appropriately decompresses it.

```julia
function load(f::File{F}, args...; kwargs...) where F # <- registered format
    open(f) do s
        load(s, args...; kwargs...)
    end
end

load(f::Stream{F}, args...; kwargs...) where F = read(ImageStream(f))


```

### Loading Metadata


The `metadata` function is implemented in `FileIO.jl`. `ImageFormats` uses this convention
to construct the `ImageInfo` type. Packages

`metadata(f::File{F}) -> open(f) -->metadata(open(f))

1. We know it's not in the final state and decompress
 `ImageInfo(s::Stream{DataFormat{:GZIP}} =  unsafe_metadata(subformat(s))`
2. No idea if in final state, superficially looks like:
  `ImageInfo(s::Stream{DataFormat{F}} =  unsafe_metadata(subformat(s))`
   but we check extensions internally for consistancy with Stream
## Introduction

## Roadmap

* [] Common interface to read various image formats
* [] Ease reading to specific array types with a `load(file, sink)` interface
* [] Chunked reading of on disk images
* [] Accelerated reading of images that share important properties
    * The biggest obstacle here will be designing and interface that bypasses needlessly reading header information from images the are already processed and have the same info, but still is safe for most users to use without accidently causing out of memory read errors
    * Could be mediated by only checking file size, which should be the same for all if they share dimensions and element types
* [] Iterator interface for accessing and manipulating sets of images
