## Get started

```@example demo
using NIfTI
using Downloads, FileIO
using GLMakie

file = "https://raw.githubusercontent.com/nipraxis/nipraxis-data/0.5/ds114_sub009_highres.nii"
ds_brain = Downloads.download(file, "ds114_sub009_highres.nii")
img = niread(ds_brain);

cmap = :linear_protanopic_deuteranopic_kbw_5_98_c40_n256
n = 202
g(x) = x^2
alphas = [g(x) for x in range(0.0, 1, length = n)]
cmap_alpha = resample_cmap(cmap, n; alpha = alphas)

x = LinRange(1, size(img, 1), size(img, 1))
y = LinRange(1, size(img, 2), size(img, 2))
z = LinRange(1, size(img, 3), size(img, 3))

with_theme(theme_light()) do
    fig = Figure(; size = (1200, 1200), backgroundcolor = :black)
    g = GridLayout(fig[1,1])
    ax = LScene(g[1, 1]; show_axis = false)
    # slices grid
    sgrid = SliderGrid(
        g[2, 1],
        (label = "yz plane - x axis", range = 1:length(x)),
        (label = "xz plane - y axis", range = 1:length(y)),
        (label = "xy plane - z axis", range = 1:length(z)),
    )

    lo = sgrid.layout
    nc = ncols(lo)

    plt = volume!(ax, 1..size(img,1), 1..size(img,2), 1..size(img,3), Array(img);
        colorrange = (1, 3000), lowclip=:transparent, highclip=:white,
        colormap = cmap_alpha, transparency = true)
    # add volume slices    
    plt_slices = volumeslices!(ax, x, y, z, Array(img); colormap = cmap,
        colorrange = (10, 3000), highclip=:white,
        # transparency = true
        )

    # connect sliders to `volumeslices` update methods
    sl_yz, sl_xz, sl_xy = sgrid.sliders

    on(sl_yz.value) do v; plt_slices[:update_yz][](v) end
    on(sl_xz.value) do v; plt_slices[:update_xz][](v) end
    on(sl_xy.value) do v; plt_slices[:update_xy][](v) end

    set_close_to!(sl_yz, .5length(x))
    set_close_to!(sl_xz, .5length(y))
    set_close_to!(sl_xy, .5length(z))

    # add toggles to show/hide heatmaps
    hmaps = [plt_slices[Symbol(:heatmap_, s)][] for s ∈ (:yz, :xz, :xy)]
    toggles = [Toggle(lo[i, nc + 1], active = true) for i ∈ 1:length(hmaps)]
    # labels = ["Show yz slice", "Show xz slice", "Show xy slice"]
    hm1 = plot(fig[1,2], hmaps[1][3]; colormap = cmap, colorrange = (1, 3000),
        lowclip=:transparent, highclip=:white, transparency = true)
    plot(fig[2,2], hmaps[2][3]; colormap = cmap, colorrange = (1, 3000),
        lowclip=:transparent, highclip=:white, transparency = true)
    plot(fig[2,1], hmaps[3][3]; colormap = cmap, colorrange = (1, 3000),
        lowclip=:transparent, highclip=:white, transparency = true)

    map(zip(hmaps, toggles)) do (h, t)
        on(t.active) do active
            h.visible = active
        end
    end
    fig
end
```