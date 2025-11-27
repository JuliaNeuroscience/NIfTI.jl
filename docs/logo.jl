using NIfTI
using Downloads, FileIO
using GLMakie
GLMakie.activate!()

file = "https://raw.githubusercontent.com/nipraxis/nipraxis-data/0.5/ds114_sub009_highres.nii"
ds_brain = Downloads.download(file, "ds114_sub009_highres.nii")
img = niread(ds_brain);

# cmap = :linear_protanopic_deuteranopic_kbw_5_98_c40_n256
cmap = :cmr_cosmic
n = 202
g(x) = x^2
alphas = [g(x) for x in range(0.0, 1, length = n)]
cmap_alpha = resample_cmap(cmap, n; alpha = alphas)
# cmap_alpha = cmap_alpha[80:end]
slice_d = Array(img)[:, :, size(img, 3) ÷ 2]
x = LinRange(1, size(img, 1), size(img, 1))
y = LinRange(1, size(img, 2), size(img, 2))
z = LinRange(1, size(img, 3), size(img, 3))

with_theme(theme_light()) do
    fig = Figure(; size = (600, 600), backgroundcolor = :black)
    ax = LScene(fig[1, 1]; show_axis = false)
    plt = volume!(ax, 1..size(img,1), 1..size(img,2), 1..size(img,3), Array(img);
        colorrange = (1, 3000),
        lowclip=:transparent, highclip=:white, colormap = cmap_alpha, transparency = true)
    plt_slices = volumeslices!(ax, x, y, z, Array(img); colormap = cmap,
        colorrange = (120, 3000), highclip=:white, lowclip=(:grey25, 0.1),
        bbox_visible = false,
        transparency = false
        )

    plt_slices[:update_yz][](length(x)÷2)
    plt_slices[:update_xz][](length(y)÷2)
    plt_slices[:update_xy][](length(z)÷2)

    fig
end
save(joinpath(@__DIR__, "src/assets/logo.png"), current_figure(), px_per_unit=4, update=false)
