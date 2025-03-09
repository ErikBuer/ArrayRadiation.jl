using Plots;
gr();

using ArrayRadiation


function plot_window(r, W, plot_title::String)
    scatter(r, W, 
        marker=:circle, 
        linecolor=:blue, 
        markersize=4, 
        xlabel="Element Position [λ]", 
        ylabel="Weight", 
        title=plot_title,
        legend=false, 
        grid=true,
    )

    savefig("plots/"*plot_title*".png")
end


N = 32
element_separation_λ = 1/2;
r = ArrayRadiation.linear_array(N, element_separation_λ)


W = Window.taylor(N,4,-25)
plot_window(r, W, "Taylor Window")

W = Window.split_taylor(N,4,-25)
plot_window(r, W, "Split Taylor Window")