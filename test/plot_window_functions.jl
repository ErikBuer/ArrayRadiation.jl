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


function multiplot_window(r, W::AbstractVector, legends::Vector{String}, plot_title::String)
    scatter(r, W[1], 
        marker=:circle, 
        linecolor=:blue, 
        markersize=4, 
        xlabel="Element Position [λ]", 
        ylabel="Weight", 
        title=plot_title,
        legend=true,
        grid=true,
        label = legends[1],
    )

    for i = 2:length(W)
        scatter!(r, W[i], 
        marker=:circle, 
        linecolor=:auto,
        markersize=4, 
        label=legends[i]
    )
    end

    savefig("plots/"*plot_title*".png")
end

N = 32
element_separation_λ = 1/2;
r = ArrayRadiation.linear_array(N, element_separation_λ)


W = Window.taylor(N,4,-25)
plot_window(r, W, "Taylor Window")

W = Window.split_window(W)
plot_window(r, W, "Split Taylor Window")


W0 = Window.cosine_q(N, 0)
W1 = Window.cosine_q(N, 1)
W2 = Window.cosine_q(N, 2)

multiplot_window(r, [W0, W1, W2], ["q=0", "q=1", "q=2"], "Cosine q Window")
