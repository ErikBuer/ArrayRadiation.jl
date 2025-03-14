using CairoMakie

using ArrayRadiation


angleRad = LinRange(π / 2, -π / 2, 501);
angleDeg = rad2deg.(angleRad);

k_xyz = 2π*Kspace.elevation2k_hat.(angleRad)

element_gain_approximation = Kspace.cos_taper.(angleRad)


# Polar plot
fig = Figure()
ax = PolarAxis(fig[1, 1], thetalimits=(-π/2, π/2), radius_at_origin=-12, title="Radiation Pattern")

# Plot the data
lines!(ax, plotAngleRad, DspUtility.pow2db.(abs.(element_gain_approximation)))
save("plots/antenna_radiation_pattern.png", fig)