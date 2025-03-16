using GLMakie
using ArrayRadiation


angleRad = LinRange(π / 2, -π / 2, 81);
angleDeg = rad2deg.(angleRad);

element_gain_approximation = Kspace.cos_taper.(angleRad)
element_gain_approximation_dB = DspUtility.pow2db(element_gain_approximation)

y_lower_limit = -30 # dB
element_gain_approximation_dB = clamp.(element_gain_approximation_dB, y_lower_limit, Inf)

# Create figure and polar axis
f = Figure()
ax = PolarAxis(f[1, 1], 
    title = "Element Radiation Pattern",
    thetalimits = (-pi/2, pi/2),
    radius_at_origin = -30,
    theta_0 = -pi/2,
    direction = -1,
)

# Plot radiation pattern
lines!(ax, angleRad, element_gain_approximation_dB, color = :blue, linewidth = 2)

# Save the figure (optional)
save("plots/antenna_radiation_pattern.png", f)

f  # Display figure