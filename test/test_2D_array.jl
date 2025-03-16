using Plots; gr()
using ArrayRadiation

element_separation_λ = 1/2;

# Place elements symmetrically around zero
axis_element_count = 32;

r_xyz = ArrayRadiation.antenna_matrix(axis_element_count, axis_element_count, element_separation_λ)


# Flatten and extract x and y coordinates
x_positions = [r[1] for r in vec(r_xyz)]
y_positions = [r[2] for r in vec(r_xyz)]

# Plot
scatter(x_positions, y_positions, 
    markershape=:circle, 
    markersize=4, 
    xlabel="X Position [λ]", 
    ylabel="Y Position [λ]", 
    title="Antenna Element Positions", 
    legend=false, 
    grid=true
)

#_______________________________________________________________
angleRad = LinRange(π / 2, -π / 2, 501);
angleDeg = rad2deg.(angleRad);

element_gain_approximation = Kspace.cos_taper.(angleRad)

#_______________________________________________________________


# Antenna element weight
W = ones(axis_element_count, axis_element_count)

#_______________________________________________________________
# Calculate K-space vectors
k_xyz = 2π*Kspace.elevation2k_hat.(angleRad)

# Map K-space gain calculation function.
GΩ(k) = Kspace.gain(k, 1, r_xyz, W)
GΩ_lin = broadcast(GΩ, k_xyz).*element_gain_approximation

GΩ_dB = DspUtility.pow2db.(abs.(GΩ_lin))

plot(angleDeg, GΩ_dB,
    xlabel = "Elevation [deg]",
    ylabel = "G_Ω [dB]",
    title  = "Array gain",
    ylims  = (-20, 30),
    reuse  = true,
    legend = false
)
#_______________________________________________________________
