# Linear Array Polar Gain Pattern

``` @example StaticArray
using GLMakie
using ArrayRadiation

element_separation_λ = 1/2;

# Place elements symmetrically around zero
element_count = 32;

r = ArrayRadiation.linear_array(element_count, element_separation_λ)

angleRad = LinRange(π / 2, -π / 2, 501);
angleDeg = rad2deg.(angleRad);

element_gain_approximation = AntennaElement.cos_taper.(angleRad)

# Uniform weights
W = ones(element_count)


# Calculate K-space vectors
k_xyz = 2π*Kspace.elevation2k_hat.(angleRad)

# Map K-space gain calculation function.
GΩ(k) = Kspace.gain_1D(k, 1, r, W)
GΩ_lin = broadcast(GΩ, k_xyz).*element_gain_approximation

GΩ_dB = DspUtility.pow2db.(abs.(GΩ_lin))

# Clamp dB values for readability
y_lower_limit = -30
GΩ_dB = clamp.(GΩ_dB, y_lower_limit, Inf)

f = Figure()
ax = PolarAxis(f[1, 1], 
    title = "Linear Array, Uniform Weights, $element_count Elements",
    thetalimits = (-pi/2, pi/2),
    radius_at_origin = y_lower_limit,
    theta_0 = -pi/2,
    direction = -1,
)
lines!(ax, angleRad, GΩ_dB, color = :blue, linewidth = 2)
f  # Display figure
```
