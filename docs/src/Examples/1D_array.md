# 1D Array

## 1D Radiation Pattern

Lets create an array and look at its radiation pattern.

First we must place the antenna elements. Lets give them λ/2 spacing and spread them linearly.

``` @example StaticArray
using Plots;
gr();
using LaTeXStrings

using ArrayRadiation

element_separation_λ = 1/2;

# Place elements symmetrically around zero
element_count = 32;

r = ArrayRadiation.linear_array(element_count, element_separation_λ)

scatter(r, zeros(length(r)), 
    markershape=:circle, 
    # markercolor=:blue, 
    markersize=4, 
    xlabel="Element Position [λ]", 
    title="Antenna Element Positions", 
    legend=false, 
    #grid=false,
    yticks= false,
    ylims  = (-1, 1),
)
```

These elements now represent the antenna placements.
For now, we use the `cos_taper` funciton to approzimate the element gain.

``` @example StaticArray
angleRad = LinRange(π / 2, -π / 2, 501);
angleDeg = rad2deg.(angleRad);

element_gain_approximation = Kspace.cos_taper.(angleRad)
```

We can apply different weights to each element so lets look at how this will affect our array.

Lets start by giving each element a uniform weight:

``` @example StaticArray
# Antenna element weigth
W = ones(element_count)

scatter(r, W, 
    marker=:circle, 
    linecolor=:blue, 
    markersize=4, 
    xlabel="Element Position [λ]", 
    ylabel="Weight", 
    title="Antenna Element Weights", 
    legend=false, 
    grid=true,
    ylims  = (0, 1.1),
)
```

With this defined, we can calculate the radiation pattern of the array.

``` @example StaticArray
# Calculate K-space vectors
k_xyz = 2π*Kspace.elevation2k_hat.(angleRad)

# Map K-space gain calculation function.
GΩ(k) = Kspace.gain_1D(k, 1, r, W)
GΩ_lin = broadcast(GΩ, k_xyz).*element_gain_approximation

GΩ_dB = DspUtility.pow2db.(abs.(GΩ_lin))

plot(angleDeg, GΩ_dB,
    xlabel = "Elevation [deg]",
    ylabel = L"G_Ω\; [dB]",
    title  = "Array gain",
    ylims  = (-30, 18),
    reuse  = true,
    legend = false
)

```

We can now inspect the radiation pattern in one dimension. This is useful to get a sense of the performance.

However one is often interrested in the complete radiaiton pattern in 3 dimensions.

## 2D Radiation Pattern

``` @example StaticArray

resolution = 201

# Create a 2D grid for x and y values
x_vals = LinRange(-1, 1, resolution)
y_vals = LinRange(-1, 1, resolution)

# Initialize matrices to store the k-values
k_z = zeros(resolution, resolution)
GΩ_lin = zeros(resolution, resolution)

for m in 1:resolution
    for n in 1:resolution
        k_z[m, n] = sqrt(max(0, 1 - x_vals[n]^2 - y_vals[m]^2))
    end
end

heatmap(x_vals, y_vals, k_z, xlabel=L"\hat{k}_x", ylabel=L"\hat{k}_y", title=L"\hat{k}_z", color=:jet1)

```

``` @example StaticArray

element_gain(elevation) = Kspace.cos_taper.(elevation, 1.4)

# Calculate k_z values and gain values
for m in 1:resolution
    for n in 1:resolution
        k_x = x_vals[n]
        k_y = y_vals[m]
        k_z = sqrt(max(0, 1 - k_x^2 - k_y^2))
        
        # Create the k vector for each (k_x, k_y, k_z) triplet
        k_xyz = 2π.*[k_x, k_y, k_z]

        θ = Kspace.k2elevation(k_xyz)
        
        # Calculate the gain using the provided Kspace.gain_1D function
        GΩ_lin[m, n] = GΩ(k_xyz)*element_gain(θ)
        
    end
end

# Convert the gain to dB
GΩ_dB = DspUtility.pow2db.(abs.(GΩ_lin))

GΩ_dB = clamp.(GΩ_dB, -40, Inf)

# Plot the result as a heatmap
heatmap(x_vals, y_vals, GΩ_dB, xlabel=L"\hat{k}_x", ylabel=L"\hat{k}_y", title=L"G_Ω(\vec{k})\; [dB]", color=:jet1)

```

