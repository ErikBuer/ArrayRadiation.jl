# 2D Array

## 1D Radiation Pattern

Lets create a 2D antenna array and look at its radiation pattern.

First we must place the antenna elements. Lets give them λ/2 spacing and spread them equally.

``` @example StaticArray
using Plots;
gr();
using LaTeXStrings
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
    xlabel="X [λ]", 
    ylabel="Y [λ]", 
    title="Element Positions", 
    legend=false, 
    grid=true
)
```

These elements now represent the antenna placements.
We use the `cos_taper` function to approzimate the element gain.

``` @example StaticArray
angleRad = LinRange(π / 2, -π / 2, 501);
angleDeg = rad2deg.(angleRad);

element_gain_approximation = AntennaElement.cos_taper.(angleRad)
```

We can apply different weights to each element to create different beam directions.

But in this example we give all elements a uniform weight:

``` @example StaticArray
# Antenna element weight
W = ones(axis_element_count, axis_element_count)
```

With this defined, we can calculate the radiation pattern of the array.

The propagation direciton is decided by the vector ``\vec{k}``, as opposed to elevation and azimuth angles (``θ, ϕ``).

``||\vec{k}|| = \frac{2π}{λ_0}``

``` @example StaticArray
# Calculate K-space vectors
k_xyz = 2π*Kspace.elevation2k_hat.(angleRad)

# Map K-space gain calculation function.
GΩ(k) = Kspace.gain(k, 1, r_xyz, W)
GΩ_lin = broadcast(GΩ, k_xyz).*element_gain_approximation

GΩ_dB = DspUtility.pow2db.(abs.(GΩ_lin))

plot(angleDeg, GΩ_dB,
    xlabel = "Elevation [deg]",
    ylabel = L"G_Ω\; [dB]",
    title  = "Array gain",
    ylims  = (-20, 30),
    reuse  = true,
    legend = false
)
```

We can now inspect the radiation pattern in one dimension. This is useful to get a sense of the performance.

However one is often interrested in the complete radiaiton pattern in all dimensions.

## 2D Radiation Pattern

The element radiation pattern `AntennaElement.cos_taper` only radiate forward, so a 2D plot shows all details of the current array.

For simplicity, we calculate the radiation pattern from ``k_x, k_y ∈ [-2\pi, 2\pi]``.

``\vec{k} = k_x\hat{x} + k_y\hat{y} + k_z\hat{z}``

This yields invalid ``\vec{k}`` vectors with a magniture greater than 2π (in our case) in the corners of the plot.

But as you can see in the resulting radiation pattern, due to our element gain pattern, the array gain is nothing here anyways.

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

heatmap(x_vals, y_vals, k_z, 
    xlabel=L"\hat{k}_x",
    ylabel=L"\hat{k}_y",
    title="Values of "*L"\hat{k}_z",
    color=:jet1
)
```

``` @example StaticArray
element_gain(elevation) = AntennaElement.cos_taper.(elevation, 1.4)

# Calculate k_z values and gain values
for m in 1:resolution
    for n in 1:resolution
        k_x = x_vals[n]
        k_y = y_vals[m]
        k_z = sqrt(max(0, 1 - k_x^2 - k_y^2))
        
        # Create the k vector for each (k_x, k_y, k_z) triplet
        k_xyz = 2π.*[k_x, k_y, k_z]

        θ = Kspace.k2elevation(k_xyz)
        GΩ_lin[m, n] = GΩ(k_xyz)*element_gain(θ)
        
    end
end

# Convert the gain to dB
GΩ_dB = DspUtility.pow2db.(abs.(GΩ_lin))

GΩ_dB = clamp.(GΩ_dB, -40, Inf)

# Plot the result as a heatmap
heatmap(x_vals, y_vals, GΩ_dB, 
    xlabel=L"\hat{k}_x", 
    ylabel=L"\hat{k}_y", 
    title=L"G_Ω(\vec{k})\; [dB]",
    color=:jet1
)
```
