# 2D Array 3D Plot

``` @example 3dPlot
using GLMakie
using Colors

using ArrayRadiation

# Define element spacing and array size
element_separation_λ = 1/2;
axis_element_count = 32;

# Generate antenna element positions
r_xyz = ArrayRadiation.antenna_matrix(axis_element_count, axis_element_count, element_separation_λ)

# Define angles in spherical coordinates
θs = range(0, π/2, length=361)  # Elevation angles
φs = range(-π, π,  length=361)  # Azimuth angles
angles = [(φ,θ) for φ in φs, θ in θs]

# Define element gain function
_element_gain(elevation) = AntennaElement.cos_taper(elevation, 1.4)

# Compute gain for each angle
_GΩ(θ, φ) = Kspace.gain(Kspace.k_xyz(θ, φ), AntennaElement.cos_taper(θ), r_xyz, ones(axis_element_count, axis_element_count))

# Compute k-space vectors and gain values
rs = [_GΩ(θ, φ) for (φ,θ) in angles]

# Convert the gain to dB
rs_dB = DspUtility.pow2db.(abs.(rs))

min_value_dB = -20
rs_dB = clamp.(rs_dB, min_value_dB, Inf)
rs_dB = rs_dB .- min_value_dB # Ensure there only are possitive values.

spherical_mesh = [(r,φ,θ) for (r,(φ,θ)) in zip(rs_dB, angles)]

# Convert spherical to Cartesian coordinates
xs = [r * sin(θ) * cos(φ) for (r, φ, θ) in spherical_mesh]
ys = [r * sin(θ) * sin(φ) for (r, φ, θ) in spherical_mesh]
zs = [r * cos(θ)          for (r, φ, θ) in spherical_mesh]

max_range = 1 + maximum(abs, vcat(xs, ys, zs))  # Find max value across all coordinates


```

In order to plot the radiation pattern in 3D and dB, we need to clip the lowest values, and bias the gain so it is purely positive.
This makes the axis in the plot biased, but the relative scaling remain intact.

``` @example 3dPlot
# Plot the 3D gain pattern
fig = Figure()
ax = Axis3(fig[1, 1], title="Array Gain Pattern [dB]")
limits!(ax, -max_range/2, max_range/2, -max_range/2, max_range/2, 0, max_range)
plt = surface!(ax, xs, ys, zs, color=rs_dB, colormap=:jet1)

fig
```

To give some sense of the radiation pattern in relation to the array, lets draw a generic patch antenna array.

``` @example 3dPlot

visual_scale = 3.0  # Scaling factor for visualization

patch_width  = 0.3163 * visual_scale # wavelengths
patch_length = 0.2207 * visual_scale # wavelengths

substrate_width  = axis_element_count * element_separation_λ * visual_scale
substrate_length = axis_element_count * element_separation_λ * visual_scale
substrate_thickness = 0.06   * visual_scale # wavelengths


# Draw substrate with teflon color (light off-white)
substrate_color = RGB(0.9, 0.9, 0.9)  # Light gray/teflon color
# Copper color for the patches
copper_color = RGB(0.72, 0.45, 0.2)  # Copper color (reddish-brown)

function rectangle(center_xyz::Vector, width::Float64, height::Float64)
    # Extract the center coordinates
    x_center, y_center, z_center = center_xyz

    # Define the half-width and half-height
    half_width = width / 2
    half_height = height / 2

    # Define the four corners of the rectangle (in the XY plane)
    vertices = [
        x_center - half_width  y_center - half_height z_center;  # Point 1
        x_center + half_width  y_center - half_height z_center;  # Point 2
        x_center + half_width  y_center + half_height z_center;  # Point 3
        x_center - half_width  y_center + half_height z_center   # Point 4
    ]

    # Define the two triangles (faces) making up the rectangle
    faces = [
        1 2 3;  # First triangle
        3 4 1   # Second triangle
    ]

    return (vertices, faces)
end


# Plot the substrate at z = 0
substrate_center = [0.0, 0.0, 0.0]
vertices_substrate, faces_substrate = rectangle(substrate_center, substrate_width, substrate_length)
mesh!(ax, vertices_substrate, faces_substrate, color=substrate_color, shading=NoShading)

# Plot each patch at positions specified in r_xyz
for element_pos in r_xyz
    _x, _y, _ = element_pos .* visual_scale
    # Each patch is placed at (x, y, 0) from r_xyz
    patch_center = [_x, _y, substrate_thickness]
    vertices_patch, faces_patch = rectangle(patch_center, patch_width, patch_length)
    
    # Plot the patch with copper color
    mesh!(ax, vertices_patch, faces_patch, color=copper_color, shading=NoShading)
end


# Hide axis
hidespines!(ax)
hidedecorations!(ax)
hidedecorations!(ax, ticks = false)
hidedecorations!(ax, grid = false)


# Display the figure
fig
```
