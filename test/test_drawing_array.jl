using GLMakie
using Colors

using ArrayRadiation


# Define element spacing and array size
element_separation_λ = 1/2;
axis_element_count = 32;

# Generate antenna element positions
r_xyz = ArrayRadiation.antenna_matrix(axis_element_count, axis_element_count, element_separation_λ)


# Plot the 3D gain pattern
fig = Figure()
ax = Axis3(fig[1, 1], title="Array Gain Pattern [dB]", xlabel="x", ylabel="y", zlabel="z")



#_________________________________________________________________________________________
# Add rectangular patch elements on a rectangular substrate

visual_scale = 2.0  # Scaling factor for visualization

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

# Display the figure
fig