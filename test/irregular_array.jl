using GLMakie
using LaTeXStrings
using ArrayRadiation
using LinearAlgebra

# Create output directory for plots
output_dir = "test/output"
mkpath(output_dir)

# Define 5 elements with different positions and orientations
# Element positions [x, y, z] in wavelengths
element_positions = [
    [0.0, 0.0, 0.0],        # Element 1: center, facing forward (z-direction)
    [0.5, 0.0, 0.0],        # Element 2: +x direction
    [0.0, 0.5, 0.0],        # Element 3: +y direction
    [-0.5, 0.0, 0.0],       # Element 4: -x direction
    [0.0, -0.5, 0.0]        # Element 5: -y direction
]

# Convert to matrix form (5×1 matrix of 3D vectors)
r_xyz = reshape(element_positions, 5, 1)

# Define antenna pointing directions (normals) for each element
element_normals = [
    [0.0, 0.0, 1.0],        # Element 1: +z direction (forward)
    [1.0, 0.0, 0.0],        # Element 2: +x direction
    [0.0, 1.0, 0.0],        # Element 3: +y direction
    [-1.0, 0.0, 0.0],       # Element 4: -x direction
    [0.0, -1.0, 0.0]        # Element 5: -y direction
]

# Gain of element at element_index in direction k_xyz
function Ge(k_xyz::AbstractVector{<:Real}, element_index::Int)::Real
    # Normalize k-space vector to get direction
    k_hat = k_xyz / norm(k_xyz)
    normal = element_normals[element_index]

    # Calculate angle between k direction and element normal
    cos_angle = dot(k_hat, normal)
    θ = acos(clamp(cos_angle, -1.0, 1.0))  # Angle in radians

    return AntennaElement.cardioid(θ)
end

# Antenna element weights 
W = ones(5, 1)

# Plot element positions and orientations
x_pos = [r[1] for r in element_positions]
y_pos = [r[2] for r in element_positions]

# Create arrows for element normals
u_arrows = [n[1] * 0.3 for n in element_normals]
v_arrows = [n[2] * 0.3 for n in element_normals]

fig1 = Figure(size=(600, 600))
ax1 = Axis(fig1[1, 1],
    xlabel="X [λ]",
    ylabel="Y [λ]",
    title="Element Positions and Orientations",
    aspect=DataAspect()
)
scatter!(ax1, x_pos, y_pos, markersize=15, color=:blue)
arrows2d!(ax1, x_pos, y_pos, u_arrows, v_arrows, color=:red, shaftwidth=2, tipwidth=0.05, tiplength=0.05)
save(joinpath(output_dir, "element_positions.png"), fig1)

# 1D elevation scan (XZ-plane cut, azimuth = 0)
elevation_angles = LinRange(π / 2, -π / 2, 501)  # From +90° to -90°
azimuth = 0.0  # Looking in XZ-plane

# Convert elevation and azimuth angles to k-space vectors using Kspace.k_xyz
k_xyz_1D = [Kspace.k_xyz(θ, azimuth) for θ in elevation_angles]


# Calculate array gain at each k-space point using the element gain function
GΩ_1D = [Kspace.gain(k, Ge, r_xyz, W) for k in k_xyz_1D]
GΩ_dB_1D = DspUtility.pow2db.(GΩ_1D)
elevation_deg = rad2deg.(elevation_angles)

fig2 = Figure(size=(800, 500))
ax2 = Axis(fig2[1, 1],
    xlabel="Elevation [deg]",
    ylabel=L"G_Ω\; [dB]",
    title="Array Gain with Cardioid Element Patterns",
    limits=(nothing, nothing, -20, 20)
)
lines!(ax2, elevation_deg, GΩ_dB_1D, linewidth=2, color=:blue)
save(joinpath(output_dir, "elevation_gain.png"), fig2)

# Polar plot in elevation plane (XZ-plane, azimuth = 0°)
# Scan full circle: elevation from 0° to 360°
elevation_angles_polar = LinRange(0, 2π, 361)
azimuth_polar_elev = 0.0

k_xyz_elev_polar = [Kspace.k_xyz(θ, azimuth_polar_elev) for θ in elevation_angles_polar]
GΩ_elev_polar = [Kspace.gain(k, Ge, r_xyz, W) for k in k_xyz_elev_polar]
GΩ_dB_elev_polar = DspUtility.pow2db.(GΩ_elev_polar)

fig3 = Figure(size=(600, 600))
ax3 = PolarAxis(fig3[1, 1],
    title="Elevation Plane (XZ, ϕ=0°)",
    theta_0=-π / 2,  # 90° at top
    direction=-1  # Clockwise to match elevation convention
)
lines!(ax3, elevation_angles_polar, GΩ_dB_elev_polar, linewidth=2, color=:blue)
save(joinpath(output_dir, "polar_elevation.png"), fig3)
# Polar plot in azimuth plane (XY-plane, elevation = 90° = horizon)
# Scan full circle: azimuth from 0° to 360°
azimuth_angles_polar = LinRange(0, 2π, 361)
elevation_polar_azim = π / 2  # 90° = horizon

k_xyz_azim_polar = [Kspace.k_xyz(elevation_polar_azim, ϕ) for ϕ in azimuth_angles_polar]
GΩ_azim_polar = [Kspace.gain(k, Ge, r_xyz, W) for k in k_xyz_azim_polar]
GΩ_dB_azim_polar = DspUtility.pow2db.(GΩ_azim_polar)

fig4 = Figure(size=(600, 600))
ax4 = PolarAxis(fig4[1, 1],
    title="Azimuth Plane (XY, θ=90°)"
)
lines!(ax4, azimuth_angles_polar, GΩ_dB_azim_polar, linewidth=2, color=:blue)
save(joinpath(output_dir, "polar_azimuth.png"), fig4)


## BEAM SCANNING ANIMATION WITH 3D VISUALIZATION

# Scan in elevation from 90° to 0° and azimuth from 0° to 360°
scan_elevations = LinRange(π / 2, 0, 31)
scan_azimuths = LinRange(0, π / 2, 31)

# Observation angles for 3D pattern
θ_3d = LinRange(0, π, 91)  # Elevation from 0 to 180° (full sphere)
φ_3d = LinRange(-π, π, 91)   # Azimuth full circle

# Observation angles (elevation) for polar plot
θ_obs = LinRange(0, 2π, 361)  # Full circle in elevation

# Calculate expected axis range based on number of elements
num_elements = length(element_positions)
array_gain_dB = 10 * log10(num_elements) + 11
max_expected_radius = array_gain_dB
axis_range = max_expected_radius

# Create figure for animation with 3 panels
fig_anim = Figure(size=(1800, 600))

# Left plot: 3D radiation pattern
ax_3d = Axis3(fig_anim[1, 1],
    title="3D Beam Pattern",
    aspect=(1, 1, 1),
    limits=(-axis_range, axis_range, -axis_range, axis_range, -axis_range, axis_range)
)
hidespines!(ax_3d)
hidedecorations!(ax_3d)

# Middle plot: Polar radiation pattern (elevation cut)
ax_pattern = PolarAxis(fig_anim[1, 2],
    title="Beam Pattern - Elevation Plane",
    theta_0=-π / 2,  #  0° at top.
    direction=-1  # Clockwise
)

# Right plot: Element phase weights
ax_phase = Axis(fig_anim[1, 3],
    xlabel="Element Index",
    ylabel="Phase [rad]",
    title="Element Phase Weights",
    limits=(0, 6, -π, π)
)

# Observables for animation
xs_3d_observable = Observable(zeros(length(θ_3d), length(φ_3d)))
ys_3d_observable = Observable(zeros(length(θ_3d), length(φ_3d)))
zs_3d_observable = Observable(zeros(length(θ_3d), length(φ_3d)))
color_3d_observable = Observable(zeros(length(θ_3d), length(φ_3d)))
gain_observable = Observable(zeros(length(θ_obs)))
phase_observable = Observable(zeros(5))

# Initial plots
surface!(ax_3d, xs_3d_observable, ys_3d_observable, zs_3d_observable,
    color=color_3d_observable, colormap=:jet1, shading=NoShading)
lines!(ax_pattern, θ_obs, gain_observable, linewidth=2, color=:blue)
scatter!(ax_phase, 1:5, phase_observable, markersize=15, color=:red)

# Animation - record to file
record(fig_anim, joinpath(output_dir, "beam_scan_3d.gif"), eachindex(scan_elevations); framerate=10) do frame_idx
    scan_elevation = scan_elevations[frame_idx]
    scan_azimuth = scan_azimuths[frame_idx]

    # Calculate phase weights for beam steering toward scan direction
    W_phases = [DspUtility.phase_weight(scan_elevation, scan_azimuth, element_positions[i])
                for i in 1:length(element_positions)]

    # Apply uniform amplitude weights with phase steering
    W_steering = ones(5, 1) .* exp.(im .* reshape(W_phases, 5, 1))

    # Calculate 3D gain pattern
    GΩ_3d = zeros(Float64, length(θ_3d), length(φ_3d))
    for (i, θ) in enumerate(θ_3d)
        for (j, φ) in enumerate(φ_3d)
            k_vec = Kspace.k_xyz(θ, φ)
            GΩ_3d[i, j] = Kspace.gain(k_vec, Ge, r_xyz, W_steering)
        end
    end

    # Convert to dB and clamp
    GΩ_dB_3d = DspUtility.pow2db.(GΩ_3d)
    min_value_dB = -10
    GΩ_dB_3d = clamp.(GΩ_dB_3d, min_value_dB, Inf)
    GΩ_dB_3d = GΩ_dB_3d .- min_value_dB  # Shift to positive values for radius

    # Convert to Cartesian coordinates (spherical -> Cartesian)
    xs_3d = [GΩ_dB_3d[i, j] * sin(θ_3d[i]) * cos(φ_3d[j]) for i in eachindex(θ_3d), j in eachindex(φ_3d)]
    ys_3d = [GΩ_dB_3d[i, j] * sin(θ_3d[i]) * sin(φ_3d[j]) for i in eachindex(θ_3d), j in eachindex(φ_3d)]
    zs_3d = [GΩ_dB_3d[i, j] * cos(θ_3d[i]) for i in eachindex(θ_3d), j in eachindex(φ_3d)]

    # Calculate gain pattern for polar plot (elevation cut at azimuth = 45°)
    GΩ_scan = zeros(Float64, length(θ_obs))
    for (idx, θ) in enumerate(θ_obs)
        k_vec = Kspace.k_xyz(θ, scan_azimuth)
        GΩ_scan[idx] = Kspace.gain(k_vec, Ge, r_xyz, W_steering)
    end

    # Convert to dB for polar plot
    GΩ_dB_scan = DspUtility.pow2db.(GΩ_scan)
    GΩ_dB_scan = clamp.(GΩ_dB_scan, -40, 200)

    # Update observables
    xs_3d_observable[] = xs_3d
    ys_3d_observable[] = ys_3d
    zs_3d_observable[] = zs_3d
    color_3d_observable[] = GΩ_dB_3d
    gain_observable[] = GΩ_dB_scan
    phase_observable[] = W_phases

    # Update titles
    ax_3d.title = "3D Pattern - Scan: θ=$(round(rad2deg(scan_elevation), digits=1))°, φ=$(round(rad2deg(scan_azimuth), digits=1))°"
    ax_pattern.title = "Elevation Cut - Scan: θ=$(round(rad2deg(scan_elevation), digits=1))°, φ=$(round(rad2deg(scan_azimuth), digits=1))°"
    ax_phase.title = "Phase Weights - Scan: θ=$(round(rad2deg(scan_elevation), digits=1))°, φ=$(round(rad2deg(scan_azimuth), digits=1))°"
end


