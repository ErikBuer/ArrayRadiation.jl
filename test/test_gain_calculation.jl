using ArrayRadiation
using Plots
gr()

element_separation_λ = 1/2;

i = im

element_count = 32;
# Place elements symmetrically around zero
r_λ = ArrayRadiation.linear_array(element_count, element_separation_λ)


# Elevation
angleRad = LinRange(π / 2, -π / 2, 501);
angledeg = rad2deg.(angleRad);

element_gain_approximation = Kspace.cos_taper.(angleRad)

k_xyz = 2π*Kspace.elevation2k_hat.(angleRad)
kx = getindex.(k_xyz, 1)

"""
Inter-element phase shift for scanning in direction θ [rad].

"""
α(θ, d_λ) = 2π*d_λ*sin(θ);

# Angles at which to scan the array.
scan_angles = LinRange(-π / 2, π / 2, 61);


plt1 = plot();
max_gain = 18;

anim = @animate for (index, scan_angle) in enumerate(scan_angles)
	global max_gain

	phase_increment = α(scan_angle, element_separation_λ)
	W_ang = LinRange(0, (element_count - 1) * phase_increment, element_count)

	# Antenna element weigth
	W = exp.(i .* W_ang)# Example of pointing to an angle.

	# Map angular domain gain calculation function.
	gain_Ω(k) = Kspace.gain_1D(k, 1, r_λ, W)

	radPattern = broadcast(gain_Ω, k_xyz).*element_gain_approximation
	radPattern_dB = broadcast(DspUtility.pow2db, abs.(radPattern))

	plot(angledeg, radPattern_dB,
		xlabel = "Angle [deg]",
		ylabel = "Gain [dB]",
		title  = "Radiation Pattern",
		ylims  = (-20, max_gain),
		reuse  = true,
		legend = false)
end

gif(anim, "plots/array_scan.gif", fps = 10)
