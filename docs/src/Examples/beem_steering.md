# Beam Steering

``` @example BeamSteering
using Plots;
gr();

using ArrayRadiation


# Place elements symmetrically around zero
element_separation_λ = 1/2;
element_count = 64;
r = ArrayRadiation.linear_array(element_count, element_separation_λ)


angleRad = LinRange(π / 2, -π / 2, 501);
angleDeg = rad2deg.(angleRad);

k_xyz = 2π*Kspace.elevation2k_hat.(angleRad)

element_gain_approximation = Kspace.cos_taper.(angleRad)

"""
Inter-element phase shift for scanning in direction θ [rad].
"""
α(θ, d_λ) = 2π*d_λ*sin(θ);

# Angles at which to scan the array.
scan_angles = LinRange(-π / 2, π / 2, 81);

# Use Hanning element weights. 
Wp = Window.cosine_q(element_count, 2)

plt1 = plot();

anim = @animate for (index, scan_angle) in enumerate(scan_angles)

	phase_increment = α(scan_angle, 1/2)
	W_ang = LinRange(0, (element_count - 1) * phase_increment, element_count)

	# Antenna element complex weight
	W = Wp.*exp.(im .* W_ang) # Example of pointing to an angle.

	# Map K-space gain calculation function.
	GΩ(k) = ArrayRadiation.Kspace.gain_1D(k, 1, r, W)

	GΩ_lin = broadcast(GΩ, k_xyz).*element_gain_approximation
	GΩ_dB = broadcast(ArrayRadiation.DspUtility.pow2db, abs.(GΩ_lin))

	plot(angleDeg, GΩ_dB,
		xlabel = "Angle [deg]",
		ylabel = "GΩ [dB]",
		title  = "Array Gain",
		ylims  = (-80, 20),
		reuse  = true,
		legend = false
		)
end

gif(anim, "array_scan.gif", fps = 10)

```
