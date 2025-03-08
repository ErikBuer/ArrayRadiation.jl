# Beam Steering

``` @example BeamSteering
using Plots;
gr();

using ArrayRadiation

element_separation_λ = 1/2;

# Place elements symmetrically around zero

element_count = 32;

r = ArrayRadiation.linear_array(element_count, element_separation_λ)

# Frequency of interrest
freq = 12e9
# Propagation speed
c = 3e8
# Wavelength
λ = c / freq

angleRad = LinRange(π / 2, -π / 2, 501);
angleDeg = rad2deg.(angleRad);

"""
Inter-element phase shift for scanning in direction k [rad].
"""
α(θ) = -π * element_separation_λ * cos(π / 2 - θ);

# Angles at which to scan the array.
scan_angles = LinRange(-π / 2, π / 2, 81);


plt1 = plot();

anim = @animate for (index, scan_angle) in enumerate(scan_angles)

	phase_increment = α(scan_angle)
	W_ang = LinRange(0, (element_count - 1) * phase_increment, element_count)

	# Antenna element weigth
	W = exp.(im .* W_ang) # Example of pointing to an angle.

	# Map K-space gain calculation function.
	GΩ(k) = ArrayRadiation.Kspace.gain_2D(k, 0, r, W)

	GΩ_lin = broadcast(GΩ, angleRad)
	GΩ_dB = broadcast(ArrayRadiation.DspUtility.pow2db, abs.(GΩ_lin))

	plot(angleDeg, GΩ_dB,
		xlabel = "Angle [deg]",
		ylabel = "GΩ [dB]",
		title  = "Array Gain",
		ylims  = (-30, 20),
		reuse  = true,
		legend = false
		)
end

gif(anim, "array_scan.gif", fps = 10)

```
