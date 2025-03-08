using Plots;
gr();

using ArrayRadiation

# element_separation_lambda =  ArrayTools.maximum_element_separation( 50 )
element_separation_lambda = 1/2;
@info "Antenna element separation is $element_separation_lambda"

i = im

# Place elements symmetrically around zero

element_count = 24;
# TODO subarray_size = 2;

r = ArrayRadiation.linear_array(element_count, element_separation_lambda)

# Frequency of interrest
freq = 24e9
# Propagation speed
c = 3e8
# Wavelength
λ = c / freq

# TODO add input for number of parallel elements. Control it through weighting. Control it through weighting.


angleRad = LinRange(π / 2, -π / 2, 501);
angledeg = rad2deg.(angleRad);

#w_magnitude = Windows.bartlett_hann(element_count);

"""
Inter-element phase shift for scanning in direction k [rad].
Assumes λ/2 spacing.
"""
α(θ) = -π * element_separation_lambda * cos(π / 2 - θ);

# Angles at which to scan the array.
scan_angles = LinRange(-π / 4, π / 4, 61);


plt1 = plot();
max_gain = 15;

anim = @animate for (index, scan_angle) in enumerate(scan_angles)
	global max_gain

	phase_increment = α(scan_angle)
	W_ang = LinRange(0, (element_count - 1) * phase_increment, element_count)
	#@info "W_ang = %s. \n" round.( W_ang; digits=2 )

	# Antenna element weigth
	W = exp.(i .* W_ang)# Example of pointing to an angle.
	#	W	= W.*w_magnitude;

	# Account for subarrays with uniform weighting.
	# W = ArrayRadiation.subarray_weighting_2d( W, subarray_size );

	# Map K-space gain calculation function.
	gain(k) = ArrayRadiation.Kspace.k_space_gain_2D(k, 0, r, W)

	radPattern = broadcast(gain, angleRad)
	radPattern_dB = broadcast(DspUtility.pow2db, abs.(radPattern))

	#half_power_beamwidth_deg = FrequencyDomain.x_db_bandwidth( radPattern_dB, 3 )*2*π/1001;
	#@info "The RX half-power (3 dB) beamwidth is %s degrees. \n" round( half_power_beamwidth_deg; digits = 2 )

	if max_gain < maximum(radPattern_dB) + 1
		max_gain = maximum(radPattern_dB) + 1
	end


	plot(angledeg, radPattern_dB,
		xlabel = "Angle [deg]",
		ylabel = "Gain [dB]",
		title  = "Radiation Pattern",
		ylims  = (-20, max_gain),
		reuse  = true,
		legend = false)
end

gif(anim, "array_scan.gif", fps = 4)
