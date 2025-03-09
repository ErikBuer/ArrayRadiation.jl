module Kspace
using ..DspUtility

c = 3e8 # Velocity of light lin vacuum.

"""
Calculate 1D cosine taper for angle θ [rad].
Optionally provide α. α=1.4 accounts for mutual coupling between elements.

- R. A. Dana, Electronically Scanned Arrays and K-Space Gain Formulation, Springer, 2019.
"""
cos_taper(θ, α = 1.4) = (cos(θ) + 0 * im)^α

"""
Calculate the maximum array element separation for maximum look angle.
Result is fraction of wavelenth.
"""
function maximum_element_separation(max_look_angle_deg)
	return 1 / (1 + sin(deg2rad(max_look_angle_deg)))
end


"""
Calculate `GΩ(k)`, the angular domain gain at angle k [rad] for an array with specified element weights.

- R. A. Dana, Electronically Scanned Arrays and K-Space Gain Formulation, Springer, 2019.

# Arguments
- `k`					Angle [rad].
- `element_gain`		The antena gain in direction `k` [dB].
- `element_placement_λ` The placement of each antenna element in wavelength units.
- `element_weighting`	The complex weight of each element.
- `α`					Cosine taper/mutual coupling coefficient α>1.
"""
function gain_2D(k::Real, element_gain::Real, element_placement_λ::AbstractVector, element_weighting::AbstractVector, α::AbstractFloat = 1.4)

	# Element gain at angle k
	Ge = element_gain

	# Account for forshortening and mutual coupling through (cosine tapering).
	Ge = (Ge + 0 * im) * cos_taper(k, α)

	r = element_placement_λ
	W = element_weighting

	# Array normal angle k [rad]
	k0 = 0

	Σ = sum
	numerator(m) = W[m] * exp(-im * (k - k0) * r[m])
	vecIndex = trunc.(Int64, LinRange(1, size(r)[1], size(r)[1]))
	return Ge * abs.(Σ(numerator, vecIndex))^2 / Σ(abs.(W) .^ 2)
end

end
