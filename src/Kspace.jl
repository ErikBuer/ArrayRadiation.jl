module Kspace

using LinearAlgebra

"""
	k_xyz(θ::Real, ϕ::Real, λ0::Real)::Vector

Calculate k-space vector, representing antenna elevation and coordinates.
"""
k_xyz(θ::Real, ϕ::Real, λ0::Real)::Vector = 2*π/λ0 * [sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]

"""
    k2elevation(k_xyz::Vector)::Real

    Calculate the elevation angle θ from a k-space vector.
    The elevation angle is the angle between the z unit vector and the the k vector.
"""
function k2elevation(k_xyz::Vector)::Real
	kx, ky, kz = k_xyz
	return π/2-atan(kz, sqrt(kx^2 + ky^2))  # Elevation angle in radians.
end

"""
    elevation2k_hat(θ::Real)::Vector{Real}

    Convert an elevation angle θ (in radians) to a unit-length k-space vector in the XZ-plane.
    The resulting vector has the form [sin(θ), 0, cos(θ)].
"""
function elevation2k_hat(θ::Real)::Vector{Real}
    return [sin(θ), 0, cos(θ)]
end

"""
    k2azimuth(k_xyz::Vector)::Real

    Calculate the azimuth angle ϕ from a k-space vector.
    The azimuth angle is the angle between the projection of the vector on the XY-plane and the x-axis.
"""
function k2azimuth(k_xyz::Vector)::Real
    kx, ky, _ = k_xyz
    return atan(ky, kx)  # Azimuth angle in radians
end

"""
    azimuth2k_hat(ϕ::Real)::Vector{Real}

    Convert an azimuth angle ϕ (in radians) to a unit-length k-space vector in the XY-plane.
    The resulting vector has the form [cos(ϕ), sin(ϕ), 0].
"""
function azimuth2k_hat(ϕ::Real)::Vector{Real}
    return [cos(ϕ), sin(ϕ), 0]
end

"""
Calculate cosine taper for elevation angle θ [rad].
Optionally provide α. α=1.4 accounts for mutual coupling between elements.

- R. A. Dana, Electronically Scanned Arrays and K-Space Gain Formulation, Springer, 2019.
"""
cos_taper(θ::Real, α = 1.4) = cos(θ)^α

"""

Calculate cosine taper for normalized k-space vector [kx, ky, kz]/|k|.
Where |k| = 2π/λ0
Optionally provide α. α=1.4 accounts for mutual coupling between elements.

- R. A. Dana, Electronically Scanned Arrays and K-Space Gain Formulation, Springer, 2019.
"""
cos_taper_k_hat(k_hat_x::Real, k_hat_y::Real, α = 1.4) = (1-(k_hat_x[1]^2+k_hat_y[1]^2))^(α/2)

"""
Calculate the maximum array element separation for maximum look angle.
Result is fraction of wavelenth.
"""
function maximum_element_separation(max_look_angle_deg)
	return 1 / (1 + sin(deg2rad(max_look_angle_deg)))
end


"""
	gain_1D(k_xyz::AbstractVector{<:Real}, Ge::Real, r_λ::AbstractVector, element_weighting::AbstractVector)

Calculate `GΩ(k)`, the angular domain gain in direction k for an array with specified element weights.

# Arguments

- `k_xyz`			    k-space vector.
- `Ge`					The antena gain in direction `k_xyz` [dB].
- `r` 					The placement of each antenna element.
- `element_weighting`	The complex weight of each element.

## References

- R. A. Dana, Electronically Scanned Arrays and K-Space Gain Formulation, Springer, 2019.
"""
function gain_1D(k_xyz::AbstractVector{<:Real}, Ge::Real, r::AbstractVector, element_weighting::AbstractVector)
    M = length(r)
	r_xyz = [[r_i, 0, 0] for r_i in r]
	
	W = element_weighting

	# Array normal vector
	k0 = [ 0.0, 0.0, 2π ]

	Σ = sum
	numerator(m) = W[m] * exp(-im * dot(k_xyz - k0, r_xyz[m]))
	vecIndex = collect(Integer, 1:M)
	return Ge * abs.(Σ(numerator, vecIndex))^2 / Σ(abs.(W) .^ 2)
end

end
