module Kspace

using LinearAlgebra

"""
	k_xyz(θ::Real, ϕ::Real, λ0::Real=1)::Vector

Calculate k-space vector, representing antenna elevation and coordinates.
"""
k_xyz(θ::Real, ϕ::Real, λ0::Real=1)::Vector = 2 * π / λ0 * [sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ)]

"""
    k2elevation(k_xyz::Vector)::Real

Calculate the elevation angle θ from a k-space vector.

The elevation angle is the angle between the z unit vector and the the k vector.
"""
function k2elevation(k_xyz::Vector)::Real
    kx, ky, kz = k_xyz
    return π / 2 - atan(kz, sqrt(kx^2 + ky^2))  # Elevation angle in radians.
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
    cos_taper_k_hat(k_hat_x::Real, k_hat_y::Real, α = 1.4)

Calculate cosine taper for normalized k-space vector [kx, ky, kz]/|k|.

Where |k| = 2π/λ0

Optionally provide α. α=1.4 accounts for mutual coupling between elements.

- R. A. Dana, Electronically Scanned Arrays and K-Space Gain Formulation, Springer, 2019.
"""
cos_taper_k_hat(k_hat_x::Real, k_hat_y::Real, α=1.4) = (1 - (k_hat_x[1]^2 + k_hat_y[1]^2))^(α / 2)


"""
    maximum_element_separation(max_look_angle_deg)

Calculate the maximum array element separation for maximum look angle.

Result is fraction of wavelenth.
"""
function maximum_element_separation(max_look_angle_deg)
    return 1 / (1 + sin(deg2rad(max_look_angle_deg)))
end


"""
	gain_1D(k_xyz::AbstractVector{<:Real}, Ge::Real, r_x::AbstractVector, element_weights::AbstractVector)::Real

Calculate GΩ(k), the angular domain gain in direction k for an array with specified element weights.

# Arguments

- `k_xyz`			    k-space vector.
- `Ge`		            The antenna gain in direction `k_xyz`.
- `r_x` 				The placement of each antenna element. All elements are placed along the x-axis.
- `element_weights`	    The complex weight of each element.

## References

- R. A. Dana, Electronically Scanned Arrays and K-Space Gain Formulation, Springer, 2019.
"""
function gain_1D(k_xyz::AbstractVector{<:Real}, Ge::Real, r_x::AbstractVector, element_weights::AbstractVector)::Real
    if length(r_x) != length(element_weights)
        throw(ArgumentError("element_weights must have the same length as r_x. Got length $(length(element_weights)) for element_weights and length $(length(r_x)) for r_x."))
    end

    M = length(r_x)
    r_xyz = [[r_i, 0, 0] for r_i in r_x]

    W = element_weights

    # Array normal vector
    k0 = [0.0, 0.0, 2π]

    Σ = sum
    numerator(m) = W[m] * exp(-im * dot(k_xyz - k0, r_xyz[m]))
    vecIndex = collect(Integer, 1:M)
    return Ge * abs.(Σ(numerator, vecIndex))^2 / Σ(abs.(W) .^ 2)
end

"""
	gain(k_xyz::AbstractVector{<:Real}, Ge::Real, r_xyz::AbstractArray, element_weights::AbstractArray)

Calculate GΩ(k), the angular domain gain in direction k for an array with specified element weights.

# Arguments

- `k_xyz`			    k-space vector.
- `Ge`		            The antenna gain in direction `k_xyz`.
- `r_xyz` 				The placement of each antenna element in cartesian coordinates.
- `element_weights`	    The complex weight of each element. Must have the same size/shape as `r_xyz`.

## References

- R. A. Dana, Electronically Scanned Arrays and K-Space Gain Formulation, Springer, 2019.
"""
function gain(k_xyz::AbstractVector{<:Real}, Ge::Real, r_xyz::AbstractArray, element_weights::AbstractArray)::Real
    if size(r_xyz) != size(element_weights)
        throw(ArgumentError("element_weights must have the same size/shape as r_xyz. Got size $(size(element_weights)) for element_weights and size $(size(r_xyz)) for r_xyz."))
    end

    _r_xyz = reshape(r_xyz, :)
    W = reshape(element_weights, :)

    M = length(_r_xyz)

    # Array normal vector
    k0 = [0.0, 0.0, 2π]

    Σ = sum
    numerator(m) = W[m] * exp(-im * dot(k_xyz - k0, _r_xyz[m]))
    vecIndex = collect(Integer, 1:M)
    return Ge * abs.(Σ(numerator, vecIndex))^2 / Σ(abs.(W) .^ 2)
end

"""
	gain(k_xyz::AbstractVector{<:Real}, Ge::Real, r_xyz::AbstractArray, element_weights::AbstractArray)

Calculate GΩ(k), the angular domain gain in direction k for an array with specified element weights.

# Arguments

- `k_xyz`               k-space vector.
- `Ge(k_xyz, element_index)`   A function returning the antenna gain of each element in direction `k_xyz`. `Ge(k_xyz, element_index)` should return a `Real` scalar value.
- `r_xyz`               The placement of each antenna element in cartesian coordinates.
- `element_weights`     The complex weight of each element. Must have the same size/shape as `r_xyz`.

## References

- R. A. Dana, Electronically Scanned Arrays and K-Space Gain Formulation, Springer, 2019.
"""
function gain(k_xyz::AbstractVector{<:Real}, Ge::Function, r_xyz::AbstractArray, element_weights::AbstractArray)::Real
    if size(r_xyz) != size(element_weights)
        throw(ArgumentError("element_weights must have the same size/shape as r_xyz. Got size $(size(element_weights)) for element_weights and size $(size(r_xyz)) for r_xyz."))
    end

    # Create _Ge_matrix with same shape as r_xyz
    # Infer the type from a test call
    test_idx = first(eachindex(r_xyz))
    T = typeof(Ge(k_xyz, test_idx))
    _Ge_matrix = Array{T}(undef, size(r_xyz))

    # Calculate gain for each element
    for i in eachindex(r_xyz)
        _Ge_matrix[i] = Ge(k_xyz, i)
    end

    _Ge = reshape(_Ge_matrix, :)
    _r_xyz = reshape(r_xyz, :)
    W = reshape(element_weights, :)

    M = length(_r_xyz)

    # Array normal vector
    k0 = [0.0, 0.0, 2π]

    Σ = sum
    numerator(m) = sqrt(_Ge[m]) * W[m] * exp(-im * dot(k_xyz - k0, _r_xyz[m]))
    vecIndex = collect(Integer, 1:M)
    return abs.(Σ(numerator, vecIndex))^2 / Σ(abs.(W) .^ 2)
end

end
