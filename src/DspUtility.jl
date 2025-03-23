module DspUtility

export linear_array, antenna_matrix


"""
	_mean(x)::Real

Calculate the mean of a vector. 

```jldoctest
julia> using ArrayRadiation

julia> a = [1,2,3,4];

julia> ArrayRadiation.DspUtility._mean(a)
2.5
```
"""
_mean(x)::Real = sum(x)/length(x)


"""
Convert between linear scale (magnitude) and decibel 
"""
function mag2db(magnitude)
	return 20 * log10.(magnitude)
end

"""
Convert between decibal scale (magnitude) and linear
"""
function db2mag(mag_dB)
	return 10^(mag_dB ./ 20)
end

"""
Convert between linear scale (magnitude) and decibel 
"""
function pow2db(power)
	return 10 * log10.(power)
end

"""
Convert between decibal scale (power) and linear
"""
function db2pow(power_dB)
	return 10^(power_dB ./ 10)
end

"""
Return the average signal power in dBW
"""
function power_dBW(signal)
	return pow2db(_mean(abs(signal) .^ 2))
end

"""
Return the average signal power in dBm
"""
function power_dBm(signal)
	return power_dBW(signal) + 30
end

"""
Return the signal energy in Joule
"""
function energy(signal)
	return sum(abs(signal) .^ 2)
end

"""
Generate locations (1D) for a evenly spaced array of N elements. N is even.

## Arguments
- `N`          	Number of elements (assumed to be even).
- `separation` 	Distance between adjacent elements.
"""
linear_array(N, separation) = LinRange(- (N - 1) / 2 * separation, (N - 1) / 2 * separation, N)

"""
    antenna_matrix(M, N, separation)

Generate locations (2D) for an evenly spaced MxN matrix of antenna elements, 
ensuring symmetry around the origin.

All elements are placed in the XY-plane.

## Arguments
- `M`          : Number of elements along the Y-axis.
- `N`          : Number of elements along the X-axis.
- `separation` : Distance between adjacent elements.

## Example

```jldoctest
julia> using ArrayRadiation

julia> DspUtility.antenna_matrix(1, 2, 0.5)
1×2 Matrix{Tuple{Float64, Float64, Float64}}:
 (-0.25, 0.0, 0.0)  (0.25, 0.0, 0.0)

julia> DspUtility.antenna_matrix(1, 3, 1/2)
1×3 Matrix{Tuple{Float64, Float64, Float64}}:
 (-0.5, 0.0, 0.0)  (0.0, 0.0, 0.0)  (0.5, 0.0, 0.0)

julia> DspUtility.antenna_matrix(2, 1, 1/2)
2×1 Matrix{Tuple{Float64, Float64, Float64}}:
 (0.0, -0.25, 0.0)
 (0.0, 0.25, 0.0)
```
"""
function antenna_matrix(M, N, separation)
    x_positions = LinRange(-(N - 1) / 2 * separation, (N - 1) / 2 * separation, N)
    y_positions = LinRange(-(M - 1) / 2 * separation, (M - 1) / 2 * separation, M)

    # Create a matrix of 3D coordinate vectors
    r_xyz = reshape([(x, y, 0.0) for y in y_positions for x in x_positions], M, N)

    return r_xyz
end

"""
    discard_low_values(scalar_value::Real, lower_limit::Real)

This function takes a scalar value and compares it against a specified lower limit.
If the scalar value is below the `lower_limit`, the function returns `nothing`.
Otherwise, it returns the scalar value unchanged.

# Arguments
- `scalar_value::Real`: The scalar value to compare.
- `lower_limit::Real`: The threshold value below which the scalar is discarded.

## Example
``` jldoctest
julia> using ArrayRadiation

julia> ArrayRadiation.DspUtility.discard_low_values(50, 40)
50

julia> ArrayRadiation.DspUtility.discard_low_values(30, 40)

```
"""
function discard_low_values(scalar_value::Real, lower_limit::Real)
    if scalar_value < lower_limit
        return nothing
    else
        return scalar_value
    end
end


end
