using Statistics

module DspUtility

export linear_array

"""
Convert between linear scale (magnitude) and decibel 
"""
function mag2db(magnitude)
	return 20 * log10.(magnitude)
end

"""
Convert between decibal scale (magnitude) and linear
"""
function db2mag(magDb)
	return 10^(magDb ./ 20)
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
function db2pow(powerDb)
	return 10^(powerDb ./ 10)
end

"""
Return the average signal power in dbW
"""
function power_dBW(signal)
	return pow2db(mean(abs(signal) .^ 2))
end

"""
Return the average signal power in dbm
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
Generate locations (1D) for an x*λ spaced array of N elements. N is even.
"""
linear_array(N, separation_lambda) = LinRange(-N / 2 + separation_lambda, N / 2 - separation_lambda, N)
end
