module FrequencyDomain

"""
Calculate the x dB bandwidth in number of bins.

```jldoctest
julia> FrequencyDomain.x_dB_bandwidth([1,2,3,4,5,4,3,2,1], 3)
7
```

# Arguments
- `spectrum_dB` The signal in frequency domain.
- `x_dB` 		The level of reduction in the band, e.g. half powe BW (3 dB).
"""
function x_dB_bandwidth(spectrum_dB, x_dB)::Int
	max, argmax = findmax(spectrum_dB)
	spectrum_dB = spectrum_dB .- max
	spectrum_dB = x_dB .+ spectrum_dB
	upper_band  = spectrum_dB[argmax:end]
	lower_band  = spectrum_dB[begin:argmax]


	"""
	Travel down hill and return the index preciding a negative value.
	Assuming the first item to be positive.
	"""
	function bins_down_hill(sideband)
		for (index, value) in enumerate(sideband)
			if value < 0
				return index .- 1
			end
		end
	end

	upper_bins = bins_down_hill(upper_band)
	lower_bins = bins_down_hill(reverse(lower_band))

	return upper_bins + lower_bins - 1
end
end
