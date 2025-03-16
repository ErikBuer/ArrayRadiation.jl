# Copyright (c) 2020 - Norsk Datateknikk AS
using Plots
include("DspUtility.jl")

# This script calculates the gain of an array atenna with a point source antenna assumption.

# Antenna element placement vector
# Use vector for a linear array
# Place elements symmetrically around zero

# Generate locations (1D) for an λ/2 spaced array of N elements. N is even.
linearArray(N) = LinRange(-N/2+0.5, N/2-0.5, N)

r = linearArray(4)

# Antenna element weight
W = ones(size(r))
# Frequency of interrest
freq = 1e9
# Propagation speed
c = 3e8
# Wavelength
λ = c/freq
# Element gain at angle k
Ge_dB = 0   # isotropic
Ge = db2pow(Ge_dB)

# Array normal angle k [rad]
k0	= 0

# Array gain at angle k [rad]
function GΩ_dB(k, k0, Ge)
	Σ = sum
	i = im
	numerator(m)= W[m]*exp(-i*(k-k0)*r[m])
	vecIndex 	= trunc.(Int64, LinRange( 1, size(r)[1], size(r)[1]))
	return Ge*abs.( Σ( numerator, vecIndex ) )^2 / Σ(abs.(W).^2) 
end

# TODO account for forshortening
# TODO account for cosine tapering

gain(k) = GΩ_dB(k, k0, Ge)

angleRad 		= LinRange( -π, π, 1001)
angledeg		= rad2deg.(angleRad)
radPattern		= broadcast(gain, angleRad)
radPattern_dB	= broadcast( pow2db, radPattern )


plot(	angledeg, radPattern_dB,
			xlabel = "Angle [deg]",
			ylabel = "Gain [dB]",
			title  = "Radiation Pattern",
			ylims = (-30, maximum(radPattern_dB)+1	))