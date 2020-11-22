# Copyright (c) 2020 - Norsk Datateknikk AS
using PyPlot
include("DspUtility.jl")

# This script calculates the gain of an array atenna with a point source antenna assumption.

# Antenna element placement vector
# Use vector for a linear array
# Place elements symmetrically around zero
r = LinRange(-64,64,129)*0.5
# Antenna element weigth
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

# Array gain at angle k [rad]
k0	= 0

function GΩ_dB(k, k0, Ge)
	Σ = sum
	i = im
	numerator(m)= W[m]*exp(-i*(k-k0)*r[m])
	vecIndex 	= trunc.(Int64, LinRange( 1, size(r)[1], size(r)[1]))
	return Ge*broadcast( abs, Σ( numerator, vecIndex ) )^2 / Σ(Σ(broadcast(abs, W).^2)) 
end

gain(k) = GΩ_dB(k, k0, Ge)

angleRad 		= LinRange( -π, π, 513)
angledeg			= rad2deg.(angleRad)
radPattern		= broadcast(gain, angleRad)
radPattern_dB	= broadcast( pow2db, radPattern )


plot(	angledeg, radPattern_dB,
			xlabel = "Angle [deg]",
			ylabel = "Gain [dB]",
			title  = "Radiation Pattern",
			ylims = (-30, maximum(radPattern_dB)+1	))