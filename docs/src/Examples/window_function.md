# Window Function

``` @example WindowWeights
using Plots;
gr();

using ArrayRadiation

element_separation_λ = 1/2;
element_count = 32;
r = ArrayRadiation.linear_array(element_count, element_separation_λ)
```

These elements now represent antennas. For the time being, these will be isotropic radiators.
We can apply different weights to each element so lets look at how this will affect our array.

Lets compare a linear weight to a commonly used Taylor weighting scheme.

``` @example WindowWeights
# Antenna element weight
W = ones(element_count)
W2 = Window.taylor(32,4,-25)

scatter(r, W, 
    marker=:circle, 
    linecolor=:blue, 
    markersize=4, 
    xlabel="Element Position [λ]", 
    ylabel="Weight", 
    title="Antenna Element Weights", 
    legend=false, 
    grid=true,
)
scatter!(r, W2)
```

As we can see, both weight functions have the same sum.

``` @example WindowWeights
sum(W)
```

``` @example WindowWeights
sum(W2)
```

Now lets compare their radiation pattern.

``` @example WindowWeights

angleRad = LinRange(π / 2, -π / 2, 501);
angleDeg = rad2deg.(angleRad);

k_xyz = 2π*Kspace.elevation2k_hat.(angleRad)

element_gain_approximation = AntennaElement.cos_taper.(angleRad)


# Map K-space gain calculation function.
GΩ(k, _W) = Kspace.gain_1D(k, 1, r, _W)
GΩ1_lin = map(k -> GΩ(k, W), k_xyz).*element_gain_approximation
GΩ2_lin = map(k -> GΩ(k, W2), k_xyz).*element_gain_approximation

GΩ1_dB = DspUtility.pow2db.(abs.(GΩ1_lin))
GΩ2_dB = DspUtility.pow2db.(abs.(GΩ2_lin))

plot(angleDeg, GΩ1_dB,
    xlabel = "Angle [deg]",
    ylabel = "GΩ [dB]",
    title  = "Array gain",
    ylims  = (-30, 18),
    reuse  = true,
    label  = "Uniformly weighted"
)

plot!(angleDeg, GΩ2_dB, 
    label  = "Taylor weighted"
)
```

It is apparent that the Taylor weighted antenna has greater peak to sidelobe distance, comming with the penalty of a wider main beam.

## Window Comparisons

Weighting distributions (windows) have been studied extensively over the years.
The reader is encouraged to read up on the different windows and their properties.

``` @example WindowWeights

W3 = Window.cosine_q(element_count, 1)
W4 = Window.cosine_q(element_count, 2)


scatter(r, W2, 
    marker=:circle, 
    linecolor=:blue, 
    markersize=4, 
    xlabel="Element Position [λ]", 
    ylabel="Weight", 
    title="Antenna Element Weights", 
    legend=true, 
    grid=true,
    label  = "Taylor weights",
)
scatter!(r, W3, label  = "Cosine weights")
scatter!(r, W4, label  = "Hanning weights")
```

``` @example WindowWeights
GΩ3_lin = map(k -> GΩ(k, W3), k_xyz).*element_gain_approximation
GΩ4_lin = map(k -> GΩ(k, W4), k_xyz).*element_gain_approximation

GΩ3_dB = DspUtility.pow2db.(abs.(GΩ3_lin))
GΩ4_dB = DspUtility.pow2db.(abs.(GΩ4_lin))

plot(angleDeg, GΩ2_dB,
    xlabel = "Angle [deg]",
    ylabel = "GΩ [dB]",
    title  = "Array gain",
    ylims  = (-50, 18),
    reuse  = true,
    label  = "Taylor weighted"
)

plot!(angleDeg, GΩ3_dB, 
    label  = "Cosine weighted"
)

plot!(angleDeg, GΩ4_dB, 
    label  = "Hanning weighted"
)
```
