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
# Antenna element weigth
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


# Map K-space gain calculation function.
GΩ(k, _W) = Kspace.gain_2D(k, 1, r, _W)
GΩ1_lin = map(k -> GΩ(k, W), angleRad)
GΩ2_lin = map(k -> GΩ(k, W2), angleRad)

GΩ1_dB = DspUtility.pow2db.(abs.(GΩ1_lin))
GΩ2_dB = DspUtility.pow2db.(abs.(GΩ2_lin))

plot(angleDeg, GΩ1_dB,
    xlabel = "Angle [deg]",
    ylabel = "GΩ [dB]",
    title  = "Array gain",
    ylims  = (-30, 18),
    reuse  = true,
    label  = "Uniform weights"
)

plot!(angleDeg, GΩ2_dB, 
    label  = "Taylor weights"
    )

```

It is apparent that the Taylor weighted antenna has a much greater peak to sidelobe distance, comming with the penalty of a wider main beam.
