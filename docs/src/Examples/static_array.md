# Static Array

Lets create an array and look at its radiation pattern.

First we must place the antenna elements. Lets give them λ/2 spacing and spread them linearly.

``` @example StaticArray
using Plots;
gr();

using ArrayRadiation

element_separation_λ = 1/2;

# Place elements symmetrically around zero
element_count = 32;

r = ArrayRadiation.linear_array(element_count, element_separation_λ)

scatter(r, zeros(length(r)), 
    markershape=:circle, 
    # markercolor=:blue, 
    markersize=4, 
    xlabel="Element Position [λ]", 
    title="Antenna Element Positions", 
    legend=false, 
    #grid=false,
    yticks= false,
    ylims  = (-1, 1),
)
```

These elements now represent antennas. For the time being, these will be isotropic radiators.
We can apply different weights to each element so lets look at how this will affect our array.

Lets start by giving each element a uniform weight:

``` @example StaticArray
# Antenna element weigth
W = ones(element_count)

scatter(r, W, 
    marker=:circle, 
    linecolor=:blue, 
    markersize=4, 
    xlabel="Element Position [λ]", 
    ylabel="Weight", 
    title="Antenna Element Weights", 
    legend=false, 
    grid=true,
    ylims  = (0, 1.1),
)
```

With this defined, we can calculate the radiation pattern of the array.

``` @example StaticArray

angleRad = LinRange(π / 2, -π / 2, 501);
angleDeg = rad2deg.(angleRad);


# Map K-space gain calculation function.
GΩ(k) = ArrayRadiation.Kspace.gain_2D(k, 0, r, W)
GΩ_lin = broadcast(GΩ, angleRad)

GΩ_dB = ArrayRadiation.DspUtility.pow2db.(abs.(GΩ_lin))

plot(angleDeg, GΩ_dB,
    xlabel = "Angle [deg]",
    ylabel = "GΩ [dB]",
    title  = "Array gain",
    ylims  = (-30, 18),
    reuse  = true,
    legend = false
)

```
