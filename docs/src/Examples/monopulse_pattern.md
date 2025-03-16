# Monopulse Pattern

For monopulse radar processing, we use a sum beam and a difference beam (on reception) in each axis, to estimate an angle to a target with the goal of tracking it.

We can create these two beams by playing with the element weights.

``` @example WindowWeights
using Plots;
gr();

using ArrayRadiation

element_separation_λ = 1/2;
element_count = 32;
r = ArrayRadiation.linear_array(element_count, element_separation_λ)

# Antenna element weight
W = Window.taylor(32,4,-25)
W2 = Window.split_window(copy(W))

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

``` @example WindowWeights
sum(W)
```

``` @example WindowWeights
sum(W2)
```

Now lets compare their radiation pattern.

``` @example WindowWeights

elevationRad = LinRange(π / 2, -π / 2, 501);
k_xyz = 2π*Kspace.elevation2k_hat.(elevationRad)
elevationDeg = rad2deg.(elevationRad);


# Map K-space gain calculation function.
GΩ(k) = ArrayRadiation.Kspace.gain_1D(k, 1, r, W)
GΩ_lin = broadcast(GΩ, k_xyz)

# Calculate
GΩ_lin = Kspace.cos_taper.(elevationRad)


# Map K-space gain calculation function.
GΩ(k, _W) = Kspace.gain_1D(k, 1, r, _W)
GΩ1_lin = map(k -> GΩ(k, W), k_xyz)
GΩ2_lin = map(k -> GΩ(k, W2), k_xyz)

GΩ1_dB = DspUtility.pow2db.(abs.(GΩ1_lin))
GΩ2_dB = DspUtility.pow2db.(abs.(GΩ2_lin))

plot(elevationDeg, GΩ1_dB,
    xlabel = "Angle [deg]",
    ylabel = "GΩ [dB]",
    title  = "Array gain",
    ylims  = (-30, 18),
    reuse  = true,
    label  = "Sum Beam"
)

plot!(elevationDeg, GΩ2_dB, 
    label  = "Difference Beam"
)

```

As we can see, the difference beam is not great..
