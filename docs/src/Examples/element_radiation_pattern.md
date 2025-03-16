# Default Element Radiation Pattern

``` @example elementRadiationPattern
using Plots
gr()

using ArrayRadiation


angleRad = LinRange(π / 2, -π / 2, 51);
angleDeg = rad2deg.(angleRad);

k_xyz = 2π*Kspace.elevation2k_hat.(angleRad)

element_gain_approximation = Kspace.cos_taper.(angleRad)
element_gain_approximation_dB = DspUtility.pow2db(element_gain_approximation)

y_lower_limit = -30 # dB
element_gain_approximation_dB = clamp.(element_gain_approximation_dB, y_lower_limit, Inf)


# Create the polar plot
plot(angleRad, element_gain_approximation_dB, 
    proj=:polar, 
    m=:none, lw=2,
    linecolor=:blue, 
    legend=false,
    grid=true,
    gridalpha=0.4,
    gridlinewidth=1,
    ylims=(y_lower_limit, 0.5)
)

title!("Element Radiation Pattern")
```
