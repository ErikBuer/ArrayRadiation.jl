module AntennaElement

"""
    half_wave_dipole(θ::AbstractFloat, ϕ::AbstractFloat)

Calculate approximate radiation pattern for a half wave dipole.

## Arguments

- `θ`   Elevation [rad]
- `ϕ`   Azimuth [rad]

## References

- Costa, Abrao, Rego, A Method for Approximating Directivity Expressions in Generalized Radiation Patterns, Elsevier, 2025.

"""
function half_wave_dipole(θ::AbstractFloat, ϕ::AbstractFloat)
    # Variables for empirical modeling.
    c1 = 0.7796
    c2 = 0.2192
    return c1*sin(θ) + c2*sin(θ)^3
end

"""
    cardioid(θ::AbstractFloat, ϕ::AbstractFloat)

Calculate approximate radiation pattern for a cardioid antenna.

## Arguments

- `θ`   Elevation [rad]
- `ϕ`   Azimuth [rad]

## References

- Costa, Abrao, Rego, A Method for Approximating Directivity Expressions in Generalized Radiation Patterns, Elsevier, 2025.

"""
function cardioid(θ::AbstractFloat, ϕ::AbstractFloat)
    return 1 + cos(θ)
end

"""
    yagi_uda(θ::AbstractFloat, ϕ::AbstractFloat)

Calculate approximate radiation pattern for a Yagi Uda antenna.

The model is only valid for θ ∈ [−π/2, π/2].

The model is based on the following geometry:

- Three directors, spaced 0.20λ apart.
- Director lengths: [0.428λ, 0.424λ, 0.428λ]
- Reflector length of 0.482λ

## Arguments

- `θ`   Elevation [rad]
- `ϕ`   Azimuth [rad]

## References

- Costa, Abrao, Rego, A Method for Approximating Directivity Expressions in Generalized Radiation Patterns, Elsevier, 2025.

"""
function yagi_uda(θ::AbstractFloat, ϕ::AbstractFloat)
    if abs(θ)>π/2
        @warn "The model is only valid for θ ∈ [−π/2, π/2]. Current value: θ = $θ"
    end

    # Variables for empirical modeling.
    c1=  0.7535
    c2= -8.1865
    c3=  21.9141
    c4=  3.7784
    c5= -6.5321
    return c1 + c2*sin(θ)^2 + c3*sin(θ)^4 + c4*sin(θ)^6 + c5*sin(θ)^8
end

"""
    microstrip_patch(θ::AbstractFloat, ϕ::AbstractFloat)

Calculate approximate radiation pattern for a microstrip patch antenna.

The model is only valid for θ ∈ [−π/2, π/2].

## Arguments

- `θ`   Elevation [rad]
- `ϕ`   Azimuth [rad]

## References

- Costa, Abrao, Rego, A Method for Approximating Directivity Expressions in Generalized Radiation Patterns, Elsevier, 2025.

"""
function microstrip_patch(θ::AbstractFloat, ϕ::AbstractFloat)
    if abs(θ)>π/2
        @warn "The model is only valid for θ ∈ [−π/2, π/2]. Current value: θ = $θ"
    end

    # Variables for empirical modeling.
    c1=  0.9731
    c2= −0.8119
    c3=  5.7330
    return c1 + c2*sin(θ)^4 + c3*cos(θ)^4
end

end