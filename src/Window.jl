module Window

using ..DspUtility

"""
Returns `1:n` but omits m.
"""
function rangeOmit(n, m)
    return filter(x -> x != m, 1:n)
end

"""
    taylor(N::Integer, n_bar::Integer, sll_dB::Real)::AbstractVector

Create taylor window.


## Arguments

- `N`       The window length.
- `n_bar`   The number of nearly constant-level sidelobes adjacent to the main lobe.
- `sll_dB`  Peak sidelobe_level in dB.

## Example

```jldoctest
julia> using ArrayRadiation

julia> W = Window.taylor(64, 4, -35);

julia> round.( W[1:4], sigdigits=3 )
4-element Vector{Float64}:
 0.28
 0.288
 0.305
 0.329

julia> round( sum(W), sigdigits=5)
64.0
```
## References

- Spotlight Synthetic Aperture Radar: Signal Processing Algorithms, Artech House, 1995
"""
function taylor(N::Integer, n_bar::Integer, sll_dB::Real)::AbstractVector
    B = DspUtility.db2mag(abs(sll_dB))
    A = log(B + sqrt(B^2 - 1)) / π
    σ2 = n_bar^2 / (A^2 + (n_bar - 1/2)^2)

    i = collect(1 : (n_bar - 1));

    fNumerator(_m) = (-1)^(_m .+ 1) ./ 2 .* prod(1 .- (_m^2 / σ2) ./ (A^2 .+ (i .- 0.5).^2));
    fDenominator(_m) = prod(1 .- _m^2 ./ rangeOmit(n_bar - 1, _m).^2);
    
    m = i
    n = 0:(N - 1)

    _F(_m) = fNumerator(_m) ./ fDenominator(_m)
    F = _F.(m)
    
    win = 1 .+ 2 .* cos.(2π/N .* (n .- (N - 1) / 2) .* m') * F
    
    return win
end



"""
    split_taylor(N::Integer, n_bar::Integer, sll_dB::Real)::AbstractVector

Create split taylor window. For monopulse difference beam


## Arguments

- `N`       The window length.
- `n_bar`   The number of nearly constant-level sidelobes adjacent to the main lobe.
- `sll_dB`  Peak sidelobe_level in dB.

## Example

```jldoctest
julia> using ArrayRadiation

julia> N = 64;

julia> W = Window.split_taylor(N, 4, -35);

julia> round.( W[1:4], sigdigits=3 )
4-element Vector{Float64}:
 0.28
 0.288
 0.305
 0.329


julia> round.( W[N-3:N], sigdigits=3 )
4-element Vector{Float64}:
 -0.329
 -0.305
 -0.288
 -0.28

julia> round( sum(W), sigdigits=2)
-6.3e-15
```
"""
function split_taylor(N::Integer, n_bar::Integer, sll_dB::Real)::AbstractVector
    if mod(N,2) ==1
        @error "`split_taylor` only works with even number of elements!"
    end
    W = taylor(N, n_bar, sll_dB)

    # Invert last half of array.
    mid = N ÷ 2
    for i in 1:(mid)
        W[mid + i] = -W[mid + i]
    end
    
    return W
end

end