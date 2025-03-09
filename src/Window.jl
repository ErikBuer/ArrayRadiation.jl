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

Create taylor weighting (window).


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

![Taylor Window](plots/Taylor Window.png)

## References

- Carrar, Goodman and Majewski, Spotlight Synthetic Aperture Radar: Signal Processing Algorithms, Artech House, 1995
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
    split_window(W::AbstractVector)::AbstractVector

Split a window to create a difference beam

## Arguments

- `W`       The window to split.

## Example

```jldoctest
julia> using ArrayRadiation

julia> N = 64;

julia> W = Window.taylor(N, 4, -35);

julia> W = Window.split_window(W);

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
```

Below is an example of a split Taylor window.

![Split Taylor Window](plots/Split Taylor Window.png)

"""
function split_window(W::AbstractVector)::AbstractVector
    N = length(W)
   
    if mod(N,2) == 1
        @error "split_window only works with even number of elements!"
        return nothing
    end

    # Invert last half of array.
    mid = N ÷ 2
    for i in 1:(mid)
        W[mid + i] = -W[mid + i]
    end
    
    return W
end


"""
    cosine_q(M::Integer, q, scale = true)::AbstractVector

Create a coseine^q weighting (window).

## Arguments

- `M`       The length of the window.
- `q`       For `q=0` yields a uniform weight,`q=1` yelds a cosine weighting, and `q=2` yields a Hanning weighting.
- `scale`   Peak of window is 1, when not scaled. `sum( cosine_q(M, q) ) = M` when scaled.

## Example

```jldoctest
julia> using ArrayRadiation

julia> N = 64;

julia> W = Window.cosine_q(N, 1);

julia> round.( W[1:4], sigdigits=3 )
4-element Vector{Float64}:
 0.0385
 0.116
 0.192
 0.269

julia> round( sum(W), sigdigits=2)
64.0
```

![Cosine q Window](plots/Cosine q Window.png)

## References

- S. Yan, Broadband Array Processing, Springer, 2019

"""
function cosine_q(M::Integer, q, scale = true)::AbstractVector

    m = 1:M

    _W(_m) = cos(π * (_m-(M+1)/2) / (M) )^q

    W = _W.(m)

    if scale
        s = sum(W)
        W .*= M / s
    end
    
    return W
end

end
