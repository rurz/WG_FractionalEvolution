### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ 0b0902da-40f0-11eb-3e10-a1fe1bd1f364
begin
	using PyPlot
	using SpecialFunctions
	using HypergeometricFunctions
end

# ╔═╡ 0b507856-40f1-11eb-244f-f15c1038b96f
begin
	rc("text", usetex = true)
	ion()
end

# ╔═╡ eaadeff2-546b-11eb-1289-9bade75bad45
function J(α, β, n, z)
    T = eltype(z/2)
    tot = zero(T)
    for k in 0:n
        tot += binomial(α + n, k) * binomial(β + n, n - k) * ((z + 1)/2)^k * ((z - 1)/2)^(n - k)
    end
    tot
end;

# ╔═╡ Cell order:
# ╠═0b0902da-40f0-11eb-3e10-a1fe1bd1f364
# ╠═0b507856-40f1-11eb-244f-f15c1038b96f
# ╠═eaadeff2-546b-11eb-1289-9bade75bad45
