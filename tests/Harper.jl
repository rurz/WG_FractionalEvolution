### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ baf9406a-2aed-11eb-1f5b-af8d76740822
begin
	using SpecialMatrices
	using LinearAlgebra
	using SparseArrays
	using PyPlot
	using LaTeXStrings
	using BenchmarkTools
end

# ╔═╡ 9600ef3e-3d07-11eb-3ebf-810a3b29cfd7
md"#### Fractional Fourier-Harper transform

Date: December 12, 2020
"

# ╔═╡ 47f9b180-2af2-11eb-2f04-9b70cbe71785
begin
	rc("text", usetex = true)
	ion()
end

# ╔═╡ 9d898a86-3cf8-11eb-2baa-e1bc43be6816
# Global constants: j = angular momentum; N = total angular numbers; dim = space dimension
begin
	const j = 16
	const N = 2 * j
	const dim = N + 1
end;

# ╔═╡ 649a4a96-2aef-11eb-04cb-c5bba93cebb1
function Δ1()
	Δ1 = zeros(Float64, dim)
	Δ1[1] = -2
	Δ1[2] = 1
	Δ1[end] = 1
	Δ1 = SpecialMatrices.Circulant(Δ1)
end;

# ╔═╡ 6d3bff0a-2af1-11eb-0863-81b26e690aab
function Δ2()
	Δ2 = zeros(Float64,dim)
	for i in -j:j
		Δ2[i + 1 + j] = -4 * (sin(pi * i / N))^2
	end
	return Δ2 = LinearAlgebra.diagm(Δ2)
end

# ╔═╡ da00b368-2af1-11eb-3111-a1bde1b4f2dd
# Natural difference analog of the discrete harmonic oscillator

H = -(1/2.0)*(Δ1() + Δ2()); 

# ╔═╡ 0e0bcffa-2af2-11eb-173b-b50a251e6015
# Eigenvectors h(m,k) for -j ≤ m ≤ j and 0 ≤ k ≤ N

h = LinearAlgebra.eigvecs(H);

# ╔═╡ 093f752a-3ce6-11eb-37a8-9502448cb5f5
# Set of Eigenvalues {k}

k = LinearAlgebra.eigvals(H);

# ╔═╡ 80445288-3e76-11eb-0e96-0bef8c278a79
η= [l for l in 0:N]; # Imported symmetry. Real eigenvalues replaced by linear number spectrum

# ╔═╡ 1fe9fce6-3ce6-11eb-2073-ad386eec723e
# Green matrix propagator 

G(m,μ,t) = sum([h[m + j + 1, l + 1] * exp(-1im * pi * t * η[l + 1]/2.0) * h[μ + j + 1, l + 1] for l in 0:N]);

# ╔═╡ 9632d7ea-3cfa-11eb-26c1-1910c42de58b
P(t) = reshape([G(m,μ,t) for m in -j:j for μ in -j:j], (dim, dim))

# ╔═╡ 051a0a20-3cfb-11eb-2f0c-51d2a26c6158
function P2(t)
	p2 = zeros(Float64, (dim, dim))
	for μ in -j:1:j
		for m in -j:1:j
			p2[m + j + 1, μ + j + 1] = G(m, μ, t)
		end
	end
	return p2
end

# ╔═╡ b2e3b8e6-3cff-11eb-190b-b96e32278be5
begin
	pint(q) = sparsevec(Dict(q + j + 1 => 1), dim); # Vector of zeros with a unit entry at q ∈ (-j, ..., j)
	eigh(n) = h[:,n]
end;

# ╔═╡ 4446e684-3f5b-11eb-2ed3-0dd74b1e2ec9
function hdc(d)
	ne = zeros(Float64,dim)
	for l in 1:dim-d
		ne[l] = h[l+d,1]
	end
	return ne
end

# ╔═╡ 52a41998-3d00-11eb-18fb-4f3336cc36f2
α = range(0, stop = 2, length = 315);

# ╔═╡ e6398104-3d03-11eb-0030-a3580b5fcf94
evol(θ) = abs2.(P(α[θ]) * normalize(hdc(8)));

# ╔═╡ bacd0730-3d04-11eb-394f-f12bfbb73b55
evolT = [evol(l) for l in 1:315];

# ╔═╡ 72edde46-3d05-11eb-0d98-472ce0d9a01f
begin
	Ankara = zeros(Float64, (dim, 315))
	for l in 1:dim
		for m in 1:315
			Ankara[l, m] = evolT[m][l]
		end
	end
end

# ╔═╡ ca90feee-3d05-11eb-352b-bd6efaf7d6fb
begin
	fig1 = figure(figsize = (10,5))
	ax1 = gca()
	ax1.imshow(Ankara, cmap = "hot")
	ax1.set_aspect("auto")
	
	tight_layout()
end

# ╔═╡ Cell order:
# ╟─9600ef3e-3d07-11eb-3ebf-810a3b29cfd7
# ╠═baf9406a-2aed-11eb-1f5b-af8d76740822
# ╠═47f9b180-2af2-11eb-2f04-9b70cbe71785
# ╠═9d898a86-3cf8-11eb-2baa-e1bc43be6816
# ╠═649a4a96-2aef-11eb-04cb-c5bba93cebb1
# ╠═6d3bff0a-2af1-11eb-0863-81b26e690aab
# ╠═da00b368-2af1-11eb-3111-a1bde1b4f2dd
# ╠═0e0bcffa-2af2-11eb-173b-b50a251e6015
# ╠═093f752a-3ce6-11eb-37a8-9502448cb5f5
# ╠═80445288-3e76-11eb-0e96-0bef8c278a79
# ╠═1fe9fce6-3ce6-11eb-2073-ad386eec723e
# ╠═9632d7ea-3cfa-11eb-26c1-1910c42de58b
# ╠═051a0a20-3cfb-11eb-2f0c-51d2a26c6158
# ╠═b2e3b8e6-3cff-11eb-190b-b96e32278be5
# ╠═4446e684-3f5b-11eb-2ed3-0dd74b1e2ec9
# ╠═52a41998-3d00-11eb-18fb-4f3336cc36f2
# ╠═e6398104-3d03-11eb-0030-a3580b5fcf94
# ╠═bacd0730-3d04-11eb-394f-f12bfbb73b55
# ╠═72edde46-3d05-11eb-0d98-472ce0d9a01f
# ╠═ca90feee-3d05-11eb-352b-bd6efaf7d6fb
