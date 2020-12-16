### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ 710bd51c-18d7-11eb-2a05-c7f9f8f34324
begin
	using HypergeometricFunctions
	using LinearAlgebra
	using SparseArrays
	using PyPlot
	using BenchmarkTools
end

# ╔═╡ 51590b9c-3c39-11eb-0a39-6bb9f81b7d68
md"#### Fractional Fourier-Kravchuk evolution

Date: December 11, 2020
"

# ╔═╡ 02f4dd56-18d9-11eb-2d38-379eca82ec3a
begin
	rc("text", usetex = true)
	ion()
end

# ╔═╡ be283892-18d0-11eb-0154-758e466c7a40
# Global constants are fine, I suppose
begin
	const p = 2.0 # For the symmetric Kravchuk polynomial
	const j = 16 # Angular momentum number
	const N = 2*j # Dimension
end;

# ╔═╡ ffe361c6-18d0-11eb-2aa6-af741f44c0f9
# There's a problem in the _₂F₁ who have a NaN value when n = q = N, the bottom right corner of the array
function k_p(n,q)
	if n == q == N # Ad hoc condition added for completeness of the calculation
		return 0
	else
		return HypergeometricFunctions._₂F₁general2(-n, -q, -N, p)
	end
end;

# ╔═╡ a4b42958-3c10-11eb-3efa-b505bbd7df59
begin
	A(n, j) = ((-1)^n / 2^j)
	B(n, j) = binomial(2 * j, n)
end;

# ╔═╡ d3a38eca-3c10-11eb-20a8-6f0267436c21
k_f(n,q) = A(n, j) * sqrt(B(n, j) * B(j + q, j)) * k_p(n, j + q);

# ╔═╡ d86e82a6-39e1-11eb-3e5d-05a1af7f904e
K = reshape([k_f(l,m) for m in -j:j for l in 0:N], (N+1,N+1)); # Kravchuk matrix

# ╔═╡ be9d003e-3c14-11eb-3f4b-4771e5c06bec
id(α,n) = exp(-1im * pi * α * n / 2.0); # Unit roots

# ╔═╡ c9a7e70c-3c17-11eb-3038-7327c598667a
idT = [[id(α,n) for n in 0:N] for α in range(0, stop = 2, length = 315)]; # Vectors for the unit roots of dimensions N×(partition of 2π)

# ╔═╡ c635d95a-3c18-11eb-37b1-5dcccc976e80
dia(α) = Diagonal(idT[α]); # Diagonal matrix of the unit roots vectors, are bounded for α ∈ [1, ..., length of the 2π partition]

# ╔═╡ d8c7b53a-3c19-11eb-2890-759ed15872ff
begin
	pint(q) = sparsevec(Dict(q + j + 1 => 1), N+1); # Vector of zeros with a unit entry at q ∈ (-j, ..., j)
	eigK(n) = K[n,:]; # Eigenfunctions of the Kravchuk oscillator
	disK(z) = [k_f(0,m+z) for m in -j:j]; # Displaced n = 0 eigenstates of the Kravchuk oscillator
end

# ╔═╡ aa913f54-3c34-11eb-267a-6dc645a56600
vK = transpose(K) * normalize(disK(8)); # Transformation of the initial condition to the Kravchuk space

# ╔═╡ 35f1846c-3c2d-11eb-20fa-7954ca2d82a9
vFrT(α) = K * dia(α) * vK; # Propagation function in Kravchuk space and returning to configuration space

# ╔═╡ 28a9a1de-3c20-11eb-09b4-61e0b5d3c9c5
FrK = [abs2.(vFrT(α)) for α in 1:1:length(idT)]; # Propagation

# ╔═╡ 3456ff94-3c35-11eb-0000-ef9e984797bb
# Sadly the resulting matrix is 315 × N, so we need to reaccomodate to be N × 315
begin
	FK = zeros(Float64,(N+1,315))
	for l in 1:N+1
		for m in 1:315
			FK[l,m] = FrK[m][l]
		end
	end
end

# ╔═╡ 91da32e4-3c35-11eb-11ec-3fd798a5d2b2
begin
	fig1 = figure(figsize = (10,5))
	ax1 = gca()
	ax1.imshow(FK, cmap = "hot")
	ax1.set_aspect("auto")
	
	tight_layout()
end

# ╔═╡ Cell order:
# ╟─51590b9c-3c39-11eb-0a39-6bb9f81b7d68
# ╠═710bd51c-18d7-11eb-2a05-c7f9f8f34324
# ╠═02f4dd56-18d9-11eb-2d38-379eca82ec3a
# ╠═be283892-18d0-11eb-0154-758e466c7a40
# ╠═ffe361c6-18d0-11eb-2aa6-af741f44c0f9
# ╠═a4b42958-3c10-11eb-3efa-b505bbd7df59
# ╠═d3a38eca-3c10-11eb-20a8-6f0267436c21
# ╠═d86e82a6-39e1-11eb-3e5d-05a1af7f904e
# ╠═be9d003e-3c14-11eb-3f4b-4771e5c06bec
# ╠═c9a7e70c-3c17-11eb-3038-7327c598667a
# ╠═c635d95a-3c18-11eb-37b1-5dcccc976e80
# ╠═d8c7b53a-3c19-11eb-2890-759ed15872ff
# ╠═aa913f54-3c34-11eb-267a-6dc645a56600
# ╠═35f1846c-3c2d-11eb-20fa-7954ca2d82a9
# ╠═28a9a1de-3c20-11eb-09b4-61e0b5d3c9c5
# ╠═3456ff94-3c35-11eb-0000-ef9e984797bb
# ╠═91da32e4-3c35-11eb-11ec-3fd798a5d2b2
