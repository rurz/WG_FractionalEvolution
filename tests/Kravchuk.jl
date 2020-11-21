### A Pluto.jl notebook ###
# v0.12.11

using Markdown
using InteractiveUtils

# ╔═╡ 710bd51c-18d7-11eb-2a05-c7f9f8f34324
begin
	using HypergeometricFunctions
	using PyPlot
end

# ╔═╡ 02f4dd56-18d9-11eb-2d38-379eca82ec3a
begin
	rc("text", usetex = true)
	ion()
end

# ╔═╡ be283892-18d0-11eb-0154-758e466c7a40
begin
	const p = 2
	j = 16
	N = 2*j
end;

# ╔═╡ ffe361c6-18d0-11eb-2aa6-af741f44c0f9
function k_p(n,q)
	if n == q == 32
		return 0
	else
		return HypergeometricFunctions._₂F₁general2(-n,-q,-N,p)
	end
end;

# ╔═╡ f69b57c2-18da-11eb-21ea-0d37ac4811d6
function k_f(n,q)
	((-1)^n/2^j)*sqrt(binomial(2*j,n)*binomial(2*j,j+q))*k_p(n,j+q)
end;

# ╔═╡ 45185d6e-18db-11eb-2d2b-e3d194057cf4
begin
	vk_f = zeros(Float64,(N+1,N+1))
	for l in -j:j
		for m in 0:N
			vk_f[m+1,l+j+1] = k_f(m,l)
		end
	end
end

# ╔═╡ dd8123ba-18db-11eb-2c41-879bc043f44f
imshow(vk_f,cmap="gray")

# ╔═╡ Cell order:
# ╠═710bd51c-18d7-11eb-2a05-c7f9f8f34324
# ╠═02f4dd56-18d9-11eb-2d38-379eca82ec3a
# ╠═be283892-18d0-11eb-0154-758e466c7a40
# ╠═ffe361c6-18d0-11eb-2aa6-af741f44c0f9
# ╠═f69b57c2-18da-11eb-21ea-0d37ac4811d6
# ╠═45185d6e-18db-11eb-2d2b-e3d194057cf4
# ╠═dd8123ba-18db-11eb-2c41-879bc043f44f
