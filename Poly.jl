### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ a70cdae2-14c4-11eb-0f58-3bf672d63847
begin
	using SpecialFunctions
	using SpecialPolynomials
	using HypergeometricFunctions
	using PyPlot
	using Polynomials
end

# ╔═╡ 814bda9c-14c4-11eb-0ba8-a71a2d9c234b
md"### Playing with the discrete polynomials definitions"

# ╔═╡ e4232312-14c6-11eb-2578-d775876aa9f4
begin
	rc("text", usetex = true)
	ion()
end

# ╔═╡ 65817e78-14c6-11eb-2016-9d7cf09a2801
begin
	N = 16
	const p = 1/2.0
end

# ╔═╡ be834efe-14c4-11eb-0263-0303df20c496
p1 = SpecialPolynomials.Krawchouk{2,N}([1,1,1])

# ╔═╡ 7e59cd32-14ef-11eb-3eea-0d489366914d
convert(Polynomial,p1)

# ╔═╡ fbcb8686-14c5-11eb-25fd-cbb092997516
function k1(n,q)
	HypergeometricFunctions._₂F₁(-n,-q,N,-(1/p))
end

# ╔═╡ bf6bcaa8-14ec-11eb-1683-0967b3046871
function k2(n,q)
	HypergeometricFunctions._₂F₁general(-n,-q,N,-(1/p))
end

# ╔═╡ 75cd1616-14c6-11eb-2d89-5fa0db29b867
begin
	vk1 = zeros(Float64,(N+1,N+1))
	for i in 0:N
		for j in 0:N
			vk1[i+1,j+1] = k1(i,j)
		end
	end
end

# ╔═╡ f9021cb8-14ec-11eb-123f-670f888615f4
begin
	vk2 = zeros(Float64,(N+1,N+1))
	for i in 0:N
		for j in 0:N
			vk2[i+1,j+1] = k2(i,j)
		end
	end
end

# ╔═╡ e0e3744a-14c6-11eb-191b-3f127abe71cf
function K1(n,q)
	((-1)^(n))/(2^(N))*sqrt(binomial(2*N,n)*binomial(2*N,N+q))*k1(n,N+q)
end

# ╔═╡ 71fdb04a-14ee-11eb-08a9-c34df27a1892
begin
	K1v = zeros(Float64,2*N+1)
	for i in -N:N
		K1v[i+1+N] = K1(0,i)
	end
end

# ╔═╡ 938345cc-14ee-11eb-320d-a7ca774f8628
plot(K1v)

# ╔═╡ Cell order:
# ╟─814bda9c-14c4-11eb-0ba8-a71a2d9c234b
# ╠═a70cdae2-14c4-11eb-0f58-3bf672d63847
# ╠═e4232312-14c6-11eb-2578-d775876aa9f4
# ╠═be834efe-14c4-11eb-0263-0303df20c496
# ╠═7e59cd32-14ef-11eb-3eea-0d489366914d
# ╠═65817e78-14c6-11eb-2016-9d7cf09a2801
# ╠═fbcb8686-14c5-11eb-25fd-cbb092997516
# ╠═bf6bcaa8-14ec-11eb-1683-0967b3046871
# ╠═75cd1616-14c6-11eb-2d89-5fa0db29b867
# ╠═f9021cb8-14ec-11eb-123f-670f888615f4
# ╠═e0e3744a-14c6-11eb-191b-3f127abe71cf
# ╠═71fdb04a-14ee-11eb-08a9-c34df27a1892
# ╠═938345cc-14ee-11eb-320d-a7ca774f8628
