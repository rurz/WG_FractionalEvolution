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
	using PyCall
end

# ╔═╡ 814bda9c-14c4-11eb-0ba8-a71a2d9c234b
md"### Playing with the discrete polynomials definitions"

# ╔═╡ bb8453e6-18b7-11eb-22dd-519ac2f2d0f6
mp = pyimport("mpmath")

# ╔═╡ e4232312-14c6-11eb-2578-d775876aa9f4
begin
	rc("text", usetex = true)
	ion()
end

# ╔═╡ 65817e78-14c6-11eb-2016-9d7cf09a2801
begin
	j = 16
	N = 2*j
	const p = 2
end

# ╔═╡ c6c1b848-18b7-11eb-03ab-ef48179d3b5d
function kravp_mp(n,q)
	convert(Float64,mp.hyp2f1(-n,-q,-N,p,zeroprec = 250))
end

# ╔═╡ 8f1e8a68-18b9-11eb-2888-7d8b6effde6a
function nbinom(n,k)
    1/((n+1)*beta(k+1, n-k+1))
end

# ╔═╡ faa722e2-18b7-11eb-1899-e5bdcf897db3
function kravf_mp(n,q)
	((-1)^n/(2^j))*sqrt(nbinom(2*j,n)*nbinom(2*j,j+q))*kravp_mp(n,q+j)
end

# ╔═╡ 11d8fc38-18b8-11eb-19fe-c739310c212c
begin
	K1v_mp = zeros(Float64,(N+1,N+1))
	for l in -j:j
		for m in 0:N
			K1v_mp[l+j+1,m+1] = kravf_mp(m,l)
		end
	end
end

# ╔═╡ d6f280d4-18b8-11eb-3240-1fb563728533
imshow(K1v_mp,cmap="gray")

# ╔═╡ 75b2d10e-18bb-11eb-0f1e-0fc1d7d51b07
begin
	Kn = zeros(Float64,(N+1,N+1))
	for l in -j:j
		for m in 0:N
			Kn[m+1,l+j+1] = kravf_mp(m,l)
		end
	end
end

# ╔═╡ a4fcc568-18bc-11eb-0f1d-a172c97d48d2
begin
	imshow(Kn,cmap="gray")
	xlabel("q")
	xticks([0,j,N],[-j,0,j])
	ylabel("n")
	yticks([0,j,N],[N,j,0])
end

# ╔═╡ f1ed4446-18bd-11eb-27b3-6da15e29c6d5
plot(rotl90(Kn)[N+1,:])

# ╔═╡ 65e59e8e-18be-11eb-3e41-cb3799c1cc46
plot(Kn[N+1,:])

# ╔═╡ Cell order:
# ╟─814bda9c-14c4-11eb-0ba8-a71a2d9c234b
# ╠═a70cdae2-14c4-11eb-0f58-3bf672d63847
# ╟─bb8453e6-18b7-11eb-22dd-519ac2f2d0f6
# ╠═e4232312-14c6-11eb-2578-d775876aa9f4
# ╠═65817e78-14c6-11eb-2016-9d7cf09a2801
# ╠═c6c1b848-18b7-11eb-03ab-ef48179d3b5d
# ╠═8f1e8a68-18b9-11eb-2888-7d8b6effde6a
# ╠═faa722e2-18b7-11eb-1899-e5bdcf897db3
# ╠═11d8fc38-18b8-11eb-19fe-c739310c212c
# ╠═d6f280d4-18b8-11eb-3240-1fb563728533
# ╠═75b2d10e-18bb-11eb-0f1e-0fc1d7d51b07
# ╠═a4fcc568-18bc-11eb-0f1d-a172c97d48d2
# ╠═f1ed4446-18bd-11eb-27b3-6da15e29c6d5
# ╠═65e59e8e-18be-11eb-3e41-cb3799c1cc46
