### A Pluto.jl notebook ###
# v0.12.11

using Markdown
using InteractiveUtils

# ╔═╡ baf9406a-2aed-11eb-1f5b-af8d76740822
begin
	using SpecialMatrices
	using LinearAlgebra
	using PyPlot
	using LaTeXStrings
end

# ╔═╡ 47f9b180-2af2-11eb-2f04-9b70cbe71785
begin
	rc("text", usetex = true)
	ion()
end

# ╔═╡ 51a7c7c4-2aef-11eb-2e97-27ad4437b92e
begin
	j = 16
	N = 2*j
	dim = N+1
end;

# ╔═╡ 649a4a96-2aef-11eb-04cb-c5bba93cebb1
begin
	Δ1 = zeros(Float64,dim)
	Δ1[1] = -2
	Δ1[2] = 1
	Δ1[dim] = 1
	Δ1 = SpecialMatrices.Circulant(Δ1)
end;

# ╔═╡ 6d3bff0a-2af1-11eb-0863-81b26e690aab
begin
	function ars(m)
		-4*sin(pi*m/N)^2
	end
	Δ2 = zeros(Float64,dim)
	for i in -j:j
		Δ2[i+1+j] = ars(i)
	end
	Δ2 = LinearAlgebra.diagm(Δ2)
end;

# ╔═╡ da00b368-2af1-11eb-3111-a1bde1b4f2dd
H = -(1/2)*(Δ1+Δ2);

# ╔═╡ 0e0bcffa-2af2-11eb-173b-b50a251e6015
EE = LinearAlgebra.eigvecs(H);

# ╔═╡ 2f2f1e9c-2af2-11eb-0eae-b3154cf93f17
begin
	fi1 = figure()
	ax1 = gca()
	ax1.imshow(EE,cmap="gray")
	ax1.set_xlabel(L"k",fontsize=12)
	ax1.set_ylabel(L"m",fontsize=12)
	ax1.set_xticks(range(0,N,step=j))
	ax1.set_yticks(range(0,N,step=j))
	ax1.set_yticklabels([-j,0,j])
end;

# ╔═╡ Cell order:
# ╠═baf9406a-2aed-11eb-1f5b-af8d76740822
# ╠═47f9b180-2af2-11eb-2f04-9b70cbe71785
# ╠═51a7c7c4-2aef-11eb-2e97-27ad4437b92e
# ╠═649a4a96-2aef-11eb-04cb-c5bba93cebb1
# ╠═6d3bff0a-2af1-11eb-0863-81b26e690aab
# ╠═da00b368-2af1-11eb-3111-a1bde1b4f2dd
# ╠═0e0bcffa-2af2-11eb-173b-b50a251e6015
# ╠═2f2f1e9c-2af2-11eb-0eae-b3154cf93f17
