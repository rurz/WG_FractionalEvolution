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

# ╔═╡ 649a4a96-2aef-11eb-04cb-c5bba93cebb1
function Δ1(j)
	Δ1 = zeros(Float64, 2*j + 1)
	Δ1[1] = -2
	Δ1[2] = 1
	Δ1[end] = 1
	Δ1 = SpecialMatrices.Circulant(Δ1)
end;

# ╔═╡ 6d3bff0a-2af1-11eb-0863-81b26e690aab
function Δ2(j)
	Δ2 = zeros(Float64,2*j+1)
	for i in -j:j
		Δ2[i+1+j] = -4 * sin(pi * i / (2*j))^2
	end
	return Δ2 = LinearAlgebra.diagm(Δ2)
end

# ╔═╡ da00b368-2af1-11eb-3111-a1bde1b4f2dd
@time H = -(1/2.0)*(Δ1(16) + Δ2(16));

# ╔═╡ 0e0bcffa-2af2-11eb-173b-b50a251e6015
@time EE = LinearAlgebra.eigvecs(H);

# ╔═╡ 2f2f1e9c-2af2-11eb-0eae-b3154cf93f17
@time begin
	fi1 = figure()
	ax1 = gca()
	ax1.imshow(EE,cmap="gray")
	ax1.set_xlabel(L"k",fontsize=12)
	ax1.set_ylabel(L"m",fontsize=12)
	ax1.set_xticks(range(0,32,step=16))
	ax1.set_yticks(range(0,32,step=16))
	ax1.set_yticklabels([-16,0,16])
end;

# ╔═╡ Cell order:
# ╠═baf9406a-2aed-11eb-1f5b-af8d76740822
# ╠═47f9b180-2af2-11eb-2f04-9b70cbe71785
# ╠═649a4a96-2aef-11eb-04cb-c5bba93cebb1
# ╠═6d3bff0a-2af1-11eb-0863-81b26e690aab
# ╠═da00b368-2af1-11eb-3111-a1bde1b4f2dd
# ╠═0e0bcffa-2af2-11eb-173b-b50a251e6015
# ╠═2f2f1e9c-2af2-11eb-0eae-b3154cf93f17
