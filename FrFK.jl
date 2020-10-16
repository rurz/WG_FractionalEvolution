### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ 4b51367c-09ec-11eb-16d9-3dcd2910078e
begin
	using PyPlot
	#using PyCall #I'll gonna use SpecialFunctions.jl and HypergeometricFunctions.jl
	using LaTeXStrings
	using SpecialFunctions
	using HypergeometricFunctions
end

# ╔═╡ eb325bf2-09e3-11eb-3b88-9b0d15626515
md"## Fractional Fourier-Kravchuk evolution on waveguides"

# ╔═╡ 74a12d8e-0f69-11eb-0650-e79c06fd1080
md"Kravchuk matrix is defined as

$K(i,j)=\sqrt{\frac{\omega(j)}{h(i)}}k(i,j),$
with the functions $\omega$ and $h$ given by

$\omega(j)=\binom{N}{j}p^{j}(1-p)^{N-j},\quad h(i)=\frac{\frac{1-p}{p^{i}}}{\binom{N}{i}},$
and the Kravchuk polynomials obtained by the Hypergeometrical function

$k(i,j)=~_{2}F_{1}(-j,-i,-N,1/p).$

Note: When $p=1/2$, the symmetric Kravchuk polynomials are obtained, henceforth, the symmetric Kravchuk functions follows.
"

# ╔═╡ 7f8aad74-09ec-11eb-3c25-a538f18da2f9
begin
	rc("text", usetex=true)
	ion()
end

# ╔═╡ c28078de-09ec-11eb-0c25-a1bd3c7cdd8b
md" Define the parameters"

# ╔═╡ cb9e6110-09ec-11eb-29bb-b33f49eb026a
begin
	N = 32
	dim = N+1
	const p = 1/2.0
end;

# ╔═╡ bb1e9a12-0fe6-11eb-0551-137403aecef6
function bspecial(n,k)
	return 1/((n+1)*beta(k+1,n-k+1))
end;

# ╔═╡ ee7b88a2-0fe6-11eb-3465-8b70e58f07c1
function kravs(i,j)
	return HypergeometricFunctions._₂F₁maclaurin(-j,-i,-N,1/p)
end;

# ╔═╡ 4f876a76-0fe7-11eb-3f46-0f18d733b446
function omegas(j)
	return bspecial(N,j)*(p^j)*(1-p)^(N-j)
end;

# ╔═╡ 590a100a-0fe7-11eb-017c-4d4ff299ed8a
function hs(i)
	return (((1-p)/p)^i)/(bspecial(N,i))
end;

# ╔═╡ a5f2bf3c-0fe7-11eb-3377-8d5228bb6476
function Ks(i,j)
	return sqrt(omegas(j)/hs(i))*kravs(i,j)
end;

# ╔═╡ d01de678-0fe8-11eb-2831-89f652e45a9a
begin
	KTs = zeros(Float64,(dim,dim))
	@timed for l in 0:N
    	for m in 0:N
        	KTs[l+1,m+1] = Ks(l,m)
    	end
	end
end

# ╔═╡ 4692817e-0fe9-11eb-22e0-dfdab6509141
begin	
	fig1 = figure(figsize=(5,5))
	ax1 = gca()
	ax1.imshow(KTs, cmap="gray")
	ax1.set_title(L"\textbf{Kravchuk Matrix}\; K(i,j)\;\textbf{(The trilobite)}")
	ax1.set_xlabel(L"\textbf{Pseudo Energy Level:}\; 0\leq i\leq N")
	ax1.set_ylabel(L"\textbf{Position Index:}\; 0\leq j\leq N")
	ax1.tick_params(direction="out",length=5,width=2,labelsize=10)
	ax1.set_xticks(0:8:N, minor = false)
	ax1.set_yticks(0:8:N, minor = false)
	
	axins1a = ax1.inset_axes([0, 0.8, 0.2, 0.2])
	axins1a.plot(KTs[1,:],"k:",linewidth=1)
	axins1a.set_xticklabels("")
	axins1a.set_yticklabels("")
	axins1a.text(7.5,0.01,L"K(0,j)",fontsize=9)
	
	axins2a = ax1.inset_axes([0.8, 0, 0.2, 0.2])
	axins2a.plot(KTs[N+1,:],"k:",linewidth=1)
	axins2a.set_xticklabels("")
	axins2a.set_yticklabels("")
	axins2a.text(15,-0.35,L"K(N,j)",fontsize=9)
	
	tight_layout()
	savefig("trilo.png",dpi=300,transparent=true)
end

# ╔═╡ 33089668-09ed-11eb-30db-8f33d969a115
md"
October 16

Status: Kravchuk matrix completed. Changed definitions of functions using native Julia libraries

Next: Write the fractional evolution, and plotting the results
"

# ╔═╡ Cell order:
# ╠═eb325bf2-09e3-11eb-3b88-9b0d15626515
# ╟─74a12d8e-0f69-11eb-0650-e79c06fd1080
# ╠═4b51367c-09ec-11eb-16d9-3dcd2910078e
# ╠═7f8aad74-09ec-11eb-3c25-a538f18da2f9
# ╟─c28078de-09ec-11eb-0c25-a1bd3c7cdd8b
# ╠═cb9e6110-09ec-11eb-29bb-b33f49eb026a
# ╠═bb1e9a12-0fe6-11eb-0551-137403aecef6
# ╠═ee7b88a2-0fe6-11eb-3465-8b70e58f07c1
# ╠═4f876a76-0fe7-11eb-3f46-0f18d733b446
# ╠═590a100a-0fe7-11eb-017c-4d4ff299ed8a
# ╠═a5f2bf3c-0fe7-11eb-3377-8d5228bb6476
# ╠═d01de678-0fe8-11eb-2831-89f652e45a9a
# ╠═4692817e-0fe9-11eb-22e0-dfdab6509141
# ╟─33089668-09ed-11eb-30db-8f33d969a115
