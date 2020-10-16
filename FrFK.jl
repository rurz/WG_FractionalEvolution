### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ 4b51367c-09ec-11eb-16d9-3dcd2910078e
begin
	using PyPlot
	using PyCall
	using LaTeXStrings
end

# ╔═╡ eb325bf2-09e3-11eb-3b88-9b0d15626515
md"## Fractional Fourier-Kravchuk evolution on waveguides"

# ╔═╡ 54542a04-09ec-11eb-1ba0-3b5eb440ac46
mp = pyimport("mpmath")

# ╔═╡ 7f8aad74-09ec-11eb-3c25-a538f18da2f9
rc("text", usetex=true)

# ╔═╡ c28078de-09ec-11eb-0c25-a1bd3c7cdd8b
md" Define the parameters"

# ╔═╡ cb9e6110-09ec-11eb-29bb-b33f49eb026a
begin
	N = 32
	dim = N+1
	const p = 1/2
end;

# ╔═╡ 839db44e-09ec-11eb-33e3-01af1f76d1b8
function nbinom(n,k); # Binomial coefficient in terms of the Gamma function
    return 1/((n+1)*convert(Float64,mp.beta(k+1, n-k+1)))
end

# ╔═╡ a73e16c4-09ec-11eb-0f0f-d9772f75fb9c
function krav(i,j) # Kravchuk symmetric polynomial in terms of 2F1
    return convert(Float64,mp.hyp2f1(-j,-i,-N,1/p,zeroprec = 250))
end

# ╔═╡ f0f453ac-09ec-11eb-35bf-5db9aac74d95
function omega(j)
    return nbinom(N,j)*(p^j)*(1-p)^(N-j)
end

# ╔═╡ f2629280-09ec-11eb-34b9-61a08407e9c7
function h(i)
    return (((1-p)/p)^i)/(nbinom(N,i))
end

# ╔═╡ f7b116c6-09ec-11eb-2b95-61db2b1a7196
function K(i,j)
    return sqrt(omega(j)/h(i))*krav(i,j)
end

# ╔═╡ fed998ce-09ec-11eb-2a73-4f900ef78a80
begin
	KT = zeros(Float64,(dim,dim))
	for l in 0:N
    	for m in 0:N
        	KT[l+1,m+1] = K(l,m)
    	end
	end
end

# ╔═╡ 4dd56d14-0f64-11eb-03d3-ef32481e399b
begin
	fig2 = figure(figsize=(5,5))
	ax2 = gca()
	ax2.imshow(KT, cmap="gray")
	title(L"\textbf{Kravchuk Matrix}\; K(i,j)\; \textbf{(The trilobite)}")
	ax2.set_xlabel(L"\textbf{Pseudo Energy Level:}\; 0\leq i\leq N")
	ax2.set_ylabel(L"\textbf{Position Index:}\; 0\leq j\leq N")
	ax2.tick_params(direction="out",length=5,width=2,labelsize=10)
	ax2.set_xticks(0:8:N, minor = false)
	ax2.set_yticks(0:8:N, minor = false)
	
	axins1 = ax2.inset_axes([0, 0.8, 0.2, 0.2])
	axins1.plot(KT[1,:],"k:",linewidth=1)
	axins1.set_xticklabels("")
	axins1.set_yticklabels("")
	axins1.text(7.5,0.01,L"K(0,j)",fontsize=9)
	
	axins2 = ax2.inset_axes([0.8, 0, 0.2, 0.2])
	axins2.plot(KT[N+1,:],"k:",linewidth=1)
	axins2.set_xticklabels("")
	axins2.set_yticklabels("")
	axins2.text(15,-0.35,L"K(N,j)",fontsize=9)
	
	tight_layout()
	#savefig("trilo.png",dpi=300,transparent=true)
end

# ╔═╡ 74a12d8e-0f69-11eb-0650-e79c06fd1080
md"Kravchuk matrix is defined as

$K(i,j)=\sqrt{\frac{\omega(j)}{h(i)}}k(i,j),$
with the functions $\omega$ and $h$ given by

$\omega(j)=\binom{N}{j}p^{j}(1-p)^{N-j},\quad h(i)=\frac{\frac{1-p}{p^{i}}}{\binom{N}{i}},$
and the Kravchuk polynomials obtained by the Hypergeometrical function

$k(i,j)=~_{2}F_{1}(-j,-i,-N,1/p).$

Note: When $p=1/2$, the symmetric Kravchuk polynomials are obtained, henceforth, the symmetric Kravchuk functions follows.
"

# ╔═╡ 33089668-09ed-11eb-30db-8f33d969a115
md"
October 9

Status: Kravchuk matrix completed

Next: Write the fractional evolution, and plotting the results
"

# ╔═╡ Cell order:
# ╠═eb325bf2-09e3-11eb-3b88-9b0d15626515
# ╠═4b51367c-09ec-11eb-16d9-3dcd2910078e
# ╟─54542a04-09ec-11eb-1ba0-3b5eb440ac46
# ╠═7f8aad74-09ec-11eb-3c25-a538f18da2f9
# ╟─c28078de-09ec-11eb-0c25-a1bd3c7cdd8b
# ╠═cb9e6110-09ec-11eb-29bb-b33f49eb026a
# ╠═839db44e-09ec-11eb-33e3-01af1f76d1b8
# ╠═a73e16c4-09ec-11eb-0f0f-d9772f75fb9c
# ╠═f0f453ac-09ec-11eb-35bf-5db9aac74d95
# ╠═f2629280-09ec-11eb-34b9-61a08407e9c7
# ╠═f7b116c6-09ec-11eb-2b95-61db2b1a7196
# ╠═fed998ce-09ec-11eb-2a73-4f900ef78a80
# ╠═4dd56d14-0f64-11eb-03d3-ef32481e399b
# ╟─74a12d8e-0f69-11eb-0650-e79c06fd1080
# ╟─33089668-09ed-11eb-30db-8f33d969a115
