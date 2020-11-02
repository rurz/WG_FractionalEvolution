### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ 1a6379ce-0ffc-11eb-0fce-853d047f252a
begin
	using PyPlot
	using PyCall #I'll gonna use SpecialFunctions.jl and HypergeometricFunctions.jl
	using LaTeXStrings
end

# ╔═╡ cd516e2c-0ffe-11eb-2415-832b4823d2ca
md"### Legacy notebook using PyCall and mpmath library

I will keep this record for future reference"

# ╔═╡ 3e643f7a-0ffc-11eb-021f-636a88842643
mp = pyimport("mpmath")

# ╔═╡ 5ccf83de-0ffc-11eb-2b02-8f1ff10f54af
begin
	rc("text", usetex=true)
	ion()
end

# ╔═╡ 61e49ab2-0ffc-11eb-157c-db7651aee125
begin
	N = 32
	dim = N+1
	const p = 1/2.0
end;

# ╔═╡ 67fd8670-0ffc-11eb-2b6e-ffefab121f31
function nbinom(n,k); # Binomial coefficient in terms of the Gamma function
    return 1/((n+1)*convert(Float64,mp.beta(k+1, n-k+1)))
end;

# ╔═╡ 6ef5f2ac-0ffc-11eb-2dcc-39e1895c1b70
function krav(i,j) # Kravchuk symmetric polynomial in terms of 2F1
    return convert(Float64,mp.hyp2f1(-j,-i,-N,1/p,zeroprec = 250))
end;

# ╔═╡ 74c15e90-0ffc-11eb-168c-afd38ba93c9b
function omega(j)
    return nbinom(N,j)*(p^j)*(1-p)^(N-j)
end;

# ╔═╡ 7a746828-0ffc-11eb-2c9e-87213ad35743
function h(i)
    return (((1-p)/p)^i)/(nbinom(N,i))
end;

# ╔═╡ 80054320-0ffc-11eb-1837-b7f067634a26
function K(i,j)
    return sqrt(omega(j)/h(i))*krav(i,j)
end;

# ╔═╡ 8804e044-0ffc-11eb-3cf9-efd4e3a8585d
begin
	KT = zeros(Float64,(dim,dim))
	@timed for l in 0:N
    	for m in 0:N
        	KT[l+1,m+1] = K(l,m)
    	end
	end
end

# ╔═╡ 8ed657c0-0ffc-11eb-1974-3154511ac71a
begin
	fig2 = figure(figsize=(5,5))
	ax2 = gca()
	ax2.imshow(KT, cmap="gray")
	title(L"\textbf{Kravchuk Matrix}\; K(i,j)\;\textbf{(The trilobite)}")
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
end

# ╔═╡ 36310bca-0ffd-11eb-1acf-dfc63763f5df
md"Using PyCall with mpmath library, we obtain, for N=32, a Kravchuk matrix of NxN in 0.673462 s. With the native libraries, Special and Hypergeometric, this same matrix is obtained in 0.0028436 s (see FrFK.jl). he native evaluation is around 237 times faster. 😲"

# ╔═╡ Cell order:
# ╠═cd516e2c-0ffe-11eb-2415-832b4823d2ca
# ╠═1a6379ce-0ffc-11eb-0fce-853d047f252a
# ╟─3e643f7a-0ffc-11eb-021f-636a88842643
# ╠═5ccf83de-0ffc-11eb-2b02-8f1ff10f54af
# ╠═61e49ab2-0ffc-11eb-157c-db7651aee125
# ╠═67fd8670-0ffc-11eb-2b6e-ffefab121f31
# ╠═6ef5f2ac-0ffc-11eb-2dcc-39e1895c1b70
# ╠═74c15e90-0ffc-11eb-168c-afd38ba93c9b
# ╠═7a746828-0ffc-11eb-2c9e-87213ad35743
# ╠═80054320-0ffc-11eb-1837-b7f067634a26
# ╠═8804e044-0ffc-11eb-3cf9-efd4e3a8585d
# ╠═8ed657c0-0ffc-11eb-1974-3154511ac71a
# ╟─36310bca-0ffd-11eb-1acf-dfc63763f5df
