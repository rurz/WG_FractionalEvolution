### A Pluto.jl notebook ###
# v0.12.11

using Markdown
using InteractiveUtils

# ╔═╡ 1619652c-2b63-11eb-2ec4-9bdbd10c1539
begin
	using SpecialFunctions
	using PyPlot
	using LaTeXStrings
	rc("text", usetex = true)
	ion()
end

# ╔═╡ 3f123684-2b63-11eb-23e2-9749d9d79fea
using HypergeometricFunctions

# ╔═╡ 5dc79c5c-2b88-11eb-240c-8bc0041210f6
using SpecialPolynomials

# ╔═╡ db77ccb8-2b5f-11eb-367b-2332d89b4e97
md"### Testing Julia packages and methods for producing Jacobi functions used to represent the eigenfunctions of the discrete and finite harmonic oscillator

##### Theory
In their 2016 article, Steffen Weimann _et. al_ [1], used the eigenfunctions of the angular momentum coupling $\hat{J}_{x}$ to represent the modes of the finite and discrete harmonic oscillator as

$$u_{n}^{q} = 2^{n}\sqrt{\frac{(j-n)!(j+n)!}{(j-q)!(j+q)!}}P_{n}^{q-n,-q-n}(0),$$

where $P_{n}^{\alpha,\beta}(0)$ is the Jacobi function at $z=0$; $n$ and $q$ are the _energy levels_ and _position index_, respectively such that $n, q\in\{-j,...,j\}$, for $j\in Z^{+}.$

##### The Task
For the purpose, we want to build a $d_{j}\times d_{j}$ matrix for $d_{j}=N+1=2j+1$ where the entries are labeled by the pair $(n,q)$. So, we need the numerical evaluation of the Jacobi function around $z=0$.
__________

[1] Weimann, S., Perez-Leija, A., Lebugle, M. et al. Implementation of quantum and classical discrete fractional Fourier transforms. Nat Commun 7, 11027 (2016). [https://doi.org/10.1038/ncomms11027](https://doi.org/10.1038/ncomms11027)
"

# ╔═╡ e5549ccc-2b62-11eb-24af-bd2ef6a166d7
begin
	# Since most of the representation of the discrete and finite harmonic oscillator are based on the analysis of the angular momentum algebras, we adopt the parameters
	j = 8 # Angular momentum number
	N = 2*j
	dim = N+1
end;

# ╔═╡ 64ad621c-2b61-11eb-339d-95be29513595
md"##### (I) Using `HypergeometricFunctions.jl`

Since this package hasn't implemented the function $P_{n}^{\alpha,\beta}(z)$, we use the definition in terms of the hypergeometric function $_{2}F_{1}(n,a,b;c)$ as

$$P_{n}^{\alpha,\beta}(z)\equiv\binom{n+\alpha}{n}_{2}F_{1}(-n,-n+\alpha+\beta+1,\alpha+1;\tfrac{1}{2}(1-z)).$$

Where we expect a nicely behaviour around $z=0$. The issue as you can see is that $\alpha,\beta>-1$, who is something not fullfilled when $\alpha=-q-n$.
"

# ╔═╡ 12149bc0-2b65-11eb-1152-7159f82080b3
function JHyp(n::Integer,α::Integer,β::Integer,z)
	binomial(n+α,n)*HypergeometricFunctions._₂F₁(-n,n+α+β+1,α+1,(1-z)/2)
end;

# ╔═╡ 63b06a18-2b65-11eb-02b4-0147f09abdfd
function JHypF(n::Integer,q::Integer)
	(2^Float64(n))*sqrt(gamma(big(j-n+1))*gamma(big(j+n+1))/(gamma(big(j-q+1))*gamma(big(j+q+1))))*JHyp(j+n,q-n,-q-n,0)
end;

# ╔═╡ 9946f82c-2b65-11eb-1472-2355273a3f18
@timed begin
	JHypM = zeros(Float64,(dim,dim))
	for l in -j:j
		for m in -j:j
			JHypM[l+j+1,m+j+1] = JHypF(l,m)
		end
	end
end

# ╔═╡ d8159ab8-2b65-11eb-3fb0-674f6f8d069b
begin
	imshow(JHypM,cmap="gray")
	savefig("JHypM.png",dpi=150)
end

# ╔═╡ 0808796e-2b9a-11eb-28bb-898219baca42
md"![](https://raw.githubusercontent.com/rurz/WG_FractionalEvolution/main/tests/JHypM.png)"

# ╔═╡ 47fcb3dc-2b66-11eb-0322-7160add87fd2
md"When we plot the matrix array, we see that only the upper left half of the values are obtained. Meaning that the method provided by `HypergeometricFunctions.jl` lacks on rules of parity, simmetry, and reflection. Examining the matrix we encounter NaN results, so we faced a problem of ill-defined numerical method when we don't fullfill the condition $\alpha,\beta>-1.$"

# ╔═╡ 4c6ceec6-2b88-11eb-36f8-f3992615474f
md"#### (II) Using `SpecialPolynomials.jl`
This package in particular takes two parameters, (α,β), and a vector of $n$ dimension to produce the related polynomial on $z$ who after needs to be evaluated.\
$$J_{n}^{\alpha,\beta}(z)\equiv$$ `Jacobi{α,β}(V)`, where $dim(V)=n\in\{0,...,N\}.$
"

# ╔═╡ 08fd30ea-2b8a-11eb-1474-ef09e8030734
begin
	function vn(n)
		vd = zeros(Int64,n+1)
		vd[end]=1
		return vd
	end
	function JSP(n,α,β,z)
		(SpecialPolynomials.Jacobi{α,β}(vn(n)))(z)
	end
end;

# ╔═╡ 07566168-2b8b-11eb-1067-b3327e8c3c06
function JSPF(n,q)
	(2^Float64(n))*sqrt(gamma(big(j-n+1))*gamma(big(j+n+1))/(gamma(big(j-q+1))*gamma(big(j+q+1))))*JSP(j+n,q-n,-q-n,0)
end;

# ╔═╡ 1f9a2ca8-2b8b-11eb-14a5-7f83166acde4
# Uncomment for evaluate. Carefull, it requires a lot of time to perform the calculations for j above 6
begin
	JSPFM = zeros(Float64,(dim,dim))
	for l in -j:j
		for m in -j:j
			JSPFM[l+j+1,m+j+1] = JSPF(l,m)
		end
	end
end

# ╔═╡ 80b72b96-2b8e-11eb-3150-6f7795a2d027
begin
	imshow(JSPFM,cmap="gray")
	savefig("JSPFM.png",dpi=150)
end

# ╔═╡ f17428a6-2b99-11eb-2e6f-73737235bc48
md"![](https://raw.githubusercontent.com/rurz/WG_FractionalEvolution/main/tests/JSPFM.png)"

# ╔═╡ 3813f4e0-2b94-11eb-3387-cf14cad21276
md"
In this particular implementation we see that a problem arise in the lower middle part of the matrix $(n,q)$, implying that the definition of the `Jacobi{α,β}` suffer an error for some combination of values, as we can see in the error message above. Moreover, the evaluation it's far from be fast because we need to evaluate every polynomial at the time.
"

# ╔═╡ d9e0a670-2b66-11eb-3181-e3a029a778cb
md"##### (III) Using series representation expansion around $z=0$

Provided that we can obtain a closed form of $P_{n}^{\alpha,\beta}(z)$ for $z\rightarrow 0$, we use this approach. Given in [Wolfram Function Library](https://functions.wolfram.com/Polynomials/JacobiP/06/01/04/01/0004/), the definition used is

$$P_{n}^{\alpha,\beta}(z)=2^{-n}\sum\limits_{k=0}^{n}\binom{\alpha+n}{k}\binom{\beta+n}{n-k}(z+1)^{k}(z-1)^{n-k},$$

who isn't ill-defined at $z\rightarrow 0$, and we expect to behave nicely for values $\alpha,\beta<-1$.
"

# ╔═╡ 980b12e2-2960-11eb-220a-7786e56f2261
function JSer(n::Integer,α::Integer,β::Integer,z)
	JSerM = zeros(Float64,n+1)
	for k in 0:n
		JSerM[k+1] = binomial(α+n,k)*binomial(β+n,n-k)*(z+1)^(k)*(z-1)^(n-k)
	end
	return (1/2^(n))*sum(JSerM)
end;

# ╔═╡ f2509478-2566-11eb-2aa0-0be2f2b03b3f
function JSerF(n::Integer,q::Integer)
	(2^Float64(n))*sqrt(gamma(big(j-n+1))*gamma(big(j+n+1))/(gamma(big(j-q+1))*gamma(big(j+q+1))))*JSer(j+n,q-n,-q-n,0)
end;

# ╔═╡ 49243b96-2955-11eb-1d0e-a3f83590722e
@timed begin
	JSerFM = zeros(Float64,(dim,dim))
	for l in -j:j
		for m in -j:j
			JSerFM[l+j+1,m+j+1] = JSerF(l,m)
		end
	end
end

# ╔═╡ 9b37f096-2962-11eb-0512-b702db069ddf
begin
	imshow(JSerFM,cmap="gray")
	savefig("JSerFM.png",dpi=150)
end

# ╔═╡ 824608be-2b99-11eb-3b4e-1752fd403812
md"![K](https://raw.githubusercontent.com/rurz/WG_FractionalEvolution/main/tests/JSerFM.png)"

# ╔═╡ 4fc55ac4-2b68-11eb-124e-675dece08b7c
md"This time, we obtain the matrix in complete form. Since this approach use a series representation, we have not problems dealing with symmetries and reflection changes. We then have a nicely numerical way to obtain the matrix. The cons are that we need to construct the series and stored in a function `JSer`, who need to be evaluated every time whe change the value of $j\in Z^{+}$."

# ╔═╡ 0045571e-2b69-11eb-18fa-9bb5003e3968
md"
**Aftermath**: A revision on the methods of `Hypergeometricfunctions.jl` and `SpecialPolynomials.jl` would be nice to overcome this problem and have a fully implementable function, who for sure, responds on this particular problem at the time.
"

# ╔═╡ de05e088-2b68-11eb-274c-b9be94d7f0a7
md"
Alejandro R. Urzúa (2020)

[Github](https://github.com/rurz)
"

# ╔═╡ Cell order:
# ╟─db77ccb8-2b5f-11eb-367b-2332d89b4e97
# ╠═1619652c-2b63-11eb-2ec4-9bdbd10c1539
# ╠═e5549ccc-2b62-11eb-24af-bd2ef6a166d7
# ╟─64ad621c-2b61-11eb-339d-95be29513595
# ╠═3f123684-2b63-11eb-23e2-9749d9d79fea
# ╠═12149bc0-2b65-11eb-1152-7159f82080b3
# ╠═63b06a18-2b65-11eb-02b4-0147f09abdfd
# ╠═9946f82c-2b65-11eb-1472-2355273a3f18
# ╠═d8159ab8-2b65-11eb-3fb0-674f6f8d069b
# ╟─0808796e-2b9a-11eb-28bb-898219baca42
# ╟─47fcb3dc-2b66-11eb-0322-7160add87fd2
# ╟─4c6ceec6-2b88-11eb-36f8-f3992615474f
# ╠═5dc79c5c-2b88-11eb-240c-8bc0041210f6
# ╠═08fd30ea-2b8a-11eb-1474-ef09e8030734
# ╠═07566168-2b8b-11eb-1067-b3327e8c3c06
# ╠═1f9a2ca8-2b8b-11eb-14a5-7f83166acde4
# ╠═80b72b96-2b8e-11eb-3150-6f7795a2d027
# ╟─f17428a6-2b99-11eb-2e6f-73737235bc48
# ╟─3813f4e0-2b94-11eb-3387-cf14cad21276
# ╟─d9e0a670-2b66-11eb-3181-e3a029a778cb
# ╠═980b12e2-2960-11eb-220a-7786e56f2261
# ╠═f2509478-2566-11eb-2aa0-0be2f2b03b3f
# ╠═49243b96-2955-11eb-1d0e-a3f83590722e
# ╠═9b37f096-2962-11eb-0512-b702db069ddf
# ╟─824608be-2b99-11eb-3b4e-1752fd403812
# ╟─4fc55ac4-2b68-11eb-124e-675dece08b7c
# ╟─0045571e-2b69-11eb-18fa-9bb5003e3968
# ╟─de05e088-2b68-11eb-274c-b9be94d7f0a7
