using SpecialMatrices
using LinearAlgebra
using SparseArrays
using PyPlot
using LaTeXStrings

""" SpecialMatrices run out an error when the register package is used,
so we pkg> add the #master branch """

rc("text", usetex = true)
ion()

const j = 16
const N = 2 * j
const dim = N + 1

function Δ1()
	Δ1 = zeros(Float64, dim)
	Δ1[1] = -2
	Δ1[2] = 1
	Δ1[end] = 1
	Δ1 = SpecialMatrices.Circulant(Δ1)
end

function Δ2()
	Δ2 = zeros(Float64,dim)
	for i in -j:j
		Δ2[i + 1 + j] = -4 * (sin(pi * i / N))^2
	end
	return Δ2 = LinearAlgebra.diagm(Δ2)
end

H = -(1/2.0)*(Δ1() + Δ2());

h = LinearAlgebra.eigvecs(H);

k = LinearAlgebra.eigvals(H);

η= [l for l in 0:N];

G(m,μ,t) = sum([h[m + j + 1, l + 1] * exp(-1im * pi * t * η[l + 1]/2.0) * h[μ + j + 1, l + 1] for l in 0:N]);

P(t) = reshape([G(m,μ,t) for m in -j:j for μ in -j:j], (dim, dim))

function P2(t)
	p2 = zeros(Float64, (dim, dim))
	for μ in -j:1:j
		for m in -j:1:j
			p2[m + j + 1, μ + j + 1] = G(m, μ, t)
		end
	end
	return p2
end


pint(q) = sparsevec(Dict(q + j + 1 => 1), dim);

eigh(n) = h[:,n]


function hdc(d)
	ne = zeros(Float64,dim)
	for l in 1:dim-d
		ne[l] = h[l+d,1]
	end
	return ne
end

α = range(0, stop = 2, length = 315);

evol(θ) = abs2.(P(α[θ]) * normalize(hdc(8)));

evolT = [evol(l) for l in 1:315];

Ankara = zeros(Float64, (dim, 315))

for l in 1:dim
	for m in 1:315
		Ankara[l, m] = evolT[m][l]
	end
end

begin
	fig1 = figure(figsize = (10,5))
	ax1 = gca()
	ax1.imshow(Ankara, cmap = "hot")
	ax1.set_aspect("auto")

	tight_layout()
end
