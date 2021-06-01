    
# function bw(L::Lattice)
#     n = prod(L.dims)
#     bw = -zeros(Int,n)
#     for i in L.site
#         i == 1 && (bw[i] = 0)
#         for j in L.neig[i]
#             bw[j] = bw[i] > 0 ? 0 : 1
#         end
#     end
#     idxw = findall(bw .== 0)
#     idxb = findall(bw .== 1)
#     sort(vcat(idxw,idxb)) == collect(1:n) || error("stocazzo")
#     shuffle!(idxw),shuffle!(idxb)
# end


struct Ising{T<:AbstractFloat,D}
    N::Int               # number of spin
    β::T                 # inverse temperature
    Λ::Lattice{D}        # Lattice
    spin::Vector{Int}    # spin config
    h::Vector{T}         # external field
    H::Vector{T}         # total field
end


mutable struct Measures{T<:AbstractFloat}
    it::Int
    β::T
    ene::Vector{T}
    mag::Vector{T}
    conf::Vector{BitArray{1}}
end

function Measures(n::Int,β::T) where T
    Measures(0,β,zeros(T,n),zeros(T,n),Array{BitArray{1}}(undef,n))
end

spin2bool = x -> x > 0
bool2spin = x -> x > 0 ? 1 : -1
const mts = MersenneTwister.(1:Threads.nthreads())
function mcising(I::NTuple{N,Int},β::T;
                 h::Vector{T}=T[],
                 x0::Vector{Int}=rand([-1,1],prod(I)),
                 nterm::Int=100,
                 nmeas::Int=1000,
                 nsweep::Int=100,
                 verbose::Bool=false) where {T<:Real,N}

    
    ising = Ising(I,h,β,x0)
    output = Measures(nmeas,β)
    t=0
    for t in 1:nterm
        onemcsweep!(ising)
    end


    for _ in 1:nmeas
        for _ in 1:nsweep
            onemcsweep!(ising)
        end
        measures!(output,ising)
    end
    output
end

function measures!(ou::Measures{T},x::Ising) where T
    ou.it += 1
    ou.ene[ou.it] = energy(x)
    ou.mag[ou.it] = mean(x.spin)
    ou.conf[ou.it] = BitArray(map(spin2bool,x.spin))
end

function onemcsweep!(x::Ising)
    for site in randperm(x.N)
        onemcstep!(x,site)
    end
end

function onemcsweepthreads!(x::Ising,idxw,idxb)
#    for site in randperm(x.N)
    for site in idxw
        onemcstep!(x,site)
    end
    for site in idxb
        onemcstep!(x,site)
    end
end

function onemcstep!(x::Ising, site::Int)
    @inbounds begin
        ΔE = 2.0 * x.spin[site] * x.H[site]
        if ΔE <= 0 || rand(mts[Threads.threadid()]) < exp( -x.β * ΔE )
            x.spin[site] *= -1
            nn = x.Λ.neig[site]
            for neigsite in nn
                x.H[neigsite] += 2 * x.spin[site]
            end
        end
    end
end

function Ising(I::NTuple{D,Int},h::Vector{T},β::T, x0::Vector{Int})  where D where T<:AbstractFloat
    N = prod(I)
    length(x0) != N && error("number of spin incompatible with lattice size = $(Λ.dims)")
    nh = length(h)
    if nh == 0
        h = zeros(T,N)
    elseif nh > 0 && nh !=N
        error("size of h should be $N ≠ $nh")
    end
    Λ = Lattice(I)
    H = zeros(T,N)
    @inbounds for i in 1:N
        _s = 0
        nn = Λ.neig[i]
        @simd for siteneig in nn
            _s += x0[siteneig]
        end
        H[i] = h[i] + _s
    end
    return Ising(N,β,Λ,x0,h,H)
end

function energy(x::Ising{T,D}) where {T,D}
    ene = zero(T)
    @inbounds @simd for i in 1:x.N
        ene -= (x.H[i] + x.h[i]) * x.spin[i]
    end
    return 0.5 * ene
end

function energyslow(x::Ising{T,D}) where {T,D}
    ene = zero(T)
    @inbounds for i in 1:x.N
        ene -= x.h[i]*x.spin[i]
        nn = x.Λ.neig[i]
        _s = 0
        @simd for siteneig in nn
            _s += x.spin[siteneig]
        end
        ene -=  0.5 * x.spin[i] * _s
    end
    return ene
end
