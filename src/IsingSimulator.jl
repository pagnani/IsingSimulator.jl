module IsingSimulator
using Base: BufferStream
using Random, Statistics,GLMakie, QuadGK


export isingview!,mytest,βc2d,magn2d,entropy,ene,free_energy


include("lattice.jl")
include("mc.jl")
include("isingview.jl")


const βc2d = log(1+√2)/2
inline!(true)
magn2d(x) = x > βc2d ? (1-(1/sinh(2x))^4)^(1/8) : 0.0
entropy(m) = -0.5*((1+m)*log(0.5*(1+m)) + (1-m)*log(0.5(1-m)))

# K = βJ L=βJ' J = J' for isotropic ising 2d

function fe_density(β,θ)
    k = 2sinh(2β)/cosh(2β)^2
    return log(1+sqrt(1-k^2 * cos(θ)^2))
end

function free_energy(β)
    f(θ) = fe_density(β,θ)
    res = QuadGK.quadgk(f,0,π)[1]
    return -log(2)/2 - log(cosh(2β))-res/(2π)
end

function free_and_mag(binf,bsup;nstep=1024)
    vals = LinRange(binf,bsup,nstep)
    m = [magn2d(β) for β in vals]
    f = [free_energy(β) for β in vals]
    return m,f
end


end # end module