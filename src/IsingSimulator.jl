module IsingSimulator
using DataFrames
using Base: BufferStream
using Random, Statistics,GLMakie, QuadGK
using ProgressMeter
import ForwardDiff: derivative

export isingview,mytest,βc2d,magn2d,free_energy2d,observables


include("lattice.jl")
include("mc.jl")
include("isingview.jl")
include("thermodynamics.jl")

const βc2d = log(1+√2)/2 # critical temperature 2d ising FM model
inline!(true)
end # end module
