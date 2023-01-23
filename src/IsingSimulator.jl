module IsingSimulator
using DataFrames
import Base: BufferStream
using Random, Statistics,GLMakie, QuadGK
import ProgressMeter: Progress, next!
import ForwardDiff: derivative
export isingview,mytest,βc2d,magn2d,free_energy2d,observables
#import AbstractPlotting:inline!

#inline!(true)


include("lattice.jl")
include("mc.jl")
include("isingview.jl")
include("thermodynamics.jl")

const βc2d = log(1+√2)/2 # critical temperature 2d ising FM model


end # end module
