# IsingSimulator

To run

```
julia> N=64; isingview((N,N),βc2d+0.01,nsweep=100000,fps=10^8,x0=ones(Int,N^2))
```