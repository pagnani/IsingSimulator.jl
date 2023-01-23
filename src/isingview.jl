
function isingview(I::NTuple{N,Int},β::T;
    h::Vector{T}=T[],
    x0::Vector{Int}=rand([-1,1],prod(I)),
    nsweep::Int=100,
    fps = 60,
    fontsize=20               
                   
    ) where N where T<:AbstractFloat

    fontsize_theme = Theme(fontsize = fontsize)
    set_theme!(fontsize_theme)
    
    meq = magn2d(β)
    progress = Progress(nsweep, 1)
    ising = Ising(I,h,β,x0)
    nspin = ising.N
    m = Observable([mean(ising.spin)])
    energia = Observable([energy(ising)/nspin])
    a = Observable(reshape(ising.spin,I))
    fig = Figure(resolution=(2048,1024))
    ax = Vector{Axis}(undef,3)
    ax[1] = Axis(fig[1:2,1])
    ax[2] = Axis(fig[1,2],title="magnetization",ylabel="M/N",titlesize=32,xlabelsize=32,ylabelsize=32)
    ax[3] = Axis(fig[2,2],title="energy",xlabel="mc sweep",ylabel="E/N",titlesize=32,xlabelsize=32,ylabelsize=32)
    pos = Vector{Makie.GridPosition}(undef,3)
    pos[1] = fig[1:2,1]
    pos[2] = fig[1,2]
    pos[3] = fig[2,2]
    xlims!(ax[2],1,nsweep)
    ylims!(ax[2],-1,1)
    xlims!(ax[3],1,nsweep)
    ylims!(ax[3],-2.01,0)
    heatmap!(pos[1],a)
    if meq != 0
        lines!(pos[2],1:nsweep, meq*ones(nsweep),color=:black,linewidth=3,label="+m")
        lines!(pos[2],1:nsweep,-meq*ones(nsweep),color=:black,linewidth=3,label="-m")
    else
        lines!(pos[2],1:nsweep, meq*ones(nsweep),color=:black,linewidth=3,label="m")
    end
    axislegend(ax[2])
    lines!(pos[2],m)
    lines!(pos[3],energia)
    lines!(pos[3],1:nsweep, derivative(free_energy2d,β)*ones(nsweep),color=:black,linewidth=3,label="e")
    display(fig)
    for _ in 1:nsweep
        onemcsweep!(ising)
        a[] = reshape(ising.spin,I)
        m[] = push!(m[],mean(ising.spin))
        energia[] = push!(energia[],energy(ising)/nspin)
        sleep(1.0/fps)
        #println(mean(ising.spin)," ",energy(ising)/nspin)
        next!(progress)
    end
end

function mytest() 
    fig = Figure()
    ax = [Axis(fig[1,i]) for i in 1:2]
    pos = [fig[1,i] for i in 1:2]
    xlims!(ax[1],0,1000)
    plot!(pos[1],rand(100))
    heatmap!(pos[2],rand(100,100))
    #xlims!(ax1,0,1000)
    display(fig)
end
