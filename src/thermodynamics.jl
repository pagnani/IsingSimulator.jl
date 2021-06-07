magn2d(x) = x > βc2d ? (1-(1/sinh(2x))^4)^(1/8) : 0.0
entropy(m) = -0.5*((1+m)*log(0.5*(1+m)) + (1-m)*log(0.5(1-m)))

# K = βJ L=βJ' J = J' for isotropic ising 2d

function fe_density2d(β,θ)
    k = 2sinh(2β)/cosh(2β)^2
    return log(1+sqrt(1-k^2 * cos(θ)^2))
end

function free_energy2d(β)
    f(θ) = fe_density2d(β,θ)
    res = QuadGK.quadgk(f,0,π)
    return -log(2)/2 - log(cosh(2β))-res[1]/(2π)
end

average_energy2d(b) = derivative(free_energy2d,b)
specific_heat2d(b) = -b^2 * derivative(average_energy2d,b)


function observables(binf,bsup;nstep=1024)
    vals = LinRange(binf,bsup,nstep)
    m = [magn2d(β) for β in vals]
    f = [free_energy2d(β) for β in vals]
    e = [average_energy2d(β) for β in vals]
    cv = [specific_heat2d(2β) for β in vals]
    s = [ vals[i]*e[i]-f[i] for i in eachindex(f)]
    return DataFrame(:b=>collect(vals),
                     :m=>m,
                     :f=>f./vals,
                     :e=>e,
                     :s=>s,                     
                     :cv=>cv)
end